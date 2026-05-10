"""
Analog enumeration for TBXT priority scaffolds.

Three approaches combined:
  1. BRICS decomposition + recombination — fragment-shuffle the elaborated scaffolds
     using RDKit's BRICS module. Best for large compounds like Z795991852.
  2. Reaction-SMARTS bioisostere swaps — apply chemistry-aware substitutions
     (carboxylate → amide/sulfonamide/tetrazole; phenyl → pyridyl; etc.). Best for
     small fragments where we want "same pharmacophore, different atom".
  3. R-group attachment — grow small fragments by attaching curated R-groups at
     aromatic/aliphatic positions.

Filters applied to all candidates:
  • Valid SMILES, sanitized
  • Relaxed lead-like rule: Chordoma hard + HA ≤ 35 + rings ≤ 6 + fused ≤ 2 + no PAINS
  • Tanimoto < 0.85 to ALL Naar (avoid duplication)
  • Tanimoto > 0.40 to parent (preserve binding signal)
  • Onepot 7-reaction retro check is heuristic only — assume the on-day Onepot lookup
    is the authoritative filter; we just bias toward common synthesizable motifs.

Outputs:
  data/analogs/<scaffold_id>_candidates.csv — per-scaffold pool (capped ~200)
  data/analogs/all_candidates.csv          — combined master pool
  data/analogs/enumeration_log.txt         — provenance + counts
"""
import csv
import logging
import sys
import time
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import AllChem, BRICS, FilterCatalog, Lipinski, Descriptors
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("analogs")

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = DATA / "analogs"
OUT.mkdir(exist_ok=True)

# Priority scaffolds — see FINDINGS.md / PROGRESS.md for selection rationale
SCAFFOLDS = {
    "Z795991852": ("Cn1c(=O)c2ccccc2n2c(COC(=O)c3cccc(NC(=O)C4Cc5ccccc5O4)c3)nnc12",
                   "BRICS",
                   "CF Labs SPR hit, site F, 10 µM. Most novel-friendly (T 0.27 to nearest Naar)."),
    "FM001580":   ("O=C(O)c1ccccc1OC(F)(F)F",
                   "GROW",
                   "Site F TEP fragment, 2-(trifluoromethoxy)benzoic acid."),
    "FM001452":   ("Nc1cccc(OCc2ccccc2)c1",
                   "GROW",
                   "Site F TEP fragment, 3-(benzyloxy)aniline."),
    "FM002150":   ("O=C(O)COCc1ccccc1",
                   "GROW",
                   "Site F TEP fragment, (benzyloxy)acetic acid."),
}

PER_SCAFFOLD_CAP = 200
MIN_TANIMOTO_TO_PARENT = 0.40
MAX_TANIMOTO_TO_NAAR = 0.85

# Bioisostere SMARTS rules: (name, reaction_smarts, comment)
# Reaction syntax: react→prod where match-and-replace pattern is used.
# We use Chem.AllChem.ReactionFromSmarts: (reactant) >> (product)
BIOISOSTERE_RULES = [
    # Carboxylic acid bioisosteres
    ("acid_to_methylamide",     "[c:1][C:2](=[O:3])[OH]>>[c:1][C:2](=[O:3])N(C)C"),
    ("acid_to_primary_amide",   "[c:1][C:2](=[O:3])[OH]>>[c:1][C:2](=[O:3])N"),
    ("acid_to_tetrazole",       "[c:1][C:2](=O)[OH]>>[c:1]c1nnn[nH]1"),
    ("acid_to_sulfonamide",     "[c:1][C:2](=O)[OH]>>[c:1]S(=O)(=O)N"),
    ("acid_to_acylsulfonamide", "[c:1][C:2](=O)[OH]>>[c:1]C(=O)NS(=O)(=O)C"),
    # Aromatic ring scaffold-hops on benzene with a key substituent
    ("phenyl_to_pyridyl_3",     "[cH:1]1[cH:2][cH:3][c:4]([*:5])[cH:6][cH:7]1>>[n:1]1[cH:2][cH:3][c:4]([*:5])[cH:6][cH:7]1"),
    ("phenyl_to_pyridyl_4",     "[cH:1]1[cH:2][cH:3][c:4]([*:5])[cH:6][cH:7]1>>[cH:1]1[n:2][cH:3][c:4]([*:5])[cH:6][cH:7]1"),
    # Ether → thioether / amide
    ("oxy_methylene_to_amide",  "[c:1][O:2][CH2:3][c:4]>>[c:1][NH:2][C:3](=O)[c:4]"),
    ("oxy_methylene_to_thio",   "[c:1][O:2][CH2:3][c:4]>>[c:1][S:2][CH2:3][c:4]"),
    # Trifluoromethoxy variations
    ("ocf3_to_cf3",             "[c:1][O:2][C:3]([F:4])([F:5])[F:6]>>[c:1][C:3]([F:4])([F:5])[F:6]"),
    ("ocf3_to_oet",             "[c:1][O:2][C:3]([F:4])([F:5])[F:6]>>[c:1][O:2]CC"),
    ("ocf3_to_omec",            "[c:1][O:2][C:3]([F:4])([F:5])[F:6]>>[c:1][O:2]C(F)(F)F"),  # noop variant for control
    # Aniline variations
    ("aniline_to_anilide_ac",   "[c:1][NH2:2]>>[c:1][NH:2]C(=O)C"),
    ("aniline_to_anilide_pr",   "[c:1][NH2:2]>>[c:1][NH:2]C(=O)CC"),
    ("aniline_to_methylamine",  "[c:1][NH2:2]>>[c:1][NH:2]C"),
    ("aniline_to_dimethyl",     "[c:1][NH2:2]>>[c:1]N(C)C"),
    ("aniline_to_sulfonamide",  "[c:1][NH2:2]>>[c:1][NH:2]S(=O)(=O)C"),
    ("aniline_to_urea_methyl",  "[c:1][NH2:2]>>[c:1][NH:2]C(=O)NC"),
]

# Small R-groups for ring substitution / fragment growing.
# For aromatic CH positions and aniline NH2, we'll attach these via reaction SMARTS.
GROW_RGROUPS_AROMATIC_H = [
    "F", "Cl", "C(=O)O", "OC", "OCC", "C(F)(F)F",
    "S(=O)(=O)C", "S(=O)(=O)N", "C(=O)N", "CN(C)C",
    "C", "N", "C(=O)NC", "C#N", "OCCO", "CO",
    "CC", "CCC", "C(C)C", "CN", "N(C)C",
    "c1ccncc1", "c1ccccn1", "c1ccoc1", "c1ccsc1",  # heteroaryl
    "C1CCCCC1", "C1CCNCC1", "C1CCOCC1", "C1CCN(C)CC1",  # saturated rings
]


def is_valid_canonical(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    try: Chem.SanitizeMol(m)
    except Exception: return None
    return Chem.MolToSmiles(m)


def descriptors_for(mol):
    return {
        "ha": mol.GetNumHeavyAtoms(),
        "mw": round(Descriptors.MolWt(mol), 1),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "rings": mol.GetRingInfo().NumRings(),
        "fused_rings": _count_fused(mol),
        "rotb": Lipinski.NumRotatableBonds(mol),
        "tpsa": round(Descriptors.TPSA(mol), 1),
    }


def _count_fused(mol):
    rings = mol.GetRingInfo().AtomRings()
    fused = 0
    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            if len(set(rings[i]) & set(rings[j])) >= 2:
                fused += 1
                break
    return fused


# ── PAINS catalog ─────────────────────────────────────────────────────────────
_pains_params = FilterCatalog.FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
_pains_catalog = FilterCatalog.FilterCatalog(_pains_params)


def is_pains(mol):
    return _pains_catalog.HasMatch(mol)


def passes_filters(mol):
    """Relaxed lead-like rule from FINDINGS.md (calibrated to accept CF Labs hits)."""
    d = descriptors_for(mol)
    # Chordoma hard rule
    if not (d["logp"] <= 6 and d["hbd"] <= 6 and d["hba"] <= 12 and d["mw"] <= 600):
        return False, d, "chordoma_hard"
    # Relaxed lead-like
    if not (d["ha"] <= 35 and d["rings"] <= 6 and d["fused_rings"] <= 2):
        return False, d, "lead_like_relaxed"
    if is_pains(mol):
        return False, d, "pains"
    return True, d, "ok"


# ── Tanimoto utilities ─────────────────────────────────────────────────────────
mfg = GetMorganGenerator(radius=2, fpSize=2048)


def fp(mol):
    return mfg.GetFingerprint(mol)


def load_naar_fingerprints():
    """Return list of (id, fp) for all Naar-derived compounds in the inventory."""
    log.info("Loading Naar fingerprints from prior_art_canonical.csv...")
    naar_fps = []
    with open(DATA / "prior_art_canonical.csv") as f:
        for r in csv.DictReader(f):
            if r["source"] in ("naar_sheet", "naar_zenodo") and r["smiles"]:
                m = Chem.MolFromSmiles(r["smiles"])
                if m is not None:
                    naar_fps.append((r["id"], fp(m)))
    log.info(f"  Loaded {len(naar_fps)} Naar fingerprints")
    return naar_fps


def max_tanimoto_to_naar(mol, naar_fps):
    f = fp(mol)
    sims = DataStructs.BulkTanimotoSimilarity(f, [x[1] for x in naar_fps])
    bi = max(range(len(sims)), key=lambda i: sims[i])
    return float(sims[bi]), naar_fps[bi][0]


# ── Approach 1: BRICS recombination ────────────────────────────────────────────
def enumerate_brics(parent_mol, max_attempts=20000, max_yield=2000):
    """Decompose parent into BRICS fragments, then build new molecules.

    Use BRICSBuild with the parent's own fragment pool (auto-shuffles them)
    and limit the number generated.
    """
    frags = BRICS.BRICSDecompose(parent_mol, returnMols=True)
    log.info(f"  BRICS produced {len(frags)} fragments")
    if len(frags) < 2:
        return []
    # BRICSBuild is a generator
    out = []
    seen = set()
    try:
        for i, m in enumerate(BRICS.BRICSBuild(frags, scrambleReagents=True, maxDepth=4)):
            if i >= max_attempts: break
            if m is None: continue
            try: Chem.SanitizeMol(m)
            except Exception: continue
            smi = Chem.MolToSmiles(m)
            if smi in seen: continue
            seen.add(smi)
            out.append(smi)
            if len(out) >= max_yield: break
    except Exception as e:
        log.warning(f"  BRICSBuild ended early: {e}")
    return out


# ── Approach 2: Reaction-SMARTS bioisosteres ───────────────────────────────────
def apply_bioisostere_rules(parent_smi):
    out = []
    for name, smarts in BIOISOSTERE_RULES:
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            if rxn is None: continue
            m = Chem.MolFromSmiles(parent_smi)
            if m is None: continue
            products = rxn.RunReactants((m,))
            for prod_set in products:
                for p in prod_set:
                    try:
                        Chem.SanitizeMol(p)
                        out.append((Chem.MolToSmiles(p), name))
                    except Exception:
                        continue
        except Exception:
            continue
    # Dedup
    uniq, seen = [], set()
    for smi, src in out:
        if smi in seen: continue
        seen.add(smi); uniq.append((smi, src))
    return uniq


# ── Approach 3: Aromatic CH substitution / fragment growing ────────────────────
def grow_at_aromatic_ch(parent_smi, rgroups=GROW_RGROUPS_AROMATIC_H):
    """For each aromatic CH on the parent, replace H with each R-group via SMARTS."""
    out = []
    parent = Chem.MolFromSmiles(parent_smi)
    if parent is None: return out
    # Build a SMARTS that matches an aromatic CH atom; we'll iterate per-position.
    matches = parent.GetSubstructMatches(Chem.MolFromSmarts("[cH]"))
    for atom_idx_tuple in matches:
        ai = atom_idx_tuple[0]
        for rg in rgroups:
            # Use editable mol approach: replace H with the R-group fragment
            try:
                rw = Chem.RWMol(parent)
                rw.GetAtomWithIdx(ai).SetNumExplicitHs(0)
                rw.GetAtomWithIdx(ai).SetNoImplicit(True)
                rg_mol = Chem.MolFromSmiles(rg)
                if rg_mol is None: continue
                # Combine: append rg_mol then bond the parent's atom to rg's first atom
                merged = Chem.CombineMols(rw, rg_mol)
                ed = Chem.RWMol(merged)
                # rg_mol atom indices start at parent.NumAtoms()
                rg_first = parent.GetNumAtoms()
                ed.AddBond(ai, rg_first, Chem.BondType.SINGLE)
                # Sanitize
                m_out = ed.GetMol()
                Chem.SanitizeMol(m_out)
                smi = Chem.MolToSmiles(m_out)
                out.append((smi, f"grow_aromCH_{rg}"))
            except Exception:
                continue
    # Also try amine N-H attachments for anilines etc.
    nh_matches = parent.GetSubstructMatches(Chem.MolFromSmarts("[c][NH2]"))
    for c_idx, n_idx in nh_matches:
        for rg in ["C(=O)C", "C(=O)CC", "C(=O)c1ccccc1", "S(=O)(=O)C", "C(=O)N"]:
            try:
                rw = Chem.RWMol(parent)
                # convert NH2 -> NH-rg by changing the N atom and adding the rg
                rw.GetAtomWithIdx(n_idx).SetNumExplicitHs(1)
                rg_mol = Chem.MolFromSmiles(rg)
                if rg_mol is None: continue
                merged = Chem.CombineMols(rw, rg_mol)
                ed = Chem.RWMol(merged)
                ed.AddBond(n_idx, parent.GetNumAtoms(), Chem.BondType.SINGLE)
                m_out = ed.GetMol()
                Chem.SanitizeMol(m_out)
                smi = Chem.MolToSmiles(m_out)
                out.append((smi, f"grow_aniline_{rg[:10]}"))
            except Exception:
                continue
    # Dedup
    uniq, seen = [], set()
    for smi, src in out:
        if smi in seen: continue
        seen.add(smi); uniq.append((smi, src))
    return uniq


# ── Main pipeline ──────────────────────────────────────────────────────────────
def enumerate_for_scaffold(scaffold_id, smiles, method, naar_fps):
    log.info(f"\n=== Enumerating analogs for {scaffold_id} (method={method}) ===")
    parent_mol = Chem.MolFromSmiles(smiles)
    parent_smi = Chem.MolToSmiles(parent_mol)
    parent_fp = fp(parent_mol)
    parent_d = descriptors_for(parent_mol)
    log.info(f"  parent: HA={parent_d['ha']}, MW={parent_d['mw']}, rings={parent_d['rings']}")

    raw_candidates = []   # (smi, source_method)

    if method in ("BRICS", "BOTH"):
        t0 = time.time()
        brics_smis = enumerate_brics(parent_mol)
        log.info(f"  BRICS: {len(brics_smis)} unique molecules in {time.time()-t0:.1f}s")
        raw_candidates.extend([(s, "brics") for s in brics_smis])

    if method in ("GROW", "BOTH"):
        t0 = time.time()
        bio = apply_bioisostere_rules(smiles)
        log.info(f"  Bioisostere rules: {len(bio)} unique products in {time.time()-t0:.1f}s")
        raw_candidates.extend(bio)

        t0 = time.time()
        grown = grow_at_aromatic_ch(smiles)
        log.info(f"  Aromatic-CH growing: {len(grown)} unique products in {time.time()-t0:.1f}s")
        raw_candidates.extend(grown)

    # Filter
    log.info(f"  Total raw: {len(raw_candidates)}; applying filters...")
    survivors = []
    seen_canonical = set()
    n_invalid = n_filter = n_naardup = n_parentlow = 0
    for smi, src in raw_candidates:
        c = is_valid_canonical(smi)
        if c is None:
            n_invalid += 1; continue
        if c in seen_canonical: continue
        seen_canonical.add(c)
        if c == parent_smi: continue
        m = Chem.MolFromSmiles(c)
        ok, d, why = passes_filters(m)
        if not ok:
            n_filter += 1; continue
        # Tanimoto checks
        tp = float(DataStructs.TanimotoSimilarity(parent_fp, fp(m)))
        if tp < MIN_TANIMOTO_TO_PARENT:
            n_parentlow += 1; continue
        tn, tn_id = max_tanimoto_to_naar(m, naar_fps)
        if tn > MAX_TANIMOTO_TO_NAAR:
            n_naardup += 1; continue
        survivors.append({
            "id": f"{scaffold_id}_analog_{len(survivors)+1:04d}",
            "smiles": c,
            "parent_id": scaffold_id,
            "parent_smiles": parent_smi,
            "method": src,
            **d,
            "tanimoto_to_parent": round(tp, 3),
            "max_tanimoto_to_naar": round(tn, 3),
            "naar_neighbor": tn_id,
        })

    log.info(f"  Survivors: {len(survivors)} | "
             f"filtered out: {n_invalid} invalid, {n_filter} property-fail, "
             f"{n_naardup} dup-of-naar, {n_parentlow} too-distant-from-parent")

    # Cap with diverse sampling: sort by max_tanimoto_to_naar ascending (most novel first),
    # then take the top PER_SCAFFOLD_CAP
    survivors.sort(key=lambda r: (r["max_tanimoto_to_naar"], -r["tanimoto_to_parent"]))
    capped = survivors[:PER_SCAFFOLD_CAP]
    log.info(f"  Capped to {len(capped)} after sorting by novelty + parent-similarity")
    return capped, {
        "raw": len(raw_candidates),
        "invalid": n_invalid,
        "property_fail": n_filter,
        "dup_of_naar": n_naardup,
        "too_distant_from_parent": n_parentlow,
        "survivors": len(survivors),
        "kept": len(capped),
    }


def main():
    naar_fps = load_naar_fingerprints()

    all_candidates = []
    log_lines = ["# Analog Enumeration Log\n"]
    log_lines.append(f"Per-scaffold cap: {PER_SCAFFOLD_CAP}")
    log_lines.append(f"Filters: passes_relaxed_leadlike, T<{MAX_TANIMOTO_TO_NAAR} to Naar, "
                     f"T>{MIN_TANIMOTO_TO_PARENT} to parent\n")

    for sid, (smi, method, comment) in SCAFFOLDS.items():
        candidates, stats = enumerate_for_scaffold(sid, smi, method, naar_fps)

        # Write per-scaffold CSV
        if candidates:
            cols = list(candidates[0].keys())
            with open(OUT / f"{sid}_candidates.csv", "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=cols)
                w.writeheader()
                w.writerows(candidates)
        log_lines.append(f"## {sid} ({method}) — {comment}")
        log_lines.append(f"- Parent SMILES: `{smi}`")
        for k, v in stats.items():
            log_lines.append(f"- {k}: {v}")
        log_lines.append("")
        all_candidates.extend(candidates)

    # Write combined master file
    if all_candidates:
        cols = list(all_candidates[0].keys())
        with open(OUT / "all_candidates.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(all_candidates)
        log.info(f"\nWrote {len(all_candidates)} total candidates to all_candidates.csv")

    log_lines.append(f"\n## TOTAL: {len(all_candidates)} candidates across all scaffolds")
    (OUT / "enumeration_log.md").write_text("\n".join(log_lines))
    log.info(f"Wrote enumeration_log.md")


if __name__ == "__main__":
    main()
