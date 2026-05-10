"""Apply the seven-criterion strict gate to a scored candidate set.

Implements the rules documented in ``docs/filter_chain.md``. Every
candidate gets per-criterion pass/fail flag columns; the strict-pass
subset is the rows where every flag is `'✓'`.

This module does not own the *evaluation* of each criterion - those
live in their own modules under ``src/filters/`` (e.g. PAINS lives
in ``filters/pains_and_motifs.py``, onepot membership in
``filters/onepot_membership.py``). This module is the orchestrator
that runs every check and writes the per-criterion flag columns
needed by ``ranking/tier_classify.py``.

Usage::

    python -m src.ranking.strict_gate \\
        --input  data/scored/all_candidates_with_scores.csv \\
        --output data/scored/all_candidates_with_flags.csv \\
        --naar   data/naar/naar_spr_kd.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, RDConfig, rdMolDescriptors

# RDKit's bundled FilterCatalog ships PAINS A/B/C
from rdkit.Chem import FilterCatalog


PASS = "✓"
FAIL = "✗"


# --- C2 reactive / covalent SMARTS -----------------------------------
COVALENT_SMARTS = {
    "michael_acceptor": "[CX3]=[CX3]C(=O)",
    "alpha_beta_unsat_carbonyl": "C=CC(=O)",
    "boron_any": "[B]",
    "vinyl_sulfone": "C=CS(=O)(=O)",
    "aldehyde": "[CX3H1](=O)[#6]",
    "alkyl_halide": "[CX4][F,Cl,Br,I]",
    "epoxide": "C1CO1",
    "aziridine": "C1CN1",
    "isocyanate": "N=C=O",
}

# --- C5 forbidden motifs --------------------------------------------
FORBIDDEN_SMARTS = {
    "acid_halide": "C(=O)[F,Cl,Br,I]",
    "diazo": "[N+]#N",
    "aliphatic_imine": "[CX3]=[NX2;!$(N=*-*=*)]",
    "long_alkyl_chain": "CCCCCCCCC",   # > C8
}


def _smarts_match(mol: Chem.Mol, smarts: str) -> bool:
    pat = Chem.MolFromSmarts(smarts)
    return pat is not None and mol.HasSubstructMatch(pat)


def _esol_log_s(mol: Chem.Mol) -> float:
    """ESOL (Delaney 2004) log(S) - fast solubility predictor."""
    logp = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    arom = sum(1 for ring in mol.GetRingInfo().AtomRings()
               if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    n_atoms = mol.GetNumHeavyAtoms()
    aromatic_frac = arom / max(n_atoms, 1)
    return 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotb - 0.74 * aromatic_frac


def _ring_metrics(mol: Chem.Mol) -> tuple[int, int]:
    """(n_rings, n_fused_rings) where 'fused' = sharing >= 2 atoms with another ring."""
    info = mol.GetRingInfo()
    rings = info.AtomRings()
    n = len(rings)
    n_fused = 0
    for i in range(n):
        for j in range(i + 1, n):
            shared = set(rings[i]) & set(rings[j])
            if len(shared) >= 2:
                n_fused += 1
    return n, n_fused


def _evaluate_row(smiles: str, naar_max_tanimoto: float | None,
                  pains_catalog: FilterCatalog.FilterCatalog) -> dict:
    """Return a dict of per-criterion flag values for one compound."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {c: FAIL for c in (
            "C2_non_covalent",
            "C3a_MW_le_600", "C3b_LogP_le_6", "C3c_HBD_le_6", "C3d_HBA_le_12",
            "C4a_HA_10_30", "C4b_HBD_HBA_le_11", "C4c_LogP_lt_5",
            "C4d_lt_5_rings", "C4e_le_2_fused",
            "C5_no_PAINS", "C5b_no_forbidden",
            "C6_naar_tanimoto_lt_085", "C7_soluble_logS_gt_neg5",
        )}

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    ha = mol.GetNumHeavyAtoms()
    n_rings, n_fused = _ring_metrics(mol)
    log_s = _esol_log_s(mol)

    cov_hits = [k for k, sm in COVALENT_SMARTS.items() if _smarts_match(mol, sm)]
    forb_hits = [k for k, sm in FORBIDDEN_SMARTS.items() if _smarts_match(mol, sm)]
    pains_hit = pains_catalog.HasMatch(mol)

    flags = {
        "C2_non_covalent": PASS if not cov_hits else FAIL,
        "C3a_MW_le_600":  PASS if mw <= 600 else FAIL,
        "C3b_LogP_le_6":  PASS if logp <= 6 else FAIL,
        "C3c_HBD_le_6":   PASS if hbd <= 6 else FAIL,
        "C3d_HBA_le_12":  PASS if hba <= 12 else FAIL,
        "C4a_HA_10_30":   PASS if 10 <= ha <= 30 else FAIL,
        "C4b_HBD_HBA_le_11": PASS if (hbd + hba) <= 11 else FAIL,
        "C4c_LogP_lt_5":  PASS if logp < 5 else FAIL,
        "C4d_lt_5_rings": PASS if n_rings < 5 else FAIL,
        "C4e_le_2_fused": PASS if n_fused <= 2 else FAIL,
        "C5_no_PAINS":    PASS if not pains_hit else FAIL,
        "C5b_no_forbidden": PASS if not forb_hits else FAIL,
        "C6_naar_tanimoto_lt_085": (
            PASS if (naar_max_tanimoto is None or naar_max_tanimoto < 0.85) else FAIL
        ),
        "C7_soluble_logS_gt_neg5": PASS if log_s > -5 else FAIL,
        "MW_Da": round(mw, 2),
        "LogP": round(logp, 2),
        "HBD": hbd, "HBA": hba, "HA": ha,
        "rings": n_rings, "fused_rings": n_fused,
        "esol_logS": round(log_s, 2),
        "forbidden_motifs": ",".join(forb_hits) if forb_hits else "-",
    }
    return flags


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True,
                    help="Scored candidates CSV with `id` + `smiles` columns "
                         "and (optionally) the C1_onepot_100 flag pre-attached.")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--naar", type=Path, required=False,
                    help="(Optional) Naar SPR Kd CSV with a `smiles` column for "
                         "Tanimoto novelty check (C6).")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    if "smiles" not in df.columns:
        sys.exit("ERROR: input CSV must contain a 'smiles' column.")

    # --- C6 prep: precompute Naar fingerprints if the file is given ---
    from rdkit.Chem import AllChem, DataStructs
    naar_fps = []
    if args.naar and args.naar.exists():
        ndf = pd.read_csv(args.naar)
        for smi in ndf["smiles"].dropna():
            m = Chem.MolFromSmiles(smi)
            if m:
                naar_fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))

    def _max_naar_tanimoto(smi: str) -> float | None:
        if not naar_fps:
            return None
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
        return max(DataStructs.TanimotoSimilarity(fp, ref) for ref in naar_fps)

    # --- PAINS catalog (RDKit-bundled A/B/C) ---
    params = FilterCatalog.FilterCatalogParams()
    for cat in (FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A,
                FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B,
                FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C):
        params.AddCatalog(cat)
    pains_catalog = FilterCatalog.FilterCatalog(params)

    # --- evaluate row by row ---
    rows = []
    for _, r in df.iterrows():
        flags = _evaluate_row(r["smiles"], _max_naar_tanimoto(r["smiles"]), pains_catalog)
        rows.append({**r.to_dict(), **flags,
                     "naar_max_tanimoto": _max_naar_tanimoto(r["smiles"])})

    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    n_strict = (out[[c for c in out.columns if c.startswith(("C1_", "C2_", "C3", "C4", "C5", "C6_", "C7_"))]]
                .eq(PASS).all(axis=1).sum())
    print(f"  Wrote {len(out)} rows -> {args.output}")
    print(f"  Strict-pass (all 7 criteria PASS): {n_strict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
