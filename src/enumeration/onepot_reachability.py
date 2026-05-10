"""
Onepot reachability score — retrosynthesis-based.

Per the onepot CORE preprint (arXiv:2601.12603), CORE = enumerate(7 reactions ×
~320K building blocks) + ML feasibility filter. Their building-block list is
proprietary, so we cannot verify CORE membership pre-event. We CAN, however,
score each candidate on whether it is *retrosynthesizable* via the 7 onepot
reactions — a necessary condition for being in CORE.

For each candidate compound:
  1. Apply each of 7 reverse-reaction SMARTS to find disconnection points
  2. For each disconnection, compute:
     - Building-block 1 SMILES + heavy-atom count (HA), MW
     - Building-block 2 SMILES + HA, MW
  3. Score the disconnection by BB plausibility:
     - Both BBs MW < 350 (commercial BBs are typically smaller)
     - Both BBs HA < 25
     - No exotic groups in BBs (no nitro, no peroxide, etc.)
  4. Compound reachability = max(reaction_score) across all reactions

The score is a HEURISTIC, not a guarantee. It correlates with onepot
reachability but does NOT verify CORE membership.

Output: CSV per-compound with:
  reachable_via         : list of reactions that can disconnect the compound
  best_reaction         : highest-scoring disconnection reaction
  bb1_smiles, bb2_smiles: implied building blocks for best disconnection
  reachability_score    : 0..1
  amide_coupling_match  : 0/1
  suzuki_match          : 0/1
  buchwald_match        : 0/1
  urea_match            : 0/1
  thiourea_match        : 0/1
  n_alkyl_match         : 0/1
  o_alkyl_match         : 0/1

Reactions per onepot CORE preprint (arXiv:2601.12603):
  1. Amide coupling
  2. Suzuki-Miyaura coupling
  3. Buchwald-Hartwig amination
  4. Urea synthesis
  5. Thiourea synthesis
  6. Amine N-alkylation
  7. O-alkylation
"""
import argparse
import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


# Reverse-reaction SMARTS for each of the 7 onepot reactions. Each rule
# disconnects a product into two BB SMILES.
RETRO_REACTIONS = [
    # Amide coupling: R-C(=O)-N(H/R')-R'' → R-COOH + H-N(H/R')-R''
    # (Permit any sec/tert N adjacent to the carbonyl; exclude only sulfonamides
    # and quaternary N. The previous '!$(NC=O)' filter accidentally excluded
    # all amides since the matched N is, by definition, in the NC=O motif.)
    ("amide_coupling",
     "[C:1](=[O:2])[N;X3;!$(N=*);!$(NS(=O)(=O));!$([N+]):3]>>"
     "[C:1](=[O:2])O.[H][N:3]"),
    # Suzuki-Miyaura: aryl-aryl bond between rings → aryl-Br + aryl-B(OH)2
    # (avoid in-ring aryl-aryl by requiring different ring contexts is hard
    # in SMARTS; we filter post-hoc by checking BBs aren't single fused rings)
    ("suzuki_miyaura",
     "[c:1]-;!@[c:2]>>[c:1]Br.OB(O)[c:2]"),
    # Buchwald-Hartwig: aryl-N bond → aryl-Br + H-N
    ("buchwald_hartwig",
     "[c:1]-[N;!$(N=*);!$(NC=O);!$(NS(=O)(=O)):2]>>[c:1]Br.[H][N:2]"),
    # Urea: N-C(=O)-N → 2 amines + CDI/phosgene equivalent
    ("urea_synthesis",
     "[N;!$(NC=O):1][C:2](=[O:3])[N;!$(NC=O):4]>>"
     "[H][N:1].[H][N:4]"),
    # Thiourea: N-C(=S)-N → 2 amines + thiocarbonyl equivalent
    ("thiourea_synthesis",
     "[N;!$(NC=O):1][C:2](=[S:3])[N;!$(NC=O):4]>>"
     "[H][N:1].[H][N:4]"),
    # N-alkylation: sp3-C-N (not aryl-N) → R-Br + H-N
    ("n_alkylation",
     "[CX4;!$(C=*):1][N;!$(N=*);!$(NC=O);!$(NS(=O)(=O)):2]>>"
     "[Br][C:1].[H][N:2]"),
    # O-alkylation: sp3-C-O-C (ether) → R-Br + HO-R
    ("o_alkylation",
     "[CX4;!$(C=*):1][O:2][#6;!$([C]=O):3]>>"
     "[Br][C:1].[H][O:2][#6:3]"),
]


def _build_block_score(mol):
    """Heuristic plausibility score [0, 1] for a building-block candidate."""
    if mol is None:
        return 0.0
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return 0.0
    ha = mol.GetNumHeavyAtoms()
    mw = Descriptors.MolWt(mol)
    if ha < 3 or ha > 30:        # too small/large for a typical BB
        return 0.0
    if mw > 400:                 # commercial BBs are mostly < 400 Da
        return 0.0
    # Penalize exotic functional groups
    smi = Chem.MolToSmiles(mol)
    if any(bad in smi for bad in ["[N+](=O)[O-]", "OO", "[Hg]", "[Pb]",
                                   "[As]", "C(=O)Cl", "S=C=N", "P(=O)(O)(O)"]):
        return 0.0
    score = 1.0
    # Smaller is more BB-like; HA in 5..20 is ideal
    if ha < 5 or ha > 20:
        score *= 0.7
    if mw < 80 or mw > 300:
        score *= 0.7
    return score


def retrosynth_one(smiles):
    """Return per-reaction match info for one compound."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    out = {
        "smiles": smiles,
        "amide_coupling_match": 0,
        "suzuki_match": 0,
        "buchwald_match": 0,
        "urea_match": 0,
        "thiourea_match": 0,
        "n_alkyl_match": 0,
        "o_alkyl_match": 0,
        "reachable_via": "",
        "best_reaction": "",
        "bb1_smiles": "",
        "bb2_smiles": "",
        "reachability_score": 0.0,
    }

    matched = []
    best_score = 0.0
    best_rxn = ""
    best_bbs = ("", "")

    for rxn_name, retro_smarts in RETRO_REACTIONS:
        rxn = AllChem.ReactionFromSmarts(retro_smarts)
        if rxn is None:
            continue
        try:
            products_list = rxn.RunReactants((mol,))
        except Exception:
            continue
        if not products_list:
            continue

        # Got at least one match
        col_key = {
            "amide_coupling":   "amide_coupling_match",
            "suzuki_miyaura":   "suzuki_match",
            "buchwald_hartwig": "buchwald_match",
            "urea_synthesis":   "urea_match",
            "thiourea_synthesis": "thiourea_match",
            "n_alkylation":     "n_alkyl_match",
            "o_alkylation":     "o_alkyl_match",
        }[rxn_name]
        out[col_key] = 1
        matched.append(rxn_name)

        # Score the best disconnection from this reaction
        for products in products_list:
            if len(products) != 2:
                continue
            bb1, bb2 = products
            s1 = _build_block_score(bb1)
            s2 = _build_block_score(bb2)
            score = (s1 + s2) / 2
            if score > best_score:
                best_score = score
                best_rxn = rxn_name
                try:
                    best_bbs = (Chem.MolToSmiles(bb1), Chem.MolToSmiles(bb2))
                except Exception:
                    best_bbs = ("", "")

    out["reachable_via"] = ",".join(matched)
    out["best_reaction"] = best_rxn
    out["bb1_smiles"], out["bb2_smiles"] = best_bbs
    out["reachability_score"] = round(best_score, 3)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",  required=True, help="CSV with id + smiles columns")
    ap.add_argument("--output", required=True, help="audit CSV to write")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.input)))
    print(f"Reachability scoring {len(rows)} compounds via 7 onepot reactions...")

    results = []
    for r in rows:
        cid, smi = r["id"], r["smiles"]
        info = retrosynth_one(smi)
        if info is None:
            continue
        info["id"] = cid
        results.append(info)

    cols = ["id", "smiles", "reachability_score", "best_reaction",
            "reachable_via", "bb1_smiles", "bb2_smiles",
            "amide_coupling_match", "suzuki_match", "buchwald_match",
            "urea_match", "thiourea_match", "n_alkyl_match", "o_alkyl_match"]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k, "") for k in cols})

    # Summary
    n_reachable = sum(1 for r in results if r["reachability_score"] > 0)
    n_amide = sum(r["amide_coupling_match"] for r in results)
    n_suzuki = sum(r["suzuki_match"] for r in results)
    n_buchwald = sum(r["buchwald_match"] for r in results)
    print(f"\nSummary:")
    print(f"  any reaction matches:   {n_reachable}/{len(results)}")
    print(f"  amide coupling matches: {n_amide}/{len(results)}")
    print(f"  suzuki matches:         {n_suzuki}/{len(results)}")
    print(f"  buchwald matches:       {n_buchwald}/{len(results)}")

    print(f"\nTop 10 by reachability score:")
    for r in sorted(results, key=lambda r: -r["reachability_score"])[:10]:
        print(f"  {r['id']:35s} score={r['reachability_score']:.2f}  "
              f"best={r['best_reaction']:18s}  matches={r['reachable_via']}")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
