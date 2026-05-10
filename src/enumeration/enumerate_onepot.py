"""
Enumerate compound products from a building-block list using onepot CORE's
7 reactions (per arXiv:2601.12603). Forward direction.

Reactions (RDKit reaction SMARTS):
  1. Amide coupling:   acid + amine -> amide
  2. Suzuki-Miyaura:   aryl-halide + aryl-boronic acid -> biaryl
  3. Buchwald-Hartwig: aryl-halide + amine -> aryl-amine
  4. Urea synthesis:   2 amines + CDI/triphosgene -> urea
  5. Thiourea:         2 amines + thiocarbonyl -> thiourea
  6. N-alkylation:     amine + alkyl halide -> alkylamine
  7. O-alkylation:     phenol/alcohol + alkyl halide -> ether

Output: SMILES CSV of unique products (deduplicated by canonical SMILES).
"""
import argparse
import csv
import itertools
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


# Forward reaction SMARTS — each maps two BB SMARTS to a single product.
REACTION_TEMPLATES = [
    # 1. Amide coupling
    ("amide_coupling",
     "[C:1](=[O:2])[OH].[N;X3;H1,H2:3]>>[C:1](=[O:2])[N:3]"),
    # 2. Suzuki-Miyaura  (aryl-halide + aryl-boronic acid -> biaryl)
    ("suzuki_miyaura",
     "[c:1][Br,Cl,I].[c:2]B(O)O>>[c:1]-[c:2]"),
    # 3. Buchwald-Hartwig
    ("buchwald_hartwig",
     "[c:1][Br,Cl,I].[N;X3;H1,H2:2]>>[c:1][N:2]"),
    # 4. Urea synthesis (two primary or secondary amines, formal CDI)
    ("urea_synthesis",
     "[N;X3;H1,H2:1].[N;X3;H1,H2:2]>>[N:1]C(=O)[N:2]"),
    # 5. Thiourea synthesis
    ("thiourea_synthesis",
     "[N;X3;H1,H2:1].[N;X3;H1,H2:2]>>[N:1]C(=S)[N:2]"),
    # 6. N-alkylation (amine + alkyl halide; require sp3 carbon)
    ("n_alkylation",
     "[N;X3;H1,H2:1].[CX4;!$(C=*):2][Br,Cl,I]>>[N:1][C:2]"),
    # 7. O-alkylation (alcohol/phenol + alkyl halide)
    ("o_alkylation",
     "[O;X2;H1:1][#6:2].[CX4;!$(C=*):3][Br,Cl,I]>>[O:1]([#6:2])[C:3]"),
]


def _matches(mol, patt_smarts):
    patt = Chem.MolFromSmarts(patt_smarts)
    return patt is not None and mol.HasSubstructMatch(patt)


def _has_acid(m):     return _matches(m, "[C](=O)[OH]")
def _has_amine(m):    return _matches(m, "[N;X3;H1,H2;!$(NC=O);!$(NS(=O)(=O))]")
def _has_aryl_x(m):   return _matches(m, "[c][Br,Cl,I]")
def _has_aryl_b(m):   return _matches(m, "[c]B(O)O")
def _has_alkyl_x(m):  return _matches(m, "[CX4;!$(C=*)][Br,Cl,I]")
def _has_oh(m):       return _matches(m, "[O;X2;H1][#6;!$([C]=O)]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bbs", required=True, help="CSV with smiles,id columns")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-products", type=int, default=50000)
    args = ap.parse_args()

    bb_rows = list(csv.DictReader(open(args.bbs)))
    bbs = []
    for r in bb_rows:
        smi = r.get("smiles") or r.get("SMILES") or ""
        if not smi:
            continue
        m = Chem.MolFromSmiles(smi)
        if m is None or m.GetNumHeavyAtoms() < 3 or m.GetNumHeavyAtoms() > 25:
            continue
        bbs.append({
            "id": r.get("id") or r.get("ID") or smi,
            "smiles": Chem.MolToSmiles(m),
            "mol": m,
            "is_acid":     _has_acid(m),
            "is_amine":    _has_amine(m),
            "is_aryl_x":   _has_aryl_x(m),
            "is_aryl_b":   _has_aryl_b(m),
            "is_alkyl_x":  _has_alkyl_x(m),
            "is_oh":       _has_oh(m),
        })

    print(f"Loaded {len(bbs)} usable BBs from {args.bbs}")

    # Quick functional-group inventory to avoid pointless enumerations
    n_acid    = sum(1 for b in bbs if b["is_acid"])
    n_amine   = sum(1 for b in bbs if b["is_amine"])
    n_aryl_x  = sum(1 for b in bbs if b["is_aryl_x"])
    n_aryl_b  = sum(1 for b in bbs if b["is_aryl_b"])
    n_alkyl_x = sum(1 for b in bbs if b["is_alkyl_x"])
    n_oh      = sum(1 for b in bbs if b["is_oh"])
    print(f"  acids={n_acid} amines={n_amine} aryl-X={n_aryl_x} aryl-B={n_aryl_b} "
          f"alkyl-X={n_alkyl_x} -OH={n_oh}")

    seen = set()
    products = []

    def _add(smi, rxn_name, bb1_id, bb2_id):
        if not smi or len(seen) >= args.max_products:
            return
        m = Chem.MolFromSmiles(smi)
        if m is None: return
        try:
            Chem.SanitizeMol(m)
        except Exception:
            return
        canonical = Chem.MolToSmiles(m)
        if canonical in seen:
            return
        seen.add(canonical)
        products.append({
            "id": f"opv1_{len(seen):06d}",
            "smiles": canonical,
            "reaction": rxn_name,
            "bb1": bb1_id,
            "bb2": bb2_id,
        })

    # Run each reaction over the relevant BB pairs
    for rxn_name, smarts in REACTION_TEMPLATES:
        if len(seen) >= args.max_products: break
        rxn = AllChem.ReactionFromSmarts(smarts)
        if rxn is None:
            print(f"  skipping {rxn_name}: invalid SMARTS")
            continue

        # Restrict the BB pairs to plausible ones (massive speed-up)
        if rxn_name == "amide_coupling":
            pairs = ((a, b) for a in bbs for b in bbs if a["is_acid"] and b["is_amine"])
        elif rxn_name == "suzuki_miyaura":
            pairs = ((a, b) for a in bbs for b in bbs if a["is_aryl_x"] and b["is_aryl_b"])
        elif rxn_name == "buchwald_hartwig":
            pairs = ((a, b) for a in bbs for b in bbs if a["is_aryl_x"] and b["is_amine"])
        elif rxn_name in ("urea_synthesis", "thiourea_synthesis"):
            pairs = ((a, b) for a in bbs for b in bbs if a["is_amine"] and b["is_amine"] and a["id"] <= b["id"])
        elif rxn_name == "n_alkylation":
            pairs = ((a, b) for a in bbs for b in bbs if a["is_amine"] and b["is_alkyl_x"])
        elif rxn_name == "o_alkylation":
            pairs = ((a, b) for a in bbs for b in bbs if a["is_oh"] and b["is_alkyl_x"])
        else:
            pairs = ()

        n_attempts = 0
        n_added_before = len(seen)
        for a, b in pairs:
            if len(seen) >= args.max_products: break
            n_attempts += 1
            try:
                outs = rxn.RunReactants((a["mol"], b["mol"]))
            except Exception:
                continue
            for prod_set in outs:
                if not prod_set: continue
                try:
                    smi = Chem.MolToSmiles(prod_set[0])
                except Exception:
                    continue
                _add(smi, rxn_name, a["id"], b["id"])
        n_added_here = len(seen) - n_added_before
        print(f"  {rxn_name:18s}  attempts={n_attempts:>6d}  new_products={n_added_here}")

    # Write
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    cols = ["id", "smiles", "reaction", "bb1", "bb2"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(products)
    print(f"\nWrote {len(products)} unique products to {args.out}")


if __name__ == "__main__":
    main()
