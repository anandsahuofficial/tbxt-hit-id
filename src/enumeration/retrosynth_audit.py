# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anand Sahu and contributors
"""
Audit a candidate pool for onepot 7-reaction synthesizability - ADVISORY ONLY.

⚠ KNOWN LIMITATION: This script is currently OVER-CONSERVATIVE because the
SMARTS templates count bonds INSIDE fused heteroaromatic ring systems as
"uncovered" - but those bonds are part of BUILDING BLOCKS in onepot's
enumeration, not products of the 7 diversification reactions. So a "red"
verdict from this script does NOT mean the compound is unsynthesizable; it
just means the script's coverage logic isn't sophisticated enough to recognize
the onepot building-block ring systems.

True onepot membership can only be checked once the organizer-provided
onepot lookup mechanism is available on May 9. Until then, this script is
useful for:

  - Visualizing which bonds in a compound MATCH onepot's 7 reactions
    (n_amide, n_suzuki, n_buchwald, ... columns)
  - Flagging compounds with NO recognizable onepot reactions at all
    (those would be the most concerning, even with limitations)
  - Identifying sulfonamides (NOT in onepot's listed 7 - flagged separately)

Onepot's 3.4B-compound enumeration uses these 7 reactions:
  1. Amide coupling
  2. Suzuki–Miyaura
  3. Buchwald–Hartwig amination
  4. Urea synthesis
  5. Thiourea synthesis
  6. N-Alkylation
  7. O-Alkylation

Usage:
  python onepot_retrosynth_audit.py \
      --input  results/top_100_consensus.csv \
      --output results/top_100_onepot_audit.csv
"""
import argparse
import csv
from collections import Counter
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


# Retrosynthetic SMARTS for each reaction's product bond.
# Each entry: (reaction_name, product_SMARTS pattern matching the formed bond).
# These are pragmatic "what does the resulting bond look like" patterns -
# imperfect but conservative enough for screening.
REACTION_PATTERNS = [
    # Amide coupling: R-C(=O)-N(R')(R'')
    ("amide_coupling",         "[#6]C(=O)[N;!$(N=*);!$(NC=O)]"),
    # Sulfonamide formation (technically NOT in onepot's listed 7 reactions -
    # but we flag separately so callers can choose to allow or reject)
    ("sulfonamide_formation",  "[#6]S(=O)(=O)N"),
    # Suzuki–Miyaura: aryl–aryl C–C bond
    ("suzuki_miyaura",         "[c;R][c;R]"),
    # Buchwald–Hartwig: aryl–N bond (heteroaryl-N also ok)
    ("buchwald_hartwig",       "[c;R][N;!$(N=*);!$(NC=O);!$(NS(=O)(=O))]"),
    # Urea: N-C(=O)-N
    ("urea_synthesis",         "N[C;X3](=O)N"),
    # Thiourea: N-C(=S)-N
    ("thiourea_synthesis",     "N[C;X3](=S)N"),
    # N-alkylation (sp3 C–N) - broad
    ("n_alkylation",           "[CX4][N;!$(N=*);!$(NC=O);!$(NS(=O)(=O))]"),
    # O-alkylation (sp3 C–O–C) - broad
    ("o_alkylation",           "[CX4]O[#6;!$(C=O)]"),
]


def classify(smiles):
    """Return (covered_bonds, total_bonds, missing_bond_atoms, reaction_counts)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0, 0, [], {}

    # Heavy-atom non-ring bonds + ring-fusion bonds (the synthesis-relevant ones)
    target_bonds = set()
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue
        target_bonds.add((min(a1.GetIdx(), a2.GetIdx()),
                          max(a1.GetIdx(), a2.GetIdx())))

    covered = set()
    rxn_counts = Counter()
    for name, pat in REACTION_PATTERNS:
        patt = Chem.MolFromSmarts(pat)
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            # The first atom-atom edge in the match represents the product bond
            for i in range(len(match) - 1):
                b = (min(match[i], match[i + 1]),
                     max(match[i], match[i + 1]))
                if b in target_bonds and b not in covered:
                    covered.add(b)
                    rxn_counts[name] += 1
                    break

    # Missing bonds = bonds not matched by any pattern, restricted to those that
    # are actually retrosynth-relevant (heteroatom or aryl-aryl, not C-C in ring)
    missing = []
    for (i, j) in target_bonds - covered:
        a1, a2 = mol.GetAtomWithIdx(i), mol.GetAtomWithIdx(j)
        # Only flag heteroatom-containing bonds and aryl-aryl bonds
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond is None: continue
        if bond.GetIsAromatic() and a1.GetIsAromatic() and a2.GetIsAromatic():
            # in-ring aromatic bond - formed by ring synthesis, not 7-rxn
            continue
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6 and not bond.GetIsAromatic():
            # sp3 C-C - not a target bond for the 7 reactions
            continue
        if not (a1.GetIsAromatic() or a2.GetIsAromatic() or
                a1.GetAtomicNum() != 6 or a2.GetAtomicNum() != 6):
            continue
        missing.append((i, j, a1.GetSymbol(), a2.GetSymbol()))

    return len(covered), len(target_bonds), missing, dict(rxn_counts)


def score_one(cid, smiles):
    cov, total, missing, rxns = classify(smiles)
    onepot_only_cov = sum(rxns.get(r, 0) for r in
                          ["amide_coupling", "suzuki_miyaura", "buchwald_hartwig",
                           "urea_synthesis", "thiourea_synthesis",
                           "n_alkylation", "o_alkylation"])
    has_sulfonamide = rxns.get("sulfonamide_formation", 0) > 0

    # Verdict:
    #   green  = all retrosynth-relevant bonds covered by the 7 onepot reactions
    #   amber  = covered, but uses sulfonamide (not in onepot's listed 7)
    #   red    = at least one heteroatom/aryl-aryl bond NOT covered by any of 7
    if missing:
        verdict = "red"
        reason = f"{len(missing)} bond(s) outside the 7 onepot reactions"
    elif has_sulfonamide:
        verdict = "amber"
        reason = "uses sulfonamide formation (not in onepot's listed 7 reactions)"
    else:
        verdict = "green"
        reason = "all bonds covered by 7 onepot reactions"

    return {
        "id": cid,
        "smiles": smiles,
        "verdict": verdict,
        "reason": reason,
        "n_bonds_total": total,
        "n_bonds_covered_by_onepot7": onepot_only_cov,
        "n_bonds_uncovered": len(missing),
        "uses_sulfonamide": int(has_sulfonamide),
        "n_amide": rxns.get("amide_coupling", 0),
        "n_suzuki": rxns.get("suzuki_miyaura", 0),
        "n_buchwald": rxns.get("buchwald_hartwig", 0),
        "n_urea": rxns.get("urea_synthesis", 0),
        "n_thiourea": rxns.get("thiourea_synthesis", 0),
        "n_n_alkyl": rxns.get("n_alkylation", 0),
        "n_o_alkyl": rxns.get("o_alkylation", 0),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",  required=True, help="CSV with 'id' + 'smiles' columns")
    ap.add_argument("--output", required=True, help="audit CSV to write")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.input)))
    print(f"Auditing {len(rows)} compounds against onepot 7-reaction set...")

    results = [score_one(r["id"], r["smiles"]) for r in rows]
    cols = list(results[0].keys())
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(results)

    counts = Counter(r["verdict"] for r in results)
    print(f"\nVerdict summary:")
    print(f"  green (clean for onepot 7):  {counts.get('green', 0)}")
    print(f"  amber (sulfonamide):         {counts.get('amber', 0)}")
    print(f"  red (uncovered bonds):       {counts.get('red', 0)}")

    print(f"\nWrote {args.output}")
    print(f"\nFirst 10 results:")
    for r in results[:10]:
        flag = {"green": "✓", "amber": "⚠", "red": "✗"}[r["verdict"]]
        print(f"  {flag} {r['id']:35s} {r['verdict']:6s} {r['reason']}")


if __name__ == "__main__":
    main()
