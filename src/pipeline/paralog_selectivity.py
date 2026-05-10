"""
T-box family selectivity analysis for the TBXT hackathon.

For each family member: pairwise-align the full sequence against TBXT DBD (residues 42-219),
walk the alignment blocks to map TBXT positions to family-equivalent positions, and read
off the residue at each TBXT site-F / site-A pocket position.

Site residues are validated against the canonical UniProt TBXT (O15178) sequence.
"""
import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner

DATA = Path(__file__).resolve().parents[1] / "data/selectivity"

TBXT_DBD_START = 42  # 1-based, inclusive
TBXT_DBD_END = 219   # 1-based, inclusive

# Site residues — validated against canonical TBXT O15178 (WT, G177)
# Note: PDB 6F59 has G→D at 177 (the variant). Canonical UniProt is G.
SITE_RESIDUES = {
    "F": {
        # primary anchors (TEP-defined)
        42:  ("L", "anchor — hydrophobic"),
        88:  ("Y", "anchor — H-bond"),
        177: ("G", "VARIANT site (G177D = chordoma SNP, our assay protein has D)"),
        # surrounding contact residues from our docking pose analysis
        81:  ("G", "loop"),
        82:  ("L", "loop hydrophobic"),
        83:  ("D", "loop polar"),
        172: ("I", "helix C-cap hydrophobic"),
        174: ("R", "salt bridge / cation"),
        181: ("M", "hydrophobic"),
        183: ("T", "polar"),
    },
    "A": {
        # site A pocket residues — verified against UniProt TBXT
        89:  ("S", "anchor polar"),
        91:  ("L", "hydrophobic"),
        120: ("P", "rigid"),
        121: ("S", "polar"),
        123: ("V", "hydrophobic"),
        125: ("I", "hydrophobic"),
        126: ("H", "polar"),
        129: ("S", "polar"),
        130: ("P", "rigid"),
        173: ("V", "hydrophobic"),
        180: ("R", "salt bridge / cation"),
    }
}

FAMILY = ["TBXT", "TBR1", "TBR2", "TBX1", "TBX2", "TBX3", "TBX5", "TBX6",
          "TBX10", "TBX15", "TBX18", "TBX19", "TBX20", "TBX21", "TBX22", "MGA"]


def load_seq(name):
    rec = next(SeqIO.parse(DATA / f"{name}.fasta", "fasta"))
    return str(rec.seq)


def map_tbxt_to_family(tbxt_dbd, family_seq, dbd_start_in_tbxt):
    """Return dict: TBXT 1-based position → (family aa, family 1-based position).

    Uses Bio.Align.PairwiseAligner.aligned blocks (gapless segments).
    """
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.match_score = 2
    aligner.mismatch_score = -1
    alignment = aligner.align(tbxt_dbd, family_seq)[0]

    # alignment.aligned = ((q_blocks), (t_blocks))
    # where each q_block = (start, end) in tbxt_dbd; each t_block = (start, end) in family_seq
    q_blocks, t_blocks = alignment.aligned
    pos_map = {}
    for (qs, qe), (ts, te) in zip(q_blocks, t_blocks):
        # qs, qe are 0-based half-open in the dbd slice
        # Map qs..qe-1 → ts..te-1 in family
        for off in range(qe - qs):
            tbxt_dbd_idx = qs + off  # 0-based in slice
            family_idx = ts + off    # 0-based in family seq
            tbxt_pos = dbd_start_in_tbxt + tbxt_dbd_idx  # 1-based TBXT position
            family_pos = family_idx + 1                  # 1-based family position
            family_aa = family_seq[family_idx]
            pos_map[tbxt_pos] = (family_aa, family_pos)
    return pos_map, alignment.score


def main():
    tbxt_seq = load_seq("TBXT")
    print(f"TBXT length: {len(tbxt_seq)} aa\n")

    # Sanity-check site residue positions
    print("Sanity-check site F + A residues against canonical TBXT:")
    for site, residues in SITE_RESIDUES.items():
        for pos, (expected, role) in residues.items():
            actual = tbxt_seq[pos - 1]
            ok = "✓" if actual == expected else "✗"
            print(f"  Site {site}  pos {pos:3d}  expected {expected}, actual {actual}  {ok}  ({role})")

    tbxt_dbd = tbxt_seq[TBXT_DBD_START - 1: TBXT_DBD_END]

    site_f_positions = sorted(SITE_RESIDUES["F"].keys())
    site_a_positions = sorted(SITE_RESIDUES["A"].keys())
    site_f_matrix = []
    site_a_matrix = []

    for name in FAMILY:
        seq = load_seq(name)
        if name == "TBXT":
            # No need to align to self — direct positional lookup
            f_row = {"family": name, "_score": "-"}
            a_row = {"family": name, "_score": "-"}
            for pos in site_f_positions:
                f_row[f"pos_{pos}"] = tbxt_seq[pos - 1]
                f_row[f"fpos_{pos}"] = pos
            for pos in site_a_positions:
                a_row[f"pos_{pos}"] = tbxt_seq[pos - 1]
                a_row[f"fpos_{pos}"] = pos
        else:
            try:
                pos_map, score = map_tbxt_to_family(tbxt_dbd, seq, TBXT_DBD_START)
            except Exception as e:
                print(f"  {name:6s}: FAIL — {e}")
                continue
            f_row = {"family": name, "_score": f"{score:.0f}"}
            a_row = {"family": name, "_score": f"{score:.0f}"}
            for pos in site_f_positions:
                if pos in pos_map:
                    aa, fp = pos_map[pos]
                    f_row[f"pos_{pos}"] = aa
                    f_row[f"fpos_{pos}"] = fp
                else:
                    f_row[f"pos_{pos}"] = "-"
                    f_row[f"fpos_{pos}"] = "-"
            for pos in site_a_positions:
                if pos in pos_map:
                    aa, fp = pos_map[pos]
                    a_row[f"pos_{pos}"] = aa
                    a_row[f"fpos_{pos}"] = fp
                else:
                    a_row[f"pos_{pos}"] = "-"
                    a_row[f"fpos_{pos}"] = "-"

        # Identity + selectivity
        f_id = sum(1 for p in site_f_positions if f_row[f"pos_{p}"] == SITE_RESIDUES["F"][p][0])
        a_id = sum(1 for p in site_a_positions if a_row[f"pos_{p}"] == SITE_RESIDUES["A"][p][0])
        f_row["identity_F"] = f"{f_id}/{len(site_f_positions)}"
        a_row["identity_A"] = f"{a_id}/{len(site_a_positions)}"
        f_row["selectivity_score"] = round(1.0 - f_id / len(site_f_positions), 3)
        a_row["selectivity_score"] = round(1.0 - a_id / len(site_a_positions), 3)

        site_f_matrix.append(f_row)
        site_a_matrix.append(a_row)

    # Print site F matrix
    print("\n" + "=" * 100)
    print("SITE F RESIDUE MATRIX  (asterisk = different from TBXT)")
    print("=" * 100)
    header = f"{'family':10s}  " + "  ".join(f"{p:>4d}" for p in site_f_positions) + "  identity  selectivity  align_score"
    print(header)
    print(f"{'TBXT-ref':10s}  " + "  ".join(f"{SITE_RESIDUES['F'][p][0]:>4s}" for p in site_f_positions) + "    --       --           --")
    for r in site_f_matrix:
        row_str = f"{r['family']:10s}  "
        for p in site_f_positions:
            aa = r[f"pos_{p}"]
            ref = SITE_RESIDUES["F"][p][0]
            if aa == ref or r["family"] == "TBXT":
                row_str += f"{aa:>4s}  "
            else:
                row_str += f"{'*'+aa:>4s}  "
        row_str += f"{r['identity_F']:>6s}    {r['selectivity_score']:>5.2f}        {r['_score']}"
        print(row_str)

    print("\n" + "=" * 100)
    print("SITE A RESIDUE MATRIX")
    print("=" * 100)
    header = f"{'family':10s}  " + "  ".join(f"{p:>4d}" for p in site_a_positions) + "  identity  selectivity"
    print(header)
    print(f"{'TBXT-ref':10s}  " + "  ".join(f"{SITE_RESIDUES['A'][p][0]:>4s}" for p in site_a_positions) + "    --       --")
    for r in site_a_matrix:
        row_str = f"{r['family']:10s}  "
        for p in site_a_positions:
            aa = r[f"pos_{p}"]
            ref = SITE_RESIDUES["A"][p][0]
            if aa == ref or r["family"] == "TBXT":
                row_str += f"{aa:>4s}  "
            else:
                row_str += f"{'*'+aa:>4s}  "
        row_str += f"{r['identity_A']:>6s}    {r['selectivity_score']:>5.2f}"
        print(row_str)

    # Save CSVs
    f_cols = ["family", "_score"] + [f"pos_{p}" for p in site_f_positions] + ["identity_F", "selectivity_score"]
    with open(DATA / "site_F_residue_matrix.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=f_cols)
        w.writeheader()
        for r in site_f_matrix:
            w.writerow({k: r.get(k, "") for k in f_cols})
    a_cols = ["family", "_score"] + [f"pos_{p}" for p in site_a_positions] + ["identity_A", "selectivity_score"]
    with open(DATA / "site_A_residue_matrix.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=a_cols)
        w.writeheader()
        for r in site_a_matrix:
            w.writerow({k: r.get(k, "") for k in a_cols})

    # Identify residues that DIFFER from TBXT in EACH family member
    print("\n" + "=" * 100)
    print("KEY SUBSTITUTIONS at site F per family member  (TBXT residue → family residue)")
    print("=" * 100)
    for r in site_f_matrix:
        if r["family"] == "TBXT": continue
        diffs = []
        for p in site_f_positions:
            aa = r[f"pos_{p}"]
            ref = SITE_RESIDUES["F"][p][0]
            if aa != ref and aa != "-":
                role = SITE_RESIDUES["F"][p][1]
                diffs.append(f"{ref}{p}{aa} ({role})")
        if diffs:
            print(f"  {r['family']:6s} ({len(diffs)}/{len(site_f_positions)} subs): {', '.join(diffs)}")
        else:
            print(f"  {r['family']:6s}: 100% identical at site F (cannot achieve selectivity)")

    print("\n" + "=" * 100)
    print("CONSERVATION SUMMARY at site F")
    print("=" * 100)
    for p in site_f_positions:
        ref = SITE_RESIDUES["F"][p][0]
        role = SITE_RESIDUES["F"][p][1]
        non_tbxt = [r for r in site_f_matrix if r["family"] != "TBXT"]
        observed = {}
        for r in non_tbxt:
            aa = r[f"pos_{p}"]
            observed.setdefault(aa, 0)
            observed[aa] += 1
        n_ref = observed.get(ref, 0)
        conservation = n_ref / len(non_tbxt)
        unique_to_tbxt = "✓ TBXT-specific" if conservation < 0.5 and ref not in [k for k, v in observed.items() if k != ref and v >= n_ref] else ""
        breakdown = ", ".join(f"{aa}={n}" for aa, n in sorted(observed.items(), key=lambda x: -x[1]))
        print(f"  {ref}{p:3d} ({role:38s})  cons={conservation:>5.0%}  observed in non-TBXT family: {breakdown}  {unique_to_tbxt}")


if __name__ == "__main__":
    main()
