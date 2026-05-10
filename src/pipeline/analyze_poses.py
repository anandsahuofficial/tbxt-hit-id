"""
Analyze docked poses: for each compound + site, parse the top pose and report
which receptor residues are within 4 Å of any ligand heavy atom.

This validates that poses actually contact the expected anchor residues
(Y88, D177, L42 for site F; S89, V173, etc. for site A).
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

DOCK = Path(__file__).resolve().parents[1] / "data/dock"
RECEPTOR = DOCK / "receptor" / "6F59_apo.pdb"
CONTACT_CUTOFF = 4.0  # Å


def parse_pdbqt_top_pose(pdbqt_path: Path):
    """Yield first MODEL block atom coords + atom names from a Vina-output PDBQT."""
    coords = []
    atom_names = []
    in_model = False
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith("MODEL "):
            if in_model: break
            in_model = True
            continue
        if line.startswith("ENDMDL"):
            break
        if not in_model: continue
        if line.startswith(("ATOM", "HETATM")):
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                element = line[12:16].strip()[0]
                coords.append((x, y, z))
                atom_names.append(element)
            except (ValueError, IndexError):
                pass
    return np.array(coords) if coords else np.zeros((0, 3)), atom_names


def get_receptor_atoms():
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("rec", str(RECEPTOR))
    atoms = []
    for r in s[0]["A"].get_residues():
        rn = r.id[1]
        rname = r.get_resname()
        for a in r.get_atoms():
            if a.element != "H":
                atoms.append((rn, rname, a.get_name(), a.coord))
    return atoms


def contacts(ligand_coords, receptor_atoms, cutoff=CONTACT_CUTOFF):
    """Return sorted list of (resnum, resname, min_dist) where min_dist <= cutoff."""
    if len(ligand_coords) == 0: return []
    rec_coords = np.array([a[3] for a in receptor_atoms])
    # pairwise distances: (N_lig, N_rec)
    d = np.linalg.norm(ligand_coords[:, None, :] - rec_coords[None, :, :], axis=-1)
    min_per_rec = d.min(axis=0)  # closest ligand atom to each receptor atom
    by_res = {}
    for i, m in enumerate(min_per_rec):
        if m > cutoff: continue
        rn, rname, aname, _ = receptor_atoms[i]
        prev = by_res.get(rn)
        if prev is None or m < prev[2]:
            by_res[rn] = (rn, rname, float(m))
    return sorted(by_res.values())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--validation-dirs", nargs="+", default=[
        str(DOCK / "validation_F"), str(DOCK / "validation_A")
    ])
    args = p.parse_args()

    rec_atoms = get_receptor_atoms()
    print(f"Receptor: {len(rec_atoms)} heavy atoms parsed\n")

    SITE_F_ANCHORS = {42, 88, 177}
    SITE_A_ANCHORS = {89, 91, 120, 123, 125, 126, 129, 130, 173, 180}

    summary_rows = []
    for vdir in args.validation_dirs:
        vdir = Path(vdir)
        site = "F" if vdir.name.endswith("_F") else "A"
        anchors = SITE_F_ANCHORS if site == "F" else SITE_A_ANCHORS
        results = list(csv.DictReader(open(vdir / "dock_results.csv")))
        print("=" * 78)
        print(f"Site {site} pose contacts (cutoff {CONTACT_CUTOFF} Å)")
        print("=" * 78)
        print(f"  Anchor residues: {sorted(anchors)}")
        for r in results:
            cid = r["id"]
            score = r["best_score"]
            pose_path = vdir / "poses" / f"{cid}_{site}.pdbqt"
            if not pose_path.exists():
                print(f"\n  {cid}: pose file missing")
                continue
            coords, _ = parse_pdbqt_top_pose(pose_path)
            cl = contacts(coords, rec_atoms)
            anchor_hit = sorted({x[0] for x in cl} & anchors)
            anchor_dists = {x[0]: x[2] for x in cl if x[0] in anchors}
            anchor_str = ", ".join(f"{rn}({anchor_dists[rn]:.2f}Å)" for rn in anchor_hit) or "—"
            other_hit = [x for x in cl if x[0] not in anchors]
            print(f"\n  {cid:14s}  score={score} kcal/mol  | anchors hit: {anchor_str}")
            print(f"    {len(cl)} contact residues; {len(other_hit)} non-anchor")
            top_other = sorted(other_hit, key=lambda x: x[2])[:5]
            if top_other:
                print(f"    nearest non-anchor: " + ", ".join(f"{rn}{rname}({d:.2f}Å)" for rn, rname, d in top_other))
            summary_rows.append({
                "id": cid, "site": site, "vina_score": score,
                "anchors_hit": ";".join(map(str, anchor_hit)),
                "n_contacts": len(cl),
            })

    # Save summary
    with open(DOCK / "pose_contacts_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "site", "vina_score", "anchors_hit", "n_contacts"])
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\nWrote {DOCK}/pose_contacts_summary.csv")


if __name__ == "__main__":
    main()
