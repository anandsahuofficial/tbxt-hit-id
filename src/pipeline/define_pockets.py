"""
Define Vina docking grids for sites F and A on the prepped 6F59 chain A receptor.

Pocket residue lists from the TEP datasheet (Fig 2 + bound-fragment captions):
  - Site F: Y88, D177, L42  (engages the variant residue D177)
  - Site A: S89, P130, V173, S129, D120, R180, L91, V123, H125, H126

For each site: grid center = centroid of pocket-residue Cα atoms; box size = 22 Å cube
(generous; covers fragment-to-leadlike compounds without being so large that scoring degrades).

Also overlay 5QS9 (G177D + bound fragments at sites A–E) onto 6F59 chain A and extract
co-crystal fragment positions to validate the grid covers them.

Output: data/dock/grid_definitions.json + grid_definitions.txt (human-readable)
"""
import json
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, Superimposer

DOCK = Path(__file__).resolve().parents[1] / "data/dock"
RCPT = DOCK / "receptor"
APO_PDB = RCPT / "6F59_apo.pdb"

POCKET_RESIDUES = {
    "F": [42, 88, 177],
    "A": [89, 91, 120, 123, 125, 126, 129, 130, 173, 180],
}

BOX_SIZE_A = 22.0  # cube edge in Å


def cα_centroid(structure, resnums):
    """Return Cα centroid of given residue numbers in chain A."""
    coords = []
    for r in structure[0]["A"].get_residues():
        if r.id[1] in resnums and "CA" in r:
            coords.append(r["CA"].coord)
    if not coords:
        raise ValueError(f"No Cα atoms found for residues {resnums}")
    return np.mean(coords, axis=0)


def main():
    parser = PDBParser(QUIET=True)
    apo = parser.get_structure("apo", str(APO_PDB))

    grids = {}
    print("=" * 70)
    print("Pocket grid definitions (6F59 chain A, G177D variant)")
    print("=" * 70)
    for site, resnums in POCKET_RESIDUES.items():
        ctr = cα_centroid(apo, resnums)
        grids[site] = {
            "anchor_residues": resnums,
            "center_xyz_A": [round(float(x), 3) for x in ctr],
            "box_size_A": [BOX_SIZE_A] * 3,
            "exhaustiveness": 16,   # Vina exhaustiveness for production runs
            "num_modes": 9,
        }
        print(f"\nSite {site}")
        print(f"  Anchor residues (chain A): {resnums}")
        print(f"  Cα centroid (Å):           x={ctr[0]:8.3f}  y={ctr[1]:8.3f}  z={ctr[2]:8.3f}")
        print(f"  Box size (Å):              {BOX_SIZE_A} × {BOX_SIZE_A} × {BOX_SIZE_A}")

    # Try overlay 5QS9 bound fragments to validate site A grid placement
    qs9_path = RCPT / "5QS9.pdb"
    if qs9_path.exists():
        print("\n" + "=" * 70)
        print("Overlay 5QS9 (G177D fragment-bound) onto 6F59 chain A")
        print("=" * 70)
        qs9 = parser.get_structure("qs9", str(qs9_path))
        # Find chain A in 5QS9 (may differ); pick the chain with most matching residues
        best_chain = None
        best_n = 0
        for ch in qs9[0]:
            n = sum(1 for r in ch.get_residues()
                    if r.get_resname() not in {"HOH"} and "CA" in r)
            if n > best_n:
                best_n = n
                best_chain = ch.id
        print(f"  Chain in 5QS9 with most Cα atoms: {best_chain} ({best_n} CAs)")

        # Match by residue id (1-to-1 across the Cα atoms common to both)
        apo_cas = {r.id[1]: r["CA"] for r in apo[0]["A"].get_residues() if "CA" in r}
        qs9_cas = {r.id[1]: r["CA"] for r in qs9[0][best_chain].get_residues() if "CA" in r and r.get_resname() not in {"HOH"}}
        common = sorted(set(apo_cas) & set(qs9_cas))
        print(f"  Common Cα residues for alignment: {len(common)}")

        sup = Superimposer()
        sup.set_atoms([apo_cas[i] for i in common], [qs9_cas[i] for i in common])
        print(f"  RMSD after Cα superposition: {sup.rms:.3f} Å")
        sup.apply(qs9.get_atoms())  # transform 5QS9 in-place

        # Find the bound fragment ligand (HETATM, not HOH/CD, with >5 atoms)
        from Bio.PDB import is_aa
        frag_centroids = []
        for ch in qs9[0]:
            for r in ch.get_residues():
                if is_aa(r) or r.get_resname() in {"HOH", "CD", "EDO", "GOL", "DMS"}:
                    continue
                atoms = [a for a in r.get_atoms()]
                if len(atoms) < 5: continue
                c = np.mean([a.coord for a in atoms], axis=0)
                frag_centroids.append((r.get_resname(), ch.id, r.id[1], c))

        for ccd, ch, rn, c in frag_centroids:
            # Distance to each pocket center
            dists = {site: float(np.linalg.norm(c - np.array(g["center_xyz_A"]))) for site, g in grids.items()}
            nearest = min(dists, key=dists.get)
            print(f"  Bound fragment {ccd} (chain {ch} resi {rn}): "
                  f"distances to pocket centers F={dists['F']:.2f} Å, A={dists['A']:.2f} Å "
                  f"→ nearest: site {nearest} ({dists[nearest]:.2f} Å)")

    # Save grid definitions
    out_json = DOCK / "grid_definitions.json"
    with open(out_json, "w") as f:
        json.dump(grids, f, indent=2)
    print(f"\nWrote: {out_json}")

    out_txt = DOCK / "grid_definitions.txt"
    with open(out_txt, "w") as f:
        for site, g in grids.items():
            f.write(f"# Site {site}\n")
            f.write(f"center_x = {g['center_xyz_A'][0]}\n")
            f.write(f"center_y = {g['center_xyz_A'][1]}\n")
            f.write(f"center_z = {g['center_xyz_A'][2]}\n")
            f.write(f"size_x = {g['box_size_A'][0]}\n")
            f.write(f"size_y = {g['box_size_A'][1]}\n")
            f.write(f"size_z = {g['box_size_A'][2]}\n\n")
    print(f"Wrote: {out_txt}")


if __name__ == "__main__":
    main()
