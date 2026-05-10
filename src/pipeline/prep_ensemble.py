"""
Prepare an ensemble of TBXT receptor conformations for docking.

Conformations:
  6F59 chain A: G177D + DNA (current default; KEEP existing prep)
  6F59 chain B: G177D + DNA, dimer partner — may have different site F geometry
  6F58 chain A: WT (G177) + DNA — for sanity comparison
  5QS9 chain A: G177D apo, fragment-bound at sites A-E
  5QSA chain A: WT apo, FM001580 (K2P) bound at site F
  5QSI chain A: WT apo, FM001452 (O1D) bound at site F

For each: extract chain A (or B for 6F59-B), strip DNA/ligands/waters,
PDBFixer + protonate at pH 7.5, write PDBQT.

Then re-define site F + site A grids on each receptor (pocket residues are the same;
only Cα coordinates differ).
"""
import json
import subprocess
from pathlib import Path

import numpy as np
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.PDB import PDBParser, Select, PDBIO, Superimposer

DOCK = Path(__file__).resolve().parents[1] / "data/dock"
RCPT = DOCK / "receptor"
ENSEMBLE = RCPT / "ensemble"
ENSEMBLE.mkdir(parents=True, exist_ok=True)

CONFS = [
    {"name": "6F59_A", "src_pdb": "6F59.pdb", "chain": "A", "label": "G177D + DNA, chain A (default)"},
    {"name": "6F59_B", "src_pdb": "6F59.pdb", "chain": "B", "label": "G177D + DNA, chain B (dimer partner)"},
    {"name": "6F58_A", "src_pdb": "6F58.pdb", "chain": "A", "label": "WT (G177) + DNA, chain A"},
    {"name": "5QS9_A", "src_pdb": "5QS9.pdb", "chain": "A", "label": "G177D apo, fragment-bound (sites A-E)"},
    {"name": "5QSA_A", "src_pdb": "5QSA.pdb", "chain": "A", "label": "WT apo, FM001580 bound at site F"},
    {"name": "5QSI_A", "src_pdb": "5QSI.pdb", "chain": "A", "label": "WT apo, FM001452 bound at site F"},
]

POCKET_RESIDUES = {
    "F": [42, 88, 177, 81, 82, 83, 172, 174, 181, 183],
    "A": [89, 91, 120, 121, 123, 125, 126, 129, 130, 173, 180],
}
BOX_SIZE = 22.0

STD_AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
          "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}


class ChainProteinSelector(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_residue(self, residue):
        return residue.parent.id == self.chain_id and residue.get_resname() in STD_AA


def prep_one(conf):
    """Returns dict with name, pdbqt_path, success."""
    src = RCPT / conf["src_pdb"]
    if not src.exists():
        return {**conf, "ok": False, "err": f"missing {src}"}

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(conf["name"], str(src))

    # Verify chain exists
    chain_ids = [c.id for c in structure[0]]
    if conf["chain"] not in chain_ids:
        return {**conf, "ok": False, "err": f"chain {conf['chain']} not found in {src.name} (have: {chain_ids})"}

    # Extract chain protein
    io = PDBIO()
    io.set_structure(structure)
    chain_pdb = ENSEMBLE / f"{conf['name']}_chain.pdb"
    io.save(str(chain_pdb), ChainProteinSelector(conf["chain"]))

    # Verify it has enough residues
    p = parser.get_structure("c", str(chain_pdb))
    n_res = sum(1 for r in p[0][conf["chain"]].get_residues())
    if n_res < 100:
        return {**conf, "ok": False, "err": f"only {n_res} residues in chain {conf['chain']}"}

    # PDBFixer
    fixer = PDBFixer(filename=str(chain_pdb))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.5)

    apo_pdb = ENSEMBLE / f"{conf['name']}_apo.pdb"
    with open(apo_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    # PDBQT
    pdbqt = ENSEMBLE / f"{conf['name']}_apo.pdbqt"
    rc = subprocess.run(
        ["obabel", str(apo_pdb), "-O", str(pdbqt), "-xr", "-p", "7.5"],
        capture_output=True, text=True, timeout=120
    )
    ok = rc.returncode == 0 and pdbqt.exists() and pdbqt.stat().st_size > 1000

    return {**conf, "n_residues": n_res, "ok": ok,
            "pdbqt": str(pdbqt), "apo_pdb": str(apo_pdb),
            "obabel_stderr": rc.stderr.strip()[:200] if rc.stderr else ""}


def compute_grids(apo_pdb, chain):
    """Compute Cα-centroid grid centers for each pocket on this conformation."""
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("c", apo_pdb)
    grids = {}
    for site, resnums in POCKET_RESIDUES.items():
        coords = []
        for r in s[0][chain].get_residues():
            if r.id[1] in resnums and "CA" in r:
                coords.append(r["CA"].coord)
        if coords:
            ctr = np.mean(coords, axis=0)
            grids[site] = {
                "anchor_residues": resnums,
                "center_xyz_A": [round(float(x), 3) for x in ctr],
                "box_size_A": [BOX_SIZE] * 3,
                "n_anchors_found": len(coords),
                "n_anchors_expected": len(resnums),
            }
    return grids


def main():
    results = []
    for conf in CONFS:
        print(f"\n=== Prepping {conf['name']}: {conf['label']} ===")
        r = prep_one(conf)
        if not r["ok"]:
            print(f"  FAIL: {r.get('err', 'unknown')}")
            results.append(r)
            continue
        print(f"  OK: {r['n_residues']} residues, PDBQT {Path(r['pdbqt']).stat().st_size} bytes")

        grids = compute_grids(r["apo_pdb"], r["chain"])
        for site, g in grids.items():
            print(f"  site {site}: center=({g['center_xyz_A'][0]:>7.3f},{g['center_xyz_A'][1]:>7.3f},{g['center_xyz_A'][2]:>7.3f})  "
                  f"anchors found {g['n_anchors_found']}/{g['n_anchors_expected']}")
        r["grids"] = grids
        results.append(r)

    # Save ensemble metadata
    out = ENSEMBLE / "ensemble_grids.json"
    with open(out, "w") as f:
        json.dump([{
            "name": r["name"], "label": r["label"], "chain": r["chain"],
            "ok": r["ok"], "pdbqt": r.get("pdbqt", ""), "grids": r.get("grids", {})
        } for r in results], f, indent=2)
    print(f"\nWrote {out}")

    # Summary
    n_ok = sum(1 for r in results if r["ok"])
    print(f"\n{n_ok}/{len(results)} receptors prepped successfully")


if __name__ == "__main__":
    main()
