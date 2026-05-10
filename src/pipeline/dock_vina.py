"""
Vina docking pipeline against TBXT G177D (6F59 chain A) at sites F or A.

Usage:
  python scripts/dock.py --smiles-csv input.csv --site F --out-dir out/
  python scripts/dock.py --smiles-csv input.csv --site A --out-dir out/

Input CSV must have columns: id, smiles
Output: out/<id>_<site>.pdbqt + out/dock_results.csv with vina scores

Pipeline:
  1. SMILES → 3D mol via RDKit (ETKDG embedding + UFF optimize)
  2. RDKit → PDBQT via Meeko (assigns Vina atom types, partial charges, rotatable bonds)
  3. Vina docking against grid box for the chosen site
  4. Parse top-9 poses, write best score + best pose to results CSV
"""
import argparse
import csv
import json
import logging
import sys
import traceback
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("dock")

DOCK = Path(__file__).resolve().parents[1] / "data/dock"
RECEPTOR_PDBQT = DOCK / "receptor" / "6F59_apo.pdbqt"
GRID_DEFS = DOCK / "grid_definitions.json"


def smiles_to_pdbqt(smiles: str, out_path: Path) -> bool:
    """SMILES → 3D conformer → PDBQT via RDKit + Meeko. Returns True on success."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        log.warning(f"  RDKit: invalid SMILES")
        return False
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True) != 0:
        log.warning(f"  ETKDG embed failed; trying random coords")
        if AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True) != 0:
            return False
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass

    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        log.warning(f"  Meeko: no setup produced")
        return False
    pdbqt_str, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        log.warning(f"  Meeko writer: {err}")
        return False
    out_path.write_text(pdbqt_str)
    return True


def dock_one(v: Vina, ligand_pdbqt: Path, out_pdbqt: Path, exhaustiveness: int = 16, n_modes: int = 9):
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_modes)
    v.write_poses(str(out_pdbqt), n_poses=n_modes, overwrite=True)
    energies = v.energies(n_poses=n_modes)
    return [float(e[0]) for e in energies]  # affinity in kcal/mol; lower = better


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--smiles-csv", required=True, help="CSV with columns id,smiles")
    p.add_argument("--site", required=True, choices=["F", "A"])
    p.add_argument("--out-dir", required=True, help="output directory")
    p.add_argument("--exhaustiveness", type=int, default=16)
    p.add_argument("--n-modes", type=int, default=9)
    p.add_argument("--limit", type=int, default=0, help="limit input rows (0 = all)")
    args = p.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "ligands").mkdir(exist_ok=True)
    (out / "poses").mkdir(exist_ok=True)

    grids = json.loads(GRID_DEFS.read_text())
    g = grids[args.site]
    cx, cy, cz = g["center_xyz_A"]
    sx, sy, sz = g["box_size_A"]
    log.info(f"Site {args.site} grid: center=({cx},{cy},{cz}) size=({sx},{sy},{sz})")

    # Initialize Vina once and reuse — receptor and maps are slow
    log.info("Initializing Vina (vina, n_cpu=auto)...")
    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(str(RECEPTOR_PDBQT))
    v.compute_vina_maps(center=[cx, cy, cz], box_size=[sx, sy, sz])
    log.info("Vina ready.")

    # Read input
    rows = list(csv.DictReader(open(args.smiles_csv)))
    if args.limit:
        rows = rows[: args.limit]
    log.info(f"Docking {len(rows)} ligands at site {args.site}")

    results = []
    for i, row in enumerate(rows, 1):
        cid = row["id"]
        smi = row["smiles"]
        log.info(f"[{i}/{len(rows)}] {cid}  ({smi[:60]}...)")

        lig_pdbqt = out / "ligands" / f"{cid}.pdbqt"
        pose_pdbqt = out / "poses" / f"{cid}_{args.site}.pdbqt"
        try:
            ok = smiles_to_pdbqt(smi, lig_pdbqt)
            if not ok:
                results.append({"id": cid, "smiles": smi, "best_score": None, "all_scores": "", "status": "prep_failed"})
                continue
            scores = dock_one(v, lig_pdbqt, pose_pdbqt,
                              exhaustiveness=args.exhaustiveness, n_modes=args.n_modes)
            best = min(scores) if scores else None
            log.info(f"  best score: {best:.2f} kcal/mol  (n_modes={len(scores)})")
            results.append({
                "id": cid, "smiles": smi, "best_score": best,
                "all_scores": ";".join(f"{s:.2f}" for s in scores),
                "status": "ok",
            })
        except Exception as e:
            log.error(f"  DOCK FAILED: {e}")
            log.debug(traceback.format_exc())
            results.append({"id": cid, "smiles": smi, "best_score": None, "all_scores": "", "status": f"error: {e}"})

    # Write results
    res_csv = out / "dock_results.csv"
    with open(res_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "smiles", "best_score", "all_scores", "status"])
        w.writeheader()
        w.writerows(results)
    log.info(f"Wrote {res_csv}")

    # Print summary
    ok_results = [r for r in results if r["best_score"] is not None]
    if ok_results:
        ok_results.sort(key=lambda r: r["best_score"])
        log.info(f"\nTop-5 by score:")
        for r in ok_results[:5]:
            log.info(f"  {r['id']:14s}  {r['best_score']:.2f} kcal/mol")


if __name__ == "__main__":
    main()
