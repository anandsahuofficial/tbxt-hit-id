"""
Ensemble docking: dock each compound against all prepped receptor conformations,
compute consensus scores.

Usage:
  python scripts/dock_ensemble.py --smiles-csv input.csv --site F --out-dir out/

Outputs:
  out/<receptor>/dock_results.csv  — per-receptor scores
  out/ensemble_summary.csv         — per-compound consensus (mean, median, min, count_negative)
"""
import argparse
import csv
import json
import logging
import statistics
import sys
import traceback
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("ensemble")

DOCK = Path(__file__).resolve().parents[1] / "data/dock"
ENSEMBLE_DIR = DOCK / "receptor" / "ensemble"
ENSEMBLE_META = ENSEMBLE_DIR / "ensemble_grids.json"


def smiles_to_pdbqt(smiles, out_path):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return False
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True) != 0:
        if AllChem.EmbedMolecule(mol, randomSeed=43, useRandomCoords=True) != 0:
            return False
    try: AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception: pass
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups: return False
    pdbqt_str, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok: return False
    out_path.write_text(pdbqt_str)
    return True


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--smiles-csv", required=True)
    p.add_argument("--site", required=True, choices=["F", "A"])
    p.add_argument("--out-dir", required=True)
    p.add_argument("--exhaustiveness", type=int, default=8)
    p.add_argument("--n-modes", type=int, default=9)
    p.add_argument("--limit", type=int, default=0)
    args = p.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "ligands").mkdir(exist_ok=True)

    # Load ensemble metadata
    receptors = json.loads(ENSEMBLE_META.read_text())
    receptors = [r for r in receptors if r["ok"] and args.site in r["grids"]]
    log.info(f"Ensemble: {len(receptors)} receptors with site {args.site} grid")
    for r in receptors:
        g = r["grids"][args.site]
        log.info(f"  {r['name']:8s}  center={g['center_xyz_A']}  ({r['label']})")

    rows = list(csv.DictReader(open(args.smiles_csv)))
    if args.limit: rows = rows[:args.limit]
    log.info(f"\nDocking {len(rows)} ligands against {len(receptors)} receptors at site {args.site}\n")

    # Prepare ligand PDBQTs once
    log.info("Preparing ligand PDBQTs...")
    valid_rows = []
    for row in rows:
        cid = row["id"]
        lig = out / "ligands" / f"{cid}.pdbqt"
        if not lig.exists():
            if not smiles_to_pdbqt(row["smiles"], lig):
                log.warning(f"  prep failed: {cid}")
                continue
        valid_rows.append(row)
    log.info(f"  {len(valid_rows)}/{len(rows)} ligands prepped")

    # Per-receptor docking
    per_receptor = {}  # cid → receptor → score
    for rec in receptors:
        rec_name = rec["name"]
        log.info(f"\n--- Docking against {rec_name} ---")
        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(rec["pdbqt"])
        g = rec["grids"][args.site]
        v.compute_vina_maps(center=g["center_xyz_A"], box_size=g["box_size_A"])
        rec_dir = out / rec_name
        rec_dir.mkdir(exist_ok=True)
        (rec_dir / "poses").mkdir(exist_ok=True)

        rec_results = []
        for row in valid_rows:
            cid = row["id"]
            try:
                v.set_ligand_from_file(str(out / "ligands" / f"{cid}.pdbqt"))
                v.dock(exhaustiveness=args.exhaustiveness, n_poses=args.n_modes)
                v.write_poses(str(rec_dir / "poses" / f"{cid}_{args.site}.pdbqt"),
                              n_poses=args.n_modes, overwrite=True)
                energies = v.energies(n_poses=args.n_modes)
                best = float(min(e[0] for e in energies))
                per_receptor.setdefault(cid, {})[rec_name] = best
                rec_results.append({"id": cid, "smiles": row["smiles"],
                                    "best_score": best, "all_scores": ";".join(f"{e[0]:.2f}" for e in energies)})
            except Exception as e:
                log.warning(f"  {cid} failed on {rec_name}: {e}")
                rec_results.append({"id": cid, "smiles": row["smiles"], "best_score": None, "all_scores": ""})

        with open(rec_dir / "dock_results.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["id", "smiles", "best_score", "all_scores"])
            w.writeheader(); w.writerows(rec_results)

    # Consensus per compound
    log.info("\n--- Consensus scores ---")
    consensus = []
    for cid, scores_by_rec in per_receptor.items():
        scores = list(scores_by_rec.values())
        valid = [s for s in scores if s is not None]
        if not valid:
            consensus.append({"id": cid, "n_receptors": 0})
            continue
        consensus.append({
            "id": cid,
            "smiles": next((r["smiles"] for r in valid_rows if r["id"] == cid), ""),
            "n_receptors": len(valid),
            "mean": round(statistics.mean(valid), 2),
            "median": round(statistics.median(valid), 2),
            "min": round(min(valid), 2),
            "max": round(max(valid), 2),
            "stdev": round(statistics.stdev(valid), 2) if len(valid) > 1 else 0,
            "count_le_minus6": sum(1 for s in valid if s <= -6.0),
            "all_scores": "; ".join(f"{r}={scores_by_rec.get(r, 'n/a')}" for r in [x["name"] for x in receptors]),
        })

    consensus.sort(key=lambda x: x.get("mean", 999))
    cols = ["id", "smiles", "n_receptors", "mean", "median", "min", "max", "stdev", "count_le_minus6", "all_scores"]
    with open(out / "ensemble_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in consensus:
            w.writerow({k: r.get(k, "") for k in cols})

    log.info(f"\nWrote {out / 'ensemble_summary.csv'}")
    log.info(f"\nTop by mean ensemble score:")
    log.info(f"  {'id':30s}  {'mean':>6}  {'median':>7}  {'min':>6}  {'max':>6}  {'stdev':>5}  {'pos≤-6':>7}")
    for r in consensus[:15]:
        if r.get("n_receptors", 0) == 0: continue
        log.info(f"  {r['id']:30s}  {r['mean']:>6.2f}  {r['median']:>7.2f}  "
                 f"{r['min']:>6.2f}  {r['max']:>6.2f}  {r['stdev']:>5.2f}  "
                 f"{r['count_le_minus6']:>3}/{r['n_receptors']}")


if __name__ == "__main__":
    main()
