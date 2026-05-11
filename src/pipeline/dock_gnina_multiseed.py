# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anand Sahu and contributors
"""
Multi-seed GNINA docking - averages CNN pose / Vina / pKd scores across N seeds
to reduce run-to-run variance.

Usage:
  python scripts/team/dock_gnina_multiseed.py \
      --smiles-csv data/full_pool_input.csv \
      --site F \
      --out-dir data/full_pool_gnina_F_multiseed \
      --seeds 10 \
      --exhaustiveness 8

For sharding across multiple GPUs/nodes:
  python scripts/team/dock_gnina_multiseed.py ... --start-idx 0   --end-idx 285
  python scripts/team/dock_gnina_multiseed.py ... --start-idx 285 --end-idx 570
"""
import argparse
import csv
import json
import logging
import os
import re
import statistics
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("multiseed")

DOCK = Path(__file__).resolve().parents[2] / "data" / "dock"
RECEPTOR_PDBQT = DOCK / "receptor" / "6F59_apo.pdbqt"
GRID_DEFS = DOCK / "grid_definitions.json"
GNINA_BIN = Path(__file__).resolve().parents[2] / "bin" / "gnina"
CONDA_PREFIX = os.environ.get("CONDA_PREFIX", os.environ.get("HOME", ""))

MODE_LINE_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$")


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
    pdbqt_str, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
    if not ok: return False
    out_path.write_text(pdbqt_str)
    return True


def parse_modes(stdout_text):
    out = []
    for line in stdout_text.splitlines():
        m = MODE_LINE_RE.match(line)
        if m:
            out.append({
                "mode": int(m.group(1)),
                "affinity": float(m.group(2)),
                "intramol": float(m.group(3)),
                "cnn_pose": float(m.group(4)),
                "cnn_pkd": float(m.group(5)),
            })
    return out


def run_gnina(receptor, ligand, center, size, out_pose, seed, exhaustiveness, num_modes=9, cnn_scoring="rescore"):
    cx, cy, cz = center
    sx, sy, sz = size
    cmd = [
        str(GNINA_BIN), "-r", str(receptor), "-l", str(ligand),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--cnn_scoring", cnn_scoring,
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        "-o", str(out_pose),
    ]
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{CONDA_PREFIX}/lib:" + env.get("LD_LIBRARY_PATH", "")
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600)
    if res.returncode != 0:
        return None
    return parse_modes(res.stdout)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--smiles-csv", required=True)
    p.add_argument("--site", required=True, choices=["F", "A"])
    p.add_argument("--out-dir", required=True)
    p.add_argument("--seeds", type=int, default=10)
    p.add_argument("--exhaustiveness", type=int, default=8)
    p.add_argument("--start-idx", type=int, default=0)
    p.add_argument("--end-idx", type=int, default=10**9)
    args = p.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "ligands").mkdir(exist_ok=True)

    grids = json.loads(GRID_DEFS.read_text())
    g = grids[args.site]
    log.info(f"Site {args.site} center={g['center_xyz_A']} size={g['box_size_A']}")
    log.info(f"Seeds: {args.seeds}  Exhaustiveness: {args.exhaustiveness}")

    rows = list(csv.DictReader(open(args.smiles_csv)))
    rows = rows[args.start_idx : args.end_idx]
    log.info(f"Processing {len(rows)} compounds (idx {args.start_idx}–{args.start_idx+len(rows)})")

    results = []
    for i, row in enumerate(rows, 1):
        cid = row["id"]; smi = row["smiles"]
        log.info(f"[{i}/{len(rows)}] {cid}")

        lig = out / "ligands" / f"{cid}.pdbqt"
        if not lig.exists():
            if not smiles_to_pdbqt(smi, lig):
                log.warning(f"  prep failed")
                results.append({"id": cid, "smiles": smi, "status": "prep_failed"})
                continue

        per_seed = []
        for seed_idx in range(args.seeds):
            seed = 42 + seed_idx * 1000
            pose = out / f"{cid}_seed{seed_idx}.pdbqt"
            modes = run_gnina(RECEPTOR_PDBQT, lig, g["center_xyz_A"], g["box_size_A"],
                              pose, seed, args.exhaustiveness)
            if not modes:
                continue
            best_pose = max(m["cnn_pose"] for m in modes)
            best_pkd = max(m["cnn_pkd"] for m in modes)
            best_vina = min(m["affinity"] for m in modes)
            per_seed.append({"cnn_pose": best_pose, "cnn_pkd": best_pkd, "vina": best_vina})
            try: pose.unlink()
            except Exception: pass

        if not per_seed:
            results.append({"id": cid, "smiles": smi, "status": "all_seeds_failed"})
            continue

        cnn_poses = [s["cnn_pose"] for s in per_seed]
        cnn_pkds = [s["cnn_pkd"] for s in per_seed]
        vinas = [s["vina"] for s in per_seed]

        r = {
            "id": cid, "smiles": smi, "status": "ok",
            "n_seeds": len(per_seed),
            "cnn_pose_mean": round(statistics.mean(cnn_poses), 4),
            "cnn_pose_stdev": round(statistics.stdev(cnn_poses) if len(cnn_poses) > 1 else 0, 4),
            "cnn_pose_min": round(min(cnn_poses), 4),
            "cnn_pose_max": round(max(cnn_poses), 4),
            "cnn_pkd_mean": round(statistics.mean(cnn_pkds), 4),
            "cnn_pkd_stdev": round(statistics.stdev(cnn_pkds) if len(cnn_pkds) > 1 else 0, 4),
            "vina_kcal_mean": round(statistics.mean(vinas), 4),
            "vina_kcal_min": round(min(vinas), 4),
            "kd_uM_from_cnn_mean": round(10**(-statistics.mean(cnn_pkds))*1e6, 2),
        }
        results.append(r)
        log.info(f"  cnn_pose mean={r['cnn_pose_mean']:.3f}±{r['cnn_pose_stdev']:.3f} | "
                 f"cnn_pKd mean={r['cnn_pkd_mean']:.2f} | vina mean={r['vina_kcal_mean']:.2f}")

    cols = ["id", "smiles", "status", "n_seeds",
            "cnn_pose_mean", "cnn_pose_stdev", "cnn_pose_min", "cnn_pose_max",
            "cnn_pkd_mean", "cnn_pkd_stdev",
            "vina_kcal_mean", "vina_kcal_min",
            "kd_uM_from_cnn_mean"]
    out_csv = out / "dock_results_multiseed.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k, "") for k in cols})
    log.info(f"\nWrote {out_csv}")

    ok = [r for r in results if r.get("status") == "ok"]
    if ok:
        ok.sort(key=lambda r: -r["cnn_pose_mean"])
        log.info(f"\nTop 10 by cnn_pose_mean (with stdev):")
        for r in ok[:10]:
            log.info(f"  {r['id']:30s}  pose={r['cnn_pose_mean']:.3f}±{r['cnn_pose_stdev']:.3f}  "
                     f"pKd={r['cnn_pkd_mean']:.2f}  vina={r['vina_kcal_mean']:.2f}")


if __name__ == "__main__":
    main()
