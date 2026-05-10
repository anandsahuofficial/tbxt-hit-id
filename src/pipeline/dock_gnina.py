"""
GNINA dock + CNN rescore pipeline.

GNINA = Vina with CNN-based pose scoring + CNN-based affinity prediction.
Per-pose outputs:
  affinity      (kcal/mol, Vina-style)
  intramol      (kcal/mol, intramolecular energy)
  CNN pose score (0–1, higher = "more native-like" pose)
  CNN affinity   (pKd = -log10(Kd in M); higher = tighter binding)

Usage:
  python scripts/dock_gnina.py --smiles-csv input.csv --site F --out-dir out/

Strategy: --cnn_scoring rescore mode → Vina finds poses, CNN rescore them.
This is faster than --cnn_scoring all (CNN-driven sampling) and validated to
correlate better with experimental Kd than Vina alone (Sunseri & Koes 2020).

Input CSV must have columns: id, smiles
Output: out/dock_results_gnina.csv with all four scores per compound + best pose PDBQT
"""
import argparse
import csv
import json
import logging
import os
import re
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("gnina")

DOCK = Path(__file__).resolve().parents[1] / "data/dock"
RECEPTOR_PDBQT = DOCK / "receptor" / "6F59_apo.pdbqt"
GRID_DEFS = DOCK / "grid_definitions.json"
# GNINA binary path. Honors $GNINA_BIN so HPC submissions can swap in a
# singularity-wrapped binary (the bare GNINA binary may be glibc-incompatible on older HPC systems).
GNINA_BIN = Path(os.environ.get("GNINA_BIN",
    str(Path(__file__).resolve().parents[1] / "bin/gnina")))

CONDA_PREFIX = os.environ.get("CONDA_PREFIX", """)")


def smiles_to_pdbqt(smiles: str, out_path: Path) -> bool:
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


# Parse the GNINA per-mode table from stdout.
# Header looks like:
#   mode |  affinity  |  intramol  |    CNN     |   CNN
#        | (kcal/mol) | (kcal/mol) | pose score | affinity
#   -----+------------+------------+------------+----------
#       1     -7.26      -0.50      0.4012      5.43
#       ...
MODE_LINE_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$")


def parse_gnina_modes(stdout_text):
    modes = []
    for line in stdout_text.splitlines():
        m = MODE_LINE_RE.match(line)
        if m:
            modes.append({
                "mode": int(m.group(1)),
                "affinity_kcal": float(m.group(2)),
                "intramol_kcal": float(m.group(3)),
                "cnn_pose_score": float(m.group(4)),
                "cnn_affinity_pkd": float(m.group(5)),
            })
    return modes


def run_gnina(receptor, ligand, center, size, out_pose, exhaustiveness=8, num_modes=9, cnn_scoring="rescore"):
    cx, cy, cz = center
    sx, sy, sz = size
    cmd = [
        str(GNINA_BIN),
        "-r", str(receptor),
        "-l", str(ligand),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx),   "--size_y", str(sy),   "--size_z", str(sz),
        "--cnn_scoring", cnn_scoring,
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "-o", str(out_pose),
    ]
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{CONDA_PREFIX}/lib:" + env.get("LD_LIBRARY_PATH", "")
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600)
    if res.returncode != 0:
        return None, res.stderr
    return parse_gnina_modes(res.stdout), res.stdout


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--smiles-csv", required=True)
    p.add_argument("--site", required=True, choices=["F", "A", "G"])
    p.add_argument("--out-dir", required=True)
    p.add_argument("--exhaustiveness", type=int, default=8)
    p.add_argument("--num-modes", type=int, default=9)
    p.add_argument("--cnn-scoring", default="rescore", choices=["rescore", "refinement", "all", "metrorescore"])
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--receptor-pdbqt", default=str(RECEPTOR_PDBQT),
                   help="Path to receptor .pdbqt (default: 6F59_apo.pdbqt)")
    args = p.parse_args()

    receptor_path = Path(args.receptor_pdbqt).resolve()
    if not receptor_path.exists():
        log.error(f"Receptor not found: {receptor_path}")
        sys.exit(2)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "ligands").mkdir(exist_ok=True)
    (out / "poses").mkdir(exist_ok=True)

    grids = json.loads(GRID_DEFS.read_text())
    g = grids[args.site]
    log.info(f"Site {args.site}: center={g['center_xyz_A']} size={g['box_size_A']}")
    log.info(f"Receptor: {receptor_path.name}")
    log.info(f"GNINA: {GNINA_BIN}  cnn_scoring={args.cnn_scoring}  exh={args.exhaustiveness}")

    rows = list(csv.DictReader(open(args.smiles_csv)))
    if args.limit:
        rows = rows[: args.limit]

    cols = ["id", "smiles", "best_vina_kcal", "best_cnn_pose_score",
            "best_cnn_affinity_pkd", "best_cnn_affinity_uM", "n_modes",
            "all_vina", "all_cnn_pose", "all_cnn_pkd", "status"]
    res_csv = out / "dock_results_gnina.csv"

    # Resume: if the CSV already exists, treat the cids in it as done.
    # New run inherits prior results; we only dock cids not in the done set.
    done_ids = set()
    results = []
    if res_csv.exists():
        for r in csv.DictReader(open(res_csv)):
            done_ids.add(r["id"])
            results.append(r)
        log.info(f"Resume: {len(done_ids)} cids already in {res_csv.name}; skipping those")

    todo = [r for r in rows if r["id"] not in done_ids]
    log.info(f"Scoring {len(todo)} ligands ({len(rows) - len(todo)} skipped)")

    # Open the CSV for incremental append. Write header if file is new.
    csv_is_new = not res_csv.exists()
    csv_fh = open(res_csv, "a", newline="")
    csv_w = csv.DictWriter(csv_fh, fieldnames=cols)
    if csv_is_new:
        csv_w.writeheader()
        csv_fh.flush()

    def emit(row_dict):
        results.append(row_dict)
        csv_w.writerow({k: row_dict.get(k, "") for k in cols})
        csv_fh.flush()
        os.fsync(csv_fh.fileno())

    for i, row in enumerate(todo, 1):
        cid = row["id"]
        smi = row["smiles"]
        log.info(f"[{i}/{len(todo)}] {cid}")

        lig_pdbqt = out / "ligands" / f"{cid}.pdbqt"
        pose_pdbqt = out / "poses" / f"{cid}_{args.site}.pdbqt"
        if not smiles_to_pdbqt(smi, lig_pdbqt):
            log.warning(f"  prep failed")
            emit({"id": cid, "smiles": smi, "status": "prep_failed"})
            continue

        modes, log_text = run_gnina(
            receptor_path, lig_pdbqt, g["center_xyz_A"], g["box_size_A"],
            pose_pdbqt, exhaustiveness=args.exhaustiveness, num_modes=args.num_modes,
            cnn_scoring=args.cnn_scoring,
        )
        if not modes:
            log.warning(f"  GNINA produced no modes")
            emit({"id": cid, "smiles": smi, "status": "no_modes"})
            continue

        # Best by CNN pose score (most native-like pose); also report best Vina + best CNN affinity
        best_cnn_pose = max(modes, key=lambda m: m["cnn_pose_score"])
        best_vina = min(modes, key=lambda m: m["affinity_kcal"])
        best_cnn_aff = max(modes, key=lambda m: m["cnn_affinity_pkd"])
        log.info(
            f"  Vina best={best_vina['affinity_kcal']:.2f}  "
            f"CNN_pose best={best_cnn_pose['cnn_pose_score']:.3f}  "
            f"CNN_pKd best={best_cnn_aff['cnn_affinity_pkd']:.2f} (={10**(-best_cnn_aff['cnn_affinity_pkd'])*1e6:.1f} µM)"
        )
        emit({
            "id": cid, "smiles": smi,
            "best_vina_kcal": best_vina["affinity_kcal"],
            "best_cnn_pose_score": best_cnn_pose["cnn_pose_score"],
            "best_cnn_affinity_pkd": best_cnn_aff["cnn_affinity_pkd"],
            "best_cnn_affinity_uM": round(10 ** (-best_cnn_aff["cnn_affinity_pkd"]) * 1e6, 2),
            "n_modes": len(modes),
            "all_vina": ";".join(f"{m['affinity_kcal']:.2f}" for m in modes),
            "all_cnn_pose": ";".join(f"{m['cnn_pose_score']:.3f}" for m in modes),
            "all_cnn_pkd": ";".join(f"{m['cnn_affinity_pkd']:.2f}" for m in modes),
            "status": "ok",
        })

    csv_fh.close()
    log.info(f"Wrote {res_csv} ({len(results)} total rows)")

    # Summary: top by CNN pKd (best Kd prediction)
    ok = [r for r in results if r.get("status") == "ok"]
    if ok:
        ok.sort(key=lambda r: -r["best_cnn_affinity_pkd"])
        log.info(f"\nTop-5 by CNN pKd (most-confident-binders):")
        for r in ok[:5]:
            log.info(f"  {r['id']:30s}  pKd={r['best_cnn_affinity_pkd']:.2f}  "
                     f"({r['best_cnn_affinity_uM']:>7.2f} µM)  vina={r['best_vina_kcal']}")


if __name__ == "__main__":
    main()
