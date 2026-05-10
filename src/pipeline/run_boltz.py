"""
Boltz-2 co-folding for TBXT G177D + ligand pairs.

Input: SMILES CSV with columns id, smiles
Output:
  data/boltz/<id>.yaml          — Boltz input
  data/boltz/runs/<id>/         — co-folded structure(s) + confidence + affinity
  data/boltz/boltz_summary.csv  — per-compound: pLDDT, ipTM, predicted affinity, etc.
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

DOCK = Path(__file__).resolve().parents[1] / "data/dock"
BOLTZ_OUT = Path(__file__).resolve().parents[1] / "data/boltz"
BOLTZ_OUT.mkdir(exist_ok=True)
(BOLTZ_OUT / "yaml").mkdir(exist_ok=True)
(BOLTZ_OUT / "runs").mkdir(exist_ok=True)

# TBXT G177D DBD — extracted directly from 6F59 chain A (PDB residues 41–224, 178 aa).
# Verified G→D at position 177.
TBXT_G177D_DBD = (
    "ELRVGLEESELWLRFKELTNEMIVTKNGRRMFPVLKVNVSGLDPNAMYSFLLDFVAADNHRWKYVNGEWVP"
    "QAPSCVYIHPDSPNFGAHWMKAPVSFSKVKLTNKLNGGGQIMLNSLHKYEPRIHIVRVGDPQRMITSHCFP"
    "ETQFIAVTAYQNEEITALKIKYNPFAKAFLDAKERS"
)

YAML_TEMPLATE = """version: 1
sequences:
  - protein:
      id: A
      sequence: {sequence}
      msa: empty
  - ligand:
      id: L
      smiles: '{smiles}'
properties:
  - affinity:
      binder: L
"""


def write_yaml(cid, smiles):
    yaml = YAML_TEMPLATE.format(sequence=TBXT_G177D_DBD, smiles=smiles)
    path = BOLTZ_OUT / "yaml" / f"{cid}.yaml"
    path.write_text(yaml)
    return path


def _detect_accelerator():
    try:
        import torch
        return ("gpu", "1") if torch.cuda.is_available() else ("cpu", "1")
    except Exception:
        return ("cpu", "1")


def run_boltz(yaml_path, out_dir, fast=False):
    accelerator, devices = _detect_accelerator()
    safeload = Path(__file__).resolve().parent / "_boltz_safeload.py"
    # Fast mode reduces samples + sampling steps for smoke-test on CPU.
    diffusion = "1" if fast else "3"
    recycling = "1" if fast else "3"
    sampling  = "10" if fast else "200"
    # CPU diffusion at 178-aa is ~3-5 min/step; allow 30 min for fast smoke-test
    # (production GPU should finish in <2 min/compound).
    timeout   = 1800 if fast else 3600
    cmd = [
        sys.executable, str(safeload), "predict", str(yaml_path),
        "--out_dir", str(out_dir),
        "--accelerator", accelerator,
        "--devices", devices,
        "--diffusion_samples", diffusion,
        "--recycling_steps", recycling,
        "--sampling_steps", sampling,
        "--output_format", "pdb",
        "--write_full_pae",
        "--seed", "42",
        "--override",
    ]
    print(f"  Running: boltz predict {yaml_path.name} ... "
          f"[accelerator={accelerator}, fast={fast}, samples={diffusion}, "
          f"recycle={recycling}, sampling={sampling}]")
    rc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return rc


def parse_results(out_dir, cid):
    """Boltz writes results to <out_dir>/boltz_results_<cid>/predictions/<cid>/."""
    pred_dir = out_dir / f"boltz_results_{cid}" / "predictions" / cid
    if not pred_dir.exists(): return None
    info = {"cid": cid, "pred_dir": str(pred_dir)}

    # Average confidence across the diffusion samples
    conf_files = sorted(pred_dir.glob(f"confidence_{cid}_model_*.json"))
    confs = []
    for f in conf_files:
        try:
            confs.append(json.loads(f.read_text()))
        except Exception:
            pass
    if confs:
        # Take model_0 (top-ranked by Boltz)
        d = confs[0]
        info["pLDDT"] = round(d.get("complex_plddt", 0), 4)
        info["pTM"] = round(d.get("ptm", 0), 4)
        info["ipTM"] = round(d.get("iptm", 0), 4)
        info["confidence"] = round(d.get("confidence_score", 0), 4)
        info["lig_iptm"] = round(d.get("ligand_iptm", 0), 4)
        info["n_models"] = len(confs)
        # Best across models for ipTM
        info["ipTM_best"] = round(max(c.get("iptm", 0) for c in confs), 4)
        info["confidence_best"] = round(max(c.get("confidence_score", 0) for c in confs), 4)

    # Affinity: <cid>_affinity.json (one per ligand, single file)
    aff_file = pred_dir / f"affinity_{cid}.json"
    if aff_file.exists():
        try:
            d = json.loads(aff_file.read_text())
            # affinity_pred_value is log10(Kd in µM). Convert to Kd µM for readability.
            log_kd_uM = d.get("affinity_pred_value")
            if log_kd_uM is not None:
                info["affinity_log_kd_uM"] = round(log_kd_uM, 3)
                info["affinity_kd_uM"] = round(10 ** log_kd_uM, 3)
                info["affinity_pkd"] = round(6.0 - log_kd_uM, 3)  # pKd in M
            info["affinity_prob_binder"] = round(d.get("affinity_probability_binary", 0), 4)
        except Exception as e:
            info["affinity_err"] = str(e)

    info["n_pdbs"] = len(list(pred_dir.glob("*.pdb")))
    return info


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--smiles-csv", required=True)
    p.add_argument("--out-dir", default=str(BOLTZ_OUT / "runs"))
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--fast", action="store_true",
                   help="reduce samples/sampling steps for smoke-test on CPU")
    args = p.parse_args()

    rows = list(csv.DictReader(open(args.smiles_csv)))
    if args.limit: rows = rows[:args.limit]
    print(f"Co-folding {len(rows)} compounds with TBXT G177D ({len(TBXT_G177D_DBD)} aa)")

    out_root = Path(args.out_dir)

    cols = ["cid", "status", "smiles", "pLDDT", "pTM", "ipTM", "ipTM_best",
            "lig_iptm", "confidence", "confidence_best",
            "affinity_log_kd_uM", "affinity_kd_uM", "affinity_pkd",
            "affinity_prob_binder", "elapsed_s", "n_pdbs", "n_models", "pred_dir"]
    out_csv = BOLTZ_OUT / "boltz_summary.csv"

    # Resume: cids already in the summary CSV with status=ok are skipped.
    # Failed/no-predictions cids are RETRIED on resume.
    done_ids = set()
    summary = []
    if out_csv.exists():
        for r in csv.DictReader(open(out_csv)):
            if r.get("status") == "ok":
                done_ids.add(r["cid"])
                summary.append(r)
        print(f"Resume: {len(done_ids)} cids already OK in {out_csv.name}; skipping those")

    todo = [r for r in rows if r["id"] not in done_ids]
    print(f"To run: {len(todo)} compounds ({len(rows) - len(todo)} skipped)")

    # Open the CSV for incremental append. Rewrite header + prior OK rows so
    # the file is well-formed in either case (new file OR existing partial).
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    csv_fh = open(out_csv, "w", newline="")
    csv_w = csv.DictWriter(csv_fh, fieldnames=cols)
    csv_w.writeheader()
    for r in summary:
        csv_w.writerow({k: r.get(k, "") for k in cols})
    csv_fh.flush()
    os.fsync(csv_fh.fileno())

    def emit(row_dict):
        summary.append(row_dict)
        csv_w.writerow({k: row_dict.get(k, "") for k in cols})
        csv_fh.flush()
        os.fsync(csv_fh.fileno())

    for i, row in enumerate(todo, 1):
        cid = row["id"]
        smiles = row["smiles"]
        print(f"\n[{i}/{len(todo)}] {cid} : {smiles}")
        yaml_path = write_yaml(cid, smiles)
        out_dir = out_root / cid
        out_dir.mkdir(parents=True, exist_ok=True)

        t0 = time.time()
        rc = run_boltz(yaml_path, out_dir, fast=args.fast)
        elapsed = time.time() - t0
        if rc.returncode != 0:
            print(f"  FAIL ({elapsed:.0f}s):")
            print(f"  stderr: {rc.stderr[-500:]}")
            emit({"cid": cid, "status": "fail", "elapsed_s": round(elapsed, 1),
                  "err": rc.stderr[-500:]})
            continue

        info = parse_results(out_dir, Path(yaml_path).stem)
        if info is None:
            print(f"  No predictions parsed ({elapsed:.0f}s)")
            emit({"cid": cid, "status": "no_predictions", "elapsed_s": round(elapsed, 1)})
            continue
        info["cid"] = cid
        info["status"] = "ok"
        info["elapsed_s"] = round(elapsed, 1)
        info["smiles"] = smiles
        print(f"  OK ({elapsed:.0f}s): pLDDT={info.get('pLDDT', 'NA')}, "
              f"ipTM={info.get('ipTM', 'NA')}, "
              f"affinity={info.get('affinity_kd_uM', 'NA')} µM, "
              f"prob_binder={info.get('affinity_prob_binder', 'NA')}")
        emit(info)

    csv_fh.close()
    print(f"\nWrote {out_csv} ({len(summary)} total rows)")


if __name__ == "__main__":
    main()
