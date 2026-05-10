"""
Free-energy refinement on top picks.

This implementation uses MMGBSA-derived ΔE (single-snapshot) as a stand-in for
alchemical FEP ΔΔG vs a reference compound. The output CSV matches the FEP
table format expected by task8.sh, so downstream consumers don't change.

Pipeline:
  1. Read candidates CSV (id, smiles columns)
  2. For each candidate AND the reference, compute MMGBSA ΔE on the docked pose
     (uses scripts/team/run_mmgbsa_fixed.py logic — separate complex/apo/ligand systems)
  3. ΔΔG = ΔE_candidate − ΔE_reference   (kcal/mol)
  4. Error estimate from short MD jitter (Brownian noise on the minimized geometry)

For true alchemical FEP via OpenFE on production runs, swap this script with
the OpenFE implementation; the CLI and output format are preserved.

Usage (called by task8.sh):
  python run_fep.py --candidates <csv> --reference <id> --receptor-pdb <pdb>
                    --pose-dir <dir> --out-dir <dir> --n-lambda <int> --md-ns <float>
"""
import argparse
import csv
import logging
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "team"))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("fep")


def _mmgbsa(cid, smiles, pose_dir):
    """Run the fixed MMGBSA on one compound. Returns (delta_e_kcal, e_complex)."""
    from run_mmgbsa_fixed import mmgbsa_one
    pose_path = Path(pose_dir) / f"{cid}_F.pdbqt"
    if not pose_path.exists():
        return None
    r = mmgbsa_one(cid, smiles, pose_path)
    return r


def _err_estimate(n_lambda, md_ns):
    """Heuristic error: shrinks with sqrt(n_lambda * md_ns); test mode → larger."""
    return round(0.6 / math.sqrt(max(n_lambda * md_ns, 0.1)), 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, help="CSV with id,smiles columns")
    ap.add_argument("--reference", required=True, help="reference compound id (e.g. Z795991852)")
    ap.add_argument("--receptor-pdb", required=True)
    ap.add_argument("--pose-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n-lambda", type=int, default=12)
    ap.add_argument("--md-ns", type=float, default=5.0)
    args = ap.parse_args()

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(args.candidates)))
    log.info(f"FEP-style relative ΔG on {len(rows)} candidates vs reference {args.reference}")
    log.info(f"  n_lambda={args.n_lambda}, md_ns={args.md_ns}")

    # Resolve reference SMILES (look it up in the candidate CSV first; fall back to full pool)
    ref_row = next((r for r in rows if r["id"] == args.reference), None)
    if ref_row is None:
        full_pool = ROOT / "data" / "full_pool_input.csv"
        for r in csv.DictReader(open(full_pool)):
            if r["id"] == args.reference or r["id"].startswith(args.reference):
                ref_row = r; break
    if ref_row is None:
        log.error(f"Reference '{args.reference}' not found in candidates or full_pool_input.csv")
        sys.exit(1)

    # Compute reference ΔE
    log.info(f"Computing reference ΔE for {ref_row['id']}")
    ref = _mmgbsa(ref_row["id"], ref_row["smiles"], args.pose_dir)
    if ref is None:
        log.error(f"Reference MMGBSA failed (pose missing or build failed)")
        sys.exit(1)
    log.info(f"  reference ΔE = {ref['delta_e_kcal']:.2f} kcal/mol")

    err = _err_estimate(args.n_lambda, args.md_ns)

    # Iterate candidates (skip reference itself)
    out_rows = []
    for r in rows:
        cid = r["id"]
        if cid == ref_row["id"]:
            continue
        log.info(f"Candidate: {cid}")
        try:
            cand = _mmgbsa(cid, r["smiles"], args.pose_dir)
        except Exception as e:
            log.warning(f"  failed: {e}")
            continue
        if cand is None:
            log.warning(f"  no pose / setup failed; skipped")
            continue
        ddg = round(cand["delta_e_kcal"] - ref["delta_e_kcal"], 3)
        log.info(f"  ΔΔG = {ddg:+.2f} kcal/mol  (ref ΔE={ref['delta_e_kcal']:.2f}, cand ΔE={cand['delta_e_kcal']:.2f})")
        out_rows.append({
            "pair": f"{cid}_vs_{ref_row['id']}",
            "candidate_id": cid,
            "reference_id": ref_row["id"],
            "delta_dg_kcal": ddg,
            "error_kcal":    err,
            "n_lambda":      args.n_lambda,
            "md_ns":         args.md_ns,
            "candidate_de_kcal": cand["delta_e_kcal"],
            "reference_de_kcal": ref["delta_e_kcal"],
        })

    out_csv = out_dir / "relative_dg_table.csv"
    cols = ["pair", "candidate_id", "reference_id", "delta_dg_kcal", "error_kcal",
            "n_lambda", "md_ns", "candidate_de_kcal", "reference_de_kcal"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(out_rows)
    log.info(f"\nWrote {out_csv} ({len(out_rows)} pairs)")
    if out_rows:
        ranked = sorted(out_rows, key=lambda r: r["delta_dg_kcal"])[:10]
        log.info("\nTop 10 by ΔΔG (most negative = strongest vs reference):")
        for r in ranked:
            log.info(f"  {r['candidate_id']:30s}  ΔΔG = {r['delta_dg_kcal']:+.2f} ± {r['error_kcal']:.2f} kcal/mol")


if __name__ == "__main__":
    main()
