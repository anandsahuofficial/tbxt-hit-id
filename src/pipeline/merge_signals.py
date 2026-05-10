"""
Strategy 9 (synthesis) — combine all orthogonal signals into a Tier-A ranking.

Inputs:
  data/analogs/all_candidates.csv          (503 enumerated analogs)
  data/generative/generative_proposals.csv (67 generative proposals)
  data/full_pool_gnina_F/dock_results_gnina.csv  (Vina + CNN_pose + CNN_pKd at site F)
  data/qsar/predictions_analogs.csv        (QSAR pKd for the 503 enumerated)
  data/generative/generative_proposals.csv (QSAR pKd for the 67 generative)
  [optional] data/boltz/boltz_summary.csv  (Boltz-2 affinity for top picks)
  data/dock/ensemble_summary.csv           (per-receptor ensemble scores — only on pre-set)
  data/selectivity/site_F_residue_matrix.csv  (TBXT-unique residue list)

Outputs:
  data/tier_a/tier_a_candidates.csv      Tier-A list ranked
  data/tier_a/all_signals.csv             every compound × every signal
  data/tier_a/TIER_A_REPORT.md           narrative summary

Tier-A scoring rule (revision 4):
  ✅ GNINA CNN pKd ≥ 4.5 AND CNN pose ≥ 0.5 AND Vina best ≤ -6.0
  ✅ QSAR pKd ≥ 4.0
  (Boltz-2: applied only to Tier-A picks, not as a filter)
  (Selectivity: rationale only, not a filter)
"""
import csv
import json
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = DATA / "tier_a"
OUT.mkdir(exist_ok=True)


def load_csv(path):
    if not Path(path).exists(): return []
    return list(csv.DictReader(open(path)))


def main():
    # 1. Combined pool (570 compounds)
    pool = load_csv(DATA / "full_pool_input.csv")
    pool_by_id = {r["id"]: r for r in pool}
    print(f"Combined pool: {len(pool)}")

    # 2. GNINA full-pool scores
    gnina = {}
    gnina_csv = DATA / "full_pool_gnina_F" / "dock_results_gnina.csv"
    if gnina_csv.exists():
        for r in load_csv(gnina_csv):
            if r.get("status") == "ok":
                gnina[r["id"]] = r
        print(f"GNINA scores: {len(gnina)}")
    else:
        print(f"GNINA scores: NOT YET (waiting on {gnina_csv})")

    # 3. QSAR scores — analogs
    qsar = {}
    for r in load_csv(DATA / "qsar" / "predictions_analogs.csv"):
        qsar[r["id"]] = r
    # 3b. Generative proposals already have QSAR scores in the same CSV
    for r in load_csv(DATA / "generative" / "generative_proposals.csv"):
        # Map column names to match
        qsar[r["id"]] = {
            "id": r["id"], "smiles": r["smiles"],
            "qsar_pkd_ens": r["qsar_pkd_ens"],
            "qsar_kd_uM_ens": r["qsar_kd_uM"],
            "tanimoto_to_parent": "",
            "max_tanimoto_to_naar": r["max_tanimoto_to_known"],
        }
    print(f"QSAR scores: {len(qsar)}")

    # 4. Boltz scores (per-compound, available for reference set + any followup)
    boltz = {}
    for r in load_csv(DATA / "boltz" / "boltz_summary.csv"):
        if r.get("status") == "ok":
            boltz[r["cid"]] = r
    print(f"Boltz scores: {len(boltz)}")

    # 5. Build the all-signals table
    rows = []
    for cid, p in pool_by_id.items():
        row = {
            "id": cid,
            "smiles": p["smiles"],
            "parent_id": p.get("parent_id", ""),
            "source": p.get("source", ""),
        }
        g = gnina.get(cid)
        if g:
            row["vina_kcal"] = g["best_vina_kcal"]
            row["cnn_pose"] = g["best_cnn_pose_score"]
            row["cnn_pkd"] = g["best_cnn_affinity_pkd"]
            row["cnn_kd_uM"] = g["best_cnn_affinity_uM"]
        q = qsar.get(cid)
        if q:
            row["qsar_pkd"] = q.get("qsar_pkd_ens", "")
            row["qsar_kd_uM"] = q.get("qsar_kd_uM_ens", "") or q.get("qsar_kd_uM", "")
            row["tanimoto_parent"] = q.get("tanimoto_to_parent", "")
            row["tanimoto_naar"] = q.get("max_tanimoto_to_naar", "") or q.get("max_tanimoto_to_known", "")
        b = boltz.get(cid)
        if b:
            row["boltz_kd_uM"] = b.get("affinity_kd_uM", "")
            row["boltz_iptm"] = b.get("ipTM", "")
            row["boltz_prob_binder"] = b.get("affinity_prob_binder", "")
        rows.append(row)

    cols = ["id", "smiles", "parent_id", "source",
            "vina_kcal", "cnn_pose", "cnn_pkd", "cnn_kd_uM",
            "qsar_pkd", "qsar_kd_uM", "tanimoto_parent", "tanimoto_naar",
            "boltz_kd_uM", "boltz_iptm", "boltz_prob_binder"]
    with open(OUT / "all_signals.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    print(f"\nWrote all_signals.csv ({len(rows)} rows)")

    # 6. Tier-A filter
    if not gnina:
        print("\nGNINA not done yet — Tier-A filtering skipped. Run again after GNINA completes.")
        return

    def f_or_none(v):
        try: return float(v)
        except (TypeError, ValueError): return None

    tier_a = []
    tier_b = []
    vina_traps = []
    for r in rows:
        vina = f_or_none(r.get("vina_kcal"))
        cnn_pose = f_or_none(r.get("cnn_pose"))
        cnn_pkd = f_or_none(r.get("cnn_pkd"))
        qsar_pkd = f_or_none(r.get("qsar_pkd"))
        if any(v is None for v in (vina, cnn_pose, cnn_pkd, qsar_pkd)):
            continue

        tier_a_pass = (
            cnn_pkd >= 4.5 and cnn_pose >= 0.5 and vina <= -6.0
            and qsar_pkd >= 4.0
        )
        tier_b_pass = (
            (cnn_pose >= 0.4 and qsar_pkd >= 3.5)
            or (cnn_pkd >= 5.0 and cnn_pose >= 0.4)
        )
        is_vina_trap = vina <= -7.0 and cnn_pose < 0.4

        if tier_a_pass:
            tier_a.append(r)
        elif tier_b_pass and not is_vina_trap:
            tier_b.append(r)
        elif is_vina_trap:
            vina_traps.append(r)

    # Composite score for Tier-A ranking: weighted sum of standardized signals
    # cnn_pkd + qsar_pkd are pKd-units; cnn_pose 0-1; we normalize and combine
    def composite(r):
        return (f_or_none(r["cnn_pkd"]) * 0.3
                + f_or_none(r["qsar_pkd"]) * 0.3
                + f_or_none(r["cnn_pose"]) * 4.0     # scale 0-1 to ~0-4
                - f_or_none(r["vina_kcal"]) * 0.2)   # negative vina is good
    tier_a.sort(key=composite, reverse=True)
    tier_b.sort(key=composite, reverse=True)

    # Save
    cols = ["id", "smiles", "parent_id", "source",
            "vina_kcal", "cnn_pose", "cnn_pkd", "cnn_kd_uM",
            "qsar_pkd", "qsar_kd_uM", "tanimoto_naar"]
    with open(OUT / "tier_a_candidates.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in tier_a:
            w.writerow({k: r.get(k, "") for k in cols})

    print(f"\n=== Tier-A passers: {len(tier_a)} ===")
    if tier_a:
        print(f"  {'id':30s}  {'cnn_pose':>8}  {'cnn_pKd':>7}  {'qsar_pKd':>8}  {'vina':>6}  parent")
        for r in tier_a[:20]:
            print(f"  {r['id']:30s}  {float(r['cnn_pose']):>8.3f}  "
                  f"{float(r['cnn_pkd']):>7.2f}  {float(r['qsar_pkd']):>8.2f}  "
                  f"{float(r['vina_kcal']):>6.2f}  {r['parent_id']}")

    print(f"\n=== Tier-B passers: {len(tier_b)} ===")
    print(f"=== Vina traps caught: {len(vina_traps)} ===")
    if vina_traps:
        for r in vina_traps[:5]:
            print(f"  ⚠ {r['id']:30s}  vina={r['vina_kcal']}  cnn_pose={r['cnn_pose']}  qsar_pkd={r['qsar_pkd']}")

    # Save Tier-B too
    with open(OUT / "tier_b_candidates.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in tier_b:
            w.writerow({k: r.get(k, "") for k in cols})


if __name__ == "__main__":
    main()
