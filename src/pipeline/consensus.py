"""
Aggregate signals from all strategies into a single consensus Tier-A list.

Inputs (all optional — uses what's available):
  --multiseed-csv        : data/full_pool_gnina_F_multiseed/dock_results_multiseed.csv
  --site-a-csv           : data/full_pool_gnina_A/dock_results_gnina.csv
  --boltz-csv            : data/boltz/full_pool_summary.csv
  --mmgbsa-csv           : data/mmgbsa/top50_results.csv
  --selectivity-csv      : data/selectivity/dock_offtarget.csv
  --fep-csv              : data/fep/relative_dg_table.csv
  --qsar-csv             : data/qsar/predictions_analogs.csv
  --pool                 : data/full_pool_input.csv (master pool)

Output:
  --out-csv              : data/tier_a/final_tier_a.csv  (sorted by composite score)
"""
import argparse
import csv
from pathlib import Path


def f_or_none(v):
    if v is None or v == "": return None
    try: return float(v)
    except (TypeError, ValueError): return None


def normalize(v, lo, hi, clip=True):
    """Map v in [lo, hi] → [0, 1]. Outside range clamped if clip=True."""
    if v is None: return 0.0
    n = (v - lo) / (hi - lo)
    if clip: n = max(0.0, min(1.0, n))
    return n


def load_csv(path):
    if not path or not Path(path).exists(): return {}
    return {r.get("id") or r.get("cid"): r for r in csv.DictReader(open(path))}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--multiseed-csv", default="data/full_pool_gnina_F_multiseed/dock_results_multiseed.csv")
    p.add_argument("--site-a-csv", default="data/full_pool_gnina_A/dock_results_gnina.csv")
    p.add_argument("--boltz-csv", default="data/boltz/full_pool_summary.csv")
    p.add_argument("--mmgbsa-csv", default="data/mmgbsa/top50_results.csv")
    p.add_argument("--selectivity-csv", default="data/selectivity/dock_offtarget.csv")
    p.add_argument("--fep-csv", default="data/fep/relative_dg_table.csv")
    p.add_argument("--qsar-csv", default="data/qsar/predictions_analogs.csv")
    p.add_argument("--pool", default="data/full_pool_input.csv")
    p.add_argument("--out-csv", default="data/tier_a/final_tier_a.csv")
    args = p.parse_args()

    # Load all signal sources
    pool = load_csv(args.pool)
    multiseed = load_csv(args.multiseed_csv)
    site_a = load_csv(args.site_a_csv)
    boltz = load_csv(args.boltz_csv)
    mmgbsa = load_csv(args.mmgbsa_csv)
    selectivity = load_csv(args.selectivity_csv)
    qsar = load_csv(args.qsar_csv)

    print(f"Pool: {len(pool)}  multiseed: {len(multiseed)}  site_a: {len(site_a)}  "
          f"boltz: {len(boltz)}  mmgbsa: {len(mmgbsa)}  qsar: {len(qsar)}")

    out = []
    for cid, p_row in pool.items():
        ms = multiseed.get(cid, {})
        sa = site_a.get(cid, {})
        bz = boltz.get(cid, {})
        mg = mmgbsa.get(cid, {})
        sel = selectivity.get(cid, {})
        q = qsar.get(cid, {})

        cnn_pose = f_or_none(ms.get("cnn_pose_mean")) if ms else f_or_none(ms.get("best_cnn_pose_score"))
        cnn_pose_stdev = f_or_none(ms.get("cnn_pose_stdev")) if ms else None
        cnn_pkd = f_or_none(ms.get("cnn_pkd_mean")) if ms else f_or_none(ms.get("best_cnn_affinity_pkd"))
        vina = f_or_none(ms.get("vina_kcal_mean")) if ms else f_or_none(ms.get("best_vina_kcal"))
        qsar_pkd = f_or_none(q.get("qsar_pkd_ens")) or f_or_none(q.get("qsar_pkd"))
        boltz_pkd = f_or_none(bz.get("affinity_pkd"))
        boltz_prob = f_or_none(bz.get("affinity_prob_binder"))
        boltz_iptm = f_or_none(bz.get("ipTM"))
        mmgbsa_dg = f_or_none(mg.get("delta_g_bind_kcal"))

        # Hard requirements (Tier-A)
        passes_hard = (
            cnn_pose is not None and cnn_pose >= 0.5
            and cnn_pkd is not None and cnn_pkd >= 4.5
            and vina is not None and vina <= -6.0
            and qsar_pkd is not None and qsar_pkd >= 4.0
        )

        # Composite score (only computed if cnn_pose / qsar_pkd present)
        composite = None
        if cnn_pose is not None and qsar_pkd is not None:
            composite = (
                0.25 * cnn_pose +
                0.20 * normalize(qsar_pkd, 4.0, 5.5) +
                0.15 * normalize(cnn_pkd or 4.5, 4.5, 6.5) +
                0.15 * normalize(boltz_prob or 0.3, 0.3, 0.7) +
                0.15 * normalize(-(mmgbsa_dg or -10), 10, 25) +
                0.10 * 0.5  # FEP placeholder
            )

        out.append({
            "id": cid,
            "smiles": p_row.get("smiles"),
            "parent_id": p_row.get("parent_id"),
            "source": p_row.get("source"),
            "cnn_pose": cnn_pose,
            "cnn_pose_stdev": cnn_pose_stdev,
            "cnn_pkd": cnn_pkd,
            "vina_kcal": vina,
            "qsar_pkd": qsar_pkd,
            "boltz_pkd": boltz_pkd,
            "boltz_prob_binder": boltz_prob,
            "boltz_iptm": boltz_iptm,
            "mmgbsa_dg_kcal": mmgbsa_dg,
            "passes_hard": passes_hard,
            "composite_score": round(composite, 4) if composite is not None else None,
            "tanimoto_naar": f_or_none(q.get("max_tanimoto_to_naar")) or f_or_none(q.get("max_tanimoto_to_known")),
        })

    # Sort: hard-passers first by composite descending, then non-hard-passers
    out.sort(key=lambda r: (
        not r["passes_hard"],
        -(r["composite_score"] or 0)
    ))

    # Write
    cols = ["id", "smiles", "parent_id", "source",
            "cnn_pose", "cnn_pose_stdev", "cnn_pkd", "vina_kcal",
            "qsar_pkd", "boltz_pkd", "boltz_prob_binder", "boltz_iptm",
            "mmgbsa_dg_kcal", "tanimoto_naar",
            "passes_hard", "composite_score"]

    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out:
            w.writerow({k: r.get(k, "") for k in cols})

    n_hard = sum(1 for r in out if r["passes_hard"])
    print(f"\n  Wrote {args.out_csv}")
    print(f"  Total: {len(out)}; hard-requirement passers: {n_hard}")
    print(f"\n  Top 10 hard-passers by composite:")
    print(f"    {'id':30s}  {'composite':>10}  {'cnn_pose':>9}  {'qsar':>5}  {'boltz_uM':>9}")
    n = 0
    for r in out:
        if not r["passes_hard"]: continue
        n += 1
        if n > 10: break
        boltz_uM = (10**(-(r["boltz_pkd"] or 0))*1e6) if r["boltz_pkd"] else None
        boltz_str = f"{boltz_uM:>9.2f}" if boltz_uM else "        -"
        print(f"    {r['id']:30s}  {r['composite_score'] or 0:>10.4f}  "
              f"{r['cnn_pose'] or 0:>9.3f}  {r['qsar_pkd'] or 0:>5.2f}  {boltz_str}")


if __name__ == "__main__":
    main()
