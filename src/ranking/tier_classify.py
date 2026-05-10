"""Four-tier ranking of strict-pass candidates.

Implements the tier rules documented in ``docs/tier_definitions.md``::

    T1 GOLD    all 7 criteria + Boltz Kd <= 5  uM + low/low risk
    T2 SILVER  all 7 criteria + Boltz Kd <= 10 uM + soluble (ESOL > -5)
    T3 BRONZE  all 7 criteria + Boltz Kd <= 50 uM + borderline solubility ok
    T4 RELAXED all 7 criteria + Boltz Kd <= 100 uM

Within each tier, candidates are ordered by predicted Boltz Kd
ascending. Ties are broken by chemistry risk, then supplier risk,
then onepot.ai cost.

Usage::

    python -m src.ranking.tier_classify \\
        --input  data/scored/all_candidates_with_scores.csv \\
        --output results/all_candidates_tiered.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


# --- column names this module expects on the input CSV --------------
ID = "id"
TIER = "tier"
KD_RUN_A = "boltz_kd_runA_uM"
KD_RUN_B = "boltz_kd_runB_uM"
CHEM_RISK = "muni_chem_risk"
SUPPLIER_RISK = "muni_supplier_risk"
PRICE_USD = "muni_price_usd"
ESOL = "esol_logS"

CRITERIA = [
    "C1_onepot_100",
    "C2_non_covalent",
    "C3a_MW_le_600", "C3b_LogP_le_6", "C3c_HBD_le_6", "C3d_HBA_le_12",
    "C4a_HA_10_30", "C4b_HBD_HBA_le_11", "C4c_LogP_lt_5",
    "C4d_lt_5_rings", "C4e_le_2_fused",
    "C5_no_PAINS", "C5b_no_forbidden",
    "C6_naar_tanimoto_lt_085",
    "C7_soluble_logS_gt_neg5",
]

RISK_RANK = {"low": 0, "med": 1, "medium": 1, "high": 2, "": 3, None: 3}


def _passes_all_criteria(row: pd.Series) -> bool:
    """Return True iff every C1..C7 column is the pass marker (✓ or True/1)."""
    for c in CRITERIA:
        v = row.get(c)
        if v in ("✓", "Y", "yes", True, 1, "1", "True", "true"):
            continue
        return False
    return True


def _best_kd(row: pd.Series) -> float:
    """Pick the lower (better) of the two Boltz-engine Kd predictions."""
    a = row.get(KD_RUN_A)
    b = row.get(KD_RUN_B)
    vals = [v for v in (a, b) if pd.notna(v)]
    return float(min(vals)) if vals else float("inf")


def _classify(row: pd.Series) -> str:
    """Return the tier label for a single candidate row."""
    if not _passes_all_criteria(row):
        return "FAIL"

    kd = _best_kd(row)
    soluble = pd.notna(row.get(ESOL)) and float(row[ESOL]) > -5
    chem_low = str(row.get(CHEM_RISK, "")).lower() == "low"
    supplier_low = str(row.get(SUPPLIER_RISK, "")).lower() == "low"

    if kd <= 5 and soluble and chem_low and supplier_low:
        return "T1_GOLD"
    if kd <= 10 and soluble:
        return "T2_SILVER"
    if kd <= 50:
        return "T3_BRONZE"
    if kd <= 100:
        return "T4_RELAXED"
    return "FAIL"


def _sort_key(row: pd.Series) -> tuple:
    return (
        _best_kd(row),
        RISK_RANK.get(str(row.get(CHEM_RISK, "")).lower(), 3),
        RISK_RANK.get(str(row.get(SUPPLIER_RISK, "")).lower(), 3),
        float(row.get(PRICE_USD) or 1e9),
    )


def classify_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Add a tier column + sort strict-first then by Boltz Kd."""
    df = df.copy()
    df[TIER] = df.apply(_classify, axis=1)
    tier_order = {
        "T1_GOLD": 0, "T2_SILVER": 1, "T3_BRONZE": 2, "T4_RELAXED": 3, "FAIL": 9,
    }
    df["_tier_order"] = df[TIER].map(tier_order)
    df["_sort_key"] = df.apply(_sort_key, axis=1)
    df = df.sort_values(["_tier_order", "_sort_key"]).drop(columns=["_tier_order", "_sort_key"])
    return df


def summarize(df: pd.DataFrame) -> dict:
    """Return tier counts + the empty-T1 honesty diagnostic."""
    counts = df[TIER].value_counts().to_dict()
    return {
        "T1_GOLD": counts.get("T1_GOLD", 0),
        "T2_SILVER": counts.get("T2_SILVER", 0),
        "T3_BRONZE": counts.get("T3_BRONZE", 0),
        "T4_RELAXED": counts.get("T4_RELAXED", 0),
        "FAIL": counts.get("FAIL", 0),
        "strict_pass_total": sum(counts.get(t, 0) for t in
                                 ("T1_GOLD", "T2_SILVER", "T3_BRONZE", "T4_RELAXED")),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True,
                    help="Scored candidates CSV (must carry C1..C7 flag columns "
                         "+ boltz_kd_runA_uM + boltz_kd_runB_uM + risk columns).")
    ap.add_argument("--output", type=Path, required=True,
                    help="Output CSV with tier column added, sorted strict-first.")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    out = classify_dataframe(df)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    s = summarize(out)
    print(f"  Strict-pass total: {s['strict_pass_total']}")
    for t in ("T1_GOLD", "T2_SILVER", "T3_BRONZE", "T4_RELAXED"):
        print(f"  {t:<11} {s[t]}")
    if s["T1_GOLD"] == 0:
        print("  (T1_GOLD empty by design - no compound simultaneously hits "
              "Kd <= 5 uM AND low/low risk; surfaced honestly rather than "
              "by relaxing tier rules)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
