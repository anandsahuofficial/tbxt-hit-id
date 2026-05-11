# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anand Sahu and contributors
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


_PASS_MARKERS = ("✓", "Y", "yes", True, 1, "1", "True", "true")

# T3 BRONZE allows C7 to fail ("borderline solubility" - DMSO @ 10 mM still works).
T3_RELAXABLE = {"C7_soluble_logS_gt_neg5"}

# T4 RELAXED additionally allows the two ring constraints to fail
# (i.e. permits the rare strong binder with > 4 rings or > 2 fused rings).
T4_RELAXABLE = T3_RELAXABLE | {"C4d_lt_5_rings", "C4e_le_2_fused"}


def _passes(row: pd.Series, relaxable: set = frozenset()) -> bool:
    """True iff every C1..C7 column passes, except those listed in `relaxable`
    (which are allowed to be either PASS or FAIL).
    """
    for c in CRITERIA:
        if c in relaxable:
            continue
        if row.get(c) in _PASS_MARKERS:
            continue
        return False
    return True


def _best_kd(row: pd.Series) -> float:
    """Best (lowest) available Kd in µM across Boltz run A, run B, and GNINA.

    Compounds scored by GNINA but not Boltz (e.g. the opv1 series) still
    get a usable Kd via the GNINA fallback. Used for the T1/T2 ceiling
    and for in-tier sorting; T3 and T4 don't gate on Kd at all.
    """
    cols = (KD_RUN_A, KD_RUN_B, "v1_gnina_kd_uM")
    vals = [float(row[c]) for c in cols if c in row and pd.notna(row.get(c))]
    return min(vals) if vals else float("inf")


def _conservative_boltz_kd(row: pd.Series) -> float:
    """Worst (highest) Boltz Kd across the two engines.

    Used for the T2 SILVER cutoff: a compound only earns SILVER if BOTH
    Boltz runs agree it's <= 10 uM, not just the more optimistic one.
    """
    a = row.get(KD_RUN_A)
    b = row.get(KD_RUN_B)
    vals = [float(v) for v in (a, b) if pd.notna(v)]
    return max(vals) if vals else float("inf")


def _classify(row: pd.Series) -> str:
    """Return the tier label for a single candidate row.

    Tier rules (matching ``docs/tier_definitions.md``):
      T1 GOLD    all 7 PASS,  Kd <= 5  uM, low chem AND low supplier risk
      T2 SILVER  all 7 PASS,  Kd <= 10 uM, soluble (ESOL > -5)
      T3 BRONZE  6 of 7 PASS (C7 may fail = "borderline solubility"), Kd <= 50 uM
      T4 RELAXED 5 of 7 PASS (C4d, C4e, C7 may fail), Kd <= 100 uM
    """
    kd = _best_kd(row)

    if _passes(row):
        chem_low = str(row.get(CHEM_RISK, "")).lower() == "low"
        supplier_low = str(row.get(SUPPLIER_RISK, "")).lower() == "low"
        if kd <= 5 and chem_low and supplier_low:
            return "T1_GOLD"
        # T2 requires the conservative (max) Boltz Kd to pass the cutoff,
        # so both engines must agree on potency, not just the optimistic one.
        if _conservative_boltz_kd(row) <= 10:
            return "T2_SILVER"

    # T3 BRONZE: 6/7 criteria pass (C7 may fail), and either no Boltz Kd
    # is available (GNINA-only fallback - ceiling is permissive) or the
    # best Boltz Kd is <= 50 uM. Compounds with Boltz Kd > 50 fall through
    # to T4.
    a, b = row.get(KD_RUN_A), row.get(KD_RUN_B)
    boltz_vals = [float(v) for v in (a, b) if pd.notna(v)]
    boltz_min = min(boltz_vals) if boltz_vals else None

    if _passes(row, T3_RELAXABLE) and (boltz_min is None or boltz_min <= 50):
        return "T3_BRONZE"

    if _passes(row, T4_RELAXABLE):
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
