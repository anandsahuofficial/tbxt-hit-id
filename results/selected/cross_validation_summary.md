# Cross-Validation Summary

Every pick in [`../top4.csv`](../top4.csv) is supported by **multiple
independent lines of evidence**. This document summarizes the
cross-validation protocols and headline agreement.

## 1. Dual-engine Boltz-2 cross-validation

Two independent Boltz-2 inference runs on physically separate compute
backends, same compounds, same receptor (`6F59:A`):

- **Run A:** local single-GPU (RTX 4090), 570-compound full pool
- **Run B:** university HPC cluster (multi-A100 nodes), independent
  re-run on the same 570 compounds

**Result:** all 4 final picks agree within **1.34×** between the two
runs on predicted Kd, with the strongest pick agreeing within **1.02×**.
Per-pick values are in [`../top4.csv`](../top4.csv), columns
`boltz_kd_runA_uM` and `boltz_kd_runB_uM`.

| Pick | Kd run A (µM) | Kd run B (µM) | Ratio |
|---|---:|---:|---:|
| FM002150_analog_0083 | 3.20 | 3.26 | 1.02× |
| FM001452_analog_0104 | 3.72 | 4.97 | 1.34× |
| FM001452_analog_0201 | 8.16 | 8.76 | 1.07× |
| FM001452_analog_0171 | 8.32 | 8.17 | 1.02× |

This controls for run-to-run variance in the inference pipeline; it
does not control for systematic Boltz training-distribution bias,
which is countered by the orthogonal Vina, GNINA, QSAR, MMGBSA, and
selectivity signals (see [`../../docs/methodology.md`](../../docs/methodology.md)).

## 2. Multi-seed GNINA pose stability

A 10-seed multi-start GNINA dock at site F across the full 570-compound
pool produces a per-compound pose-stability σ for the top-pose centroid.

- σ &lt; 0.05 Å indicates a stable binding mode across stochastic restarts
- σ &gt; 0.20 Å indicates a brittle pose that the CNN happens to like

Pick #4 (`FM001452_analog_0171`) is in the most-stable subset of site-F
picks — its top pose reproduces across all 10 seeds within < 0.05 Å.

## 3. ADMET re-rank (Rowan)

The 4 picks were independently profiled on Rowan's ADMET workflow (49
properties per compound). All 4 are within standard drug-like
windows for hERG, AMES, DILI, QED, and predicted logP. ADMET output is
not redistributed in this repo (Rowan license terms); raw output is
available on request.

## 4. Onepot.ai catalog membership (muni.bio)

All 4 picks return Tanimoto similarity = 1.000 against the onepot.ai
virtual catalog when queried via the muni.bio CLI `onepot` tool, with
list price, chemistry risk, and supplier risk attached. Per-pick
catalog metadata is in the `muni_*` columns of
[`../top4.csv`](../top4.csv).

The 18-compound pre-pick subset (all 100% onepot AND all non-covalent
in the broader 570-compound pool) is in
[`onepot_100pct_non_covalent_set.csv`](onepot_100pct_non_covalent_set.csv).

## 5. Methodological variant cross-checks

Five overnight HPC variants of the pipeline were run as additional
cross-validation. Per-variant outputs are in
[`variants/`](variants/) — see that folder's README for what each
variant tested.

## What was attempted but not delivered

See [`md_attempt_log.md`](md_attempt_log.md) — explicit-solvent
pose-analysis MD was attempted on the top 4 picks but did not complete
in the available compute window. Documented honestly rather than
silently dropped.
