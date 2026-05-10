#!/usr/bin/env bash
# Reproduce the top-4 picks end-to-end from a fresh clone.
#
# Two modes:
#   --full   (default)   Run the entire 6-signal scoring pipeline. Requires
#                        a CUDA GPU and ~6 h of wall time on an RTX 4090.
#   --demo               Skip scoring; use pre-computed signal CSVs from
#                        the data bundle. Reproduces strict-gate + tier
#                        classifier + final picks in < 2 minutes (no GPU).
#                        Requires `bash setup/fetch_data.sh --include-poses`.
#
# Usage:
#   bash examples/reproduce_top4.sh
#   bash examples/reproduce_top4.sh --demo
#
# Prereqs:
#   conda activate tbxt-hit-id
#   bash setup/fetch_receptor.sh
#   bash setup/fetch_data.sh         (--include-poses for --demo)

set -euo pipefail

MODE="full"
for a in "$@"; do
  case "$a" in
    --full) MODE="full" ;;
    --demo) MODE="demo" ;;
    --help) sed -n '2,18p' "$0"; exit 0 ;;
    *) echo "Unknown flag: $a" >&2; exit 1 ;;
  esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# ─── Sanity checks ───────────────────────────────────────────────
[ -f data/receptor/6F59_chainA.pdb ] \
  || { echo "ERROR: data/receptor/6F59_chainA.pdb missing - run: bash setup/fetch_receptor.sh"; exit 1; }
[ -f data/pool/candidate_pool_570.csv ] \
  || { echo "ERROR: data/pool/candidate_pool_570.csv missing - run: bash setup/fetch_data.sh"; exit 1; }
[ -f data/naar/naar_spr_kd.csv ] \
  || { echo "ERROR: data/naar/naar_spr_kd.csv missing - run: bash setup/fetch_data.sh"; exit 1; }

mkdir -p data/scored results

# ─── Stage 1: receptor prep + pocket definition ──────────────────
log() { printf "\n[reproduce-top4 %s] %s\n" "$(date +%H:%M:%S)" "$*"; }

log "Stage 1: receptor prep + pocket definition"
python -m src.pipeline.prep_receptor   --in data/receptor/6F59_chainA.pdb \
                                        --out data/receptor/6F59_chainA.pdbqt
python -m src.pipeline.define_pockets   --receptor data/receptor/6F59_chainA.pdbqt \
                                        --site F --site A \
                                        --out data/pockets.json

if [ "$MODE" = "demo" ]; then
  # ─── Demo path: skip scoring, use pre-computed CSVs ───────────
  log "Demo mode: using pre-computed signal CSVs from data/scored/"
  for f in vina gnina gnina_multiseed boltz mmgbsa qsar paralog; do
    [ -f data/scored/${f}_scores.csv ] \
      || { echo "ERROR: data/scored/${f}_scores.csv missing - run: bash setup/fetch_data.sh --include-poses"; exit 1; }
  done
else
  # ─── Stage 2: scoring (full mode) ─────────────────────────────
  log "Stage 2a: Vina ensemble (570 cmpds × 6 receptor confs, ~45 min)"
  python -m src.pipeline.dock_vina       --pool data/pool/candidate_pool_570.csv \
                                          --pockets data/pockets.json \
                                          --resume \
                                          --out data/scored/vina_scores.csv

  log "Stage 2b: GNINA CNN re-rank (~30 min)"
  python -m src.pipeline.dock_gnina      --pool data/pool/candidate_pool_570.csv \
                                          --vina-poses data/scored/vina_poses.sdf \
                                          --resume \
                                          --out data/scored/gnina_scores.csv

  log "Stage 2c: GNINA 10-seed pose stability (~5 h)"
  python -m src.pipeline.dock_gnina_multiseed --pool data/pool/candidate_pool_570.csv \
                                          --seeds 10 --resume \
                                          --out data/scored/gnina_multiseed.csv

  log "Stage 2d: TBXT QSAR (~2 min)"
  python -m src.pipeline.train_qsar      --naar data/naar/naar_spr_kd.csv \
                                          --predict data/pool/candidate_pool_570.csv \
                                          --out data/scored/qsar_scores.csv

  log "Stage 2e: Boltz-2 co-folding (~3 h)"
  python -m src.pipeline.run_boltz       --pool data/pool/candidate_pool_570.csv \
                                          --receptor 6F59 --resume \
                                          --out data/scored/boltz_scores.csv

  log "Stage 2f: MMGBSA on top 30 (~1.5 h)"
  python -m src.pipeline.run_mmgbsa      --pool data/pool/candidate_pool_570.csv \
                                          --top-n 30 \
                                          --out data/scored/mmgbsa_scores.csv

  log "Stage 2g: Paralog selectivity (~5 min)"
  python -m src.pipeline.paralog_selectivity --pool data/pool/candidate_pool_570.csv \
                                          --out data/scored/paralog_scores.csv
fi

# ─── Stage 3: signal merge + consensus ───────────────────────────
log "Stage 3: signal merge + consensus aggregation"
python -m src.pipeline.merge_signals  --vina  data/scored/vina_scores.csv \
                                       --gnina data/scored/gnina_scores.csv \
                                       --qsar  data/scored/qsar_scores.csv \
                                       --boltz data/scored/boltz_scores.csv \
                                       --mmgbsa data/scored/mmgbsa_scores.csv \
                                       --paralog data/scored/paralog_scores.csv \
                                       --out data/scored/all_signals_merged.csv
python -m src.pipeline.consensus      --merged data/scored/all_signals_merged.csv \
                                       --out data/scored/consensus_scored.csv

# ─── Stage 4: 7-criterion strict gate ───────────────────────────
log "Stage 4: 7-criterion strict gate"
python -m src.filters.strict_gate     --input data/scored/consensus_scored.csv \
                                       --naar  data/naar/naar_spr_kd.csv \
                                       --output data/scored/all_with_flags.csv

# ─── Stage 5: tier classification + sort ────────────────────────
log "Stage 5: tier classification"
python -m src.ranking.tier_classify   --input data/scored/all_with_flags.csv \
                                       --output results/all_candidates_tiered.csv

# Top 4 = the first 4 rows of the tiered output (header + 4 data rows)
head -1 results/all_candidates_tiered.csv > results/top4.csv
awk 'NR>1 && NR<=5' results/all_candidates_tiered.csv >> results/top4.csv
# Top 5-24 = the next 20
head -1 results/all_candidates_tiered.csv > results/top5to24.csv
awk 'NR>=6 && NR<=25' results/all_candidates_tiered.csv >> results/top5to24.csv

# ─── Done ────────────────────────────────────────────────────────
log "Done."
cat <<EOF

================================================================================
  ✅ End-to-end reproduction complete (mode: ${MODE}).

  Outputs:
    results/top4.csv                     ← 4 picks for judging
    results/top5to24.csv                 ← 20 additional for the experimental program
    results/all_candidates_tiered.csv    ← full strict-pass set (137 rows expected)

  Compare against the curated outputs already in results/ to verify that
  your local run matches the published numbers.
================================================================================
EOF
