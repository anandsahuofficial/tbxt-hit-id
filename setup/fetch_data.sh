#!/usr/bin/env bash
# Fetch the bulk data assets that are too large for git but required to
# run the pipeline end-to-end:
#   - 570-compound novelty-filtered candidate pool (SMILES CSV)
#   - 653 measured Naar SPR Kd dataset (QSAR training data)
#   - 6-conformation receptor ensemble (relaxed from 6F59:A via short MD)
#   - (optional) pre-computed Boltz / GNINA / MMGBSA outputs for skipping
#     expensive scoring steps
#
# Source: a Hugging Face dataset repo (override via HF_USER / HF_REPO env).
#
# Usage:
#   bash setup/fetch_data.sh                    # default: anandsahuofficial/tbxt-hit-id-data
#   bash setup/fetch_data.sh --include-poses    # also pull pre-computed poses (~600 MB)
#   HF_USER=foo HF_REPO=bar bash setup/fetch_data.sh
#   HF_TOKEN=hf_xxx bash setup/fetch_data.sh    # for private datasets
#
# Idempotent: re-running skips files whose SHA matches CHECKSUMS.sha256.

set -euo pipefail

# ─── Args ────────────────────────────────────────────────────────────
INCLUDE_POSES="false"
for a in "$@"; do
  case "$a" in
    --include-poses) INCLUDE_POSES="true" ;;
    --help) echo "Usage: $0 [--include-poses]"; exit 0 ;;
    *) echo "Unknown flag: $a" >&2; exit 1 ;;
  esac
done

# ─── Resolve paths ──────────────────────────────────────────────────
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$REPO_ROOT/data"
mkdir -p "$DATA_DIR"/{pool,naar,receptor,scored}

# ─── HF source config ───────────────────────────────────────────────
HF_USER="${HF_USER:-anandsahuofficial}"
HF_REPO="${HF_REPO:-tbxt-hit-id-data}"
HF_BRANCH="${HF_BRANCH:-main}"
HF_TOKEN="${HF_TOKEN:-}"
HF_BASE="https://huggingface.co/datasets/${HF_USER}/${HF_REPO}/resolve/${HF_BRANCH}"

log() { printf "[fetch_data] %s\n" "$*"; }
err() { printf "[fetch_data] ERROR: %s\n" "$*" >&2; exit 1; }

# ─── Download helper ────────────────────────────────────────────────
hf_dl() {
  local filename="$1" out="$2" expected_sha="${3:-}"
  if [ -f "$out" ] && [ -n "$expected_sha" ]; then
    local actual; actual=$(sha256sum "$out" | awk '{print $1}')
    [ "$actual" = "$expected_sha" ] && { log "  cached + verified: ${filename}"; return 0; }
  fi
  log "  downloading: ${filename}"
  local hdr=()
  [ -n "$HF_TOKEN" ] && hdr=(-H "Authorization: Bearer ${HF_TOKEN}")
  if command -v curl >/dev/null; then
    curl -sL --fail -C - -o "$out" "${HF_BASE}/${filename}" "${hdr[@]}" \
      || err "download failed: ${filename}"
  elif command -v wget >/dev/null; then
    local wgt=(); [ -n "$HF_TOKEN" ] && wgt=(--header="Authorization: Bearer ${HF_TOKEN}")
    wget -q --continue -O "$out" "${HF_BASE}/${filename}" "${wgt[@]}" \
      || err "download failed: ${filename}"
  else
    err "need curl or wget on PATH"
  fi
}

# ─── Fetch CHECKSUMS first so we can verify all subsequent files ────
log "Fetching from Hugging Face dataset: ${HF_USER}/${HF_REPO}"
if ! hf_dl "CHECKSUMS.sha256" "$DATA_DIR/CHECKSUMS.sha256" 2>/dev/null; then
  cat >&2 <<EOF

[fetch_data] ERROR: could not download CHECKSUMS.sha256 from
            https://huggingface.co/datasets/${HF_USER}/${HF_REPO}

Most likely causes:
  1. The dataset has not been published yet at that path. The bulk
     data bundle (570-cmpd pool, Naar SPR Kd, receptor ensemble) is
     uploaded separately. Until it is published you can still:
       - inspect the curated post-pipeline outputs in results/
       - read the methodology in docs/
       - render the slide deck via tools/render_slides.py
  2. The dataset is private and HF_TOKEN is not set:
       HF_TOKEN=hf_xxx bash setup/fetch_data.sh
  3. You are pointing at the wrong dataset. Override with:
       HF_USER=<user> HF_REPO=<repo> bash setup/fetch_data.sh
EOF
  exit 1
fi

sha_for() {
  grep -E "[[:space:]]+${1}\$" "$DATA_DIR/CHECKSUMS.sha256" | awk '{print $1}'
}

# ─── Required assets ────────────────────────────────────────────────
hf_dl "pool/candidate_pool_570.csv"          "$DATA_DIR/pool/candidate_pool_570.csv"          "$(sha_for pool/candidate_pool_570.csv)"
hf_dl "naar/naar_spr_kd.csv"             "$DATA_DIR/naar/naar_spr_kd.csv"             "$(sha_for naar/naar_spr_kd.csv)"
hf_dl "receptor/tbxt_pocket_ensemble.tar.gz" "$DATA_DIR/receptor/tbxt_pocket_ensemble.tar.gz" "$(sha_for receptor/tbxt_pocket_ensemble.tar.gz)"

# Unpack receptor ensemble in place
if [ ! -d "$DATA_DIR/receptor/ensemble" ]; then
  log "Extracting receptor ensemble ..."
  tar -xzf "$DATA_DIR/receptor/tbxt_pocket_ensemble.tar.gz" -C "$DATA_DIR/receptor"
fi

# ─── Optional: pre-computed pose / score outputs ─────────────────────
if [ "$INCLUDE_POSES" = "true" ]; then
  log "Including pre-computed scoring outputs (~600 MB) ..."
  hf_dl "scored/boltz_outputs.tar.gz" "$DATA_DIR/scored/boltz_outputs.tar.gz" "$(sha_for scored/boltz_outputs.tar.gz)"
  hf_dl "scored/gnina_outputs.tar.gz" "$DATA_DIR/scored/gnina_outputs.tar.gz" "$(sha_for scored/gnina_outputs.tar.gz)"
  hf_dl "scored/mmgbsa_outputs.tar.gz" "$DATA_DIR/scored/mmgbsa_outputs.tar.gz" "$(sha_for scored/mmgbsa_outputs.tar.gz)"
  for f in boltz_outputs gnina_outputs mmgbsa_outputs; do
    [ -d "$DATA_DIR/scored/$f" ] && continue
    log "  unpacking $f ..."
    tar -xzf "$DATA_DIR/scored/${f}.tar.gz" -C "$DATA_DIR/scored"
  done
fi

# ─── GPU torch upgrade if available ─────────────────────────────────
if command -v nvidia-smi >/dev/null && nvidia-smi >/dev/null 2>&1; then
  if python -c "import torch; assert torch.cuda.is_available()" 2>/dev/null; then
    log "torch CUDA already available - skipping torch upgrade"
  else
    log "Upgrading torch to CUDA 12.8 wheel (NVIDIA GPU detected) ..."
    pip install --quiet --force-reinstall --no-deps \
      "torch==2.8.0" "torchvision==0.23.0" \
      --index-url https://download.pytorch.org/whl/cu128
  fi
fi

# ─── Summary ─────────────────────────────────────────────────────────
cat <<EOF

================================================================================
  ✅ Data fetch complete.

  Pool          : $(wc -l < "$DATA_DIR/pool/candidate_pool_570.csv") lines  (data/pool/candidate_pool_570.csv)
  Naar SPR Kd   : $(wc -l < "$DATA_DIR/naar/naar_spr_kd.csv") lines  (data/naar/naar_spr_kd.csv)
  Receptor ens. : $(ls "$DATA_DIR/receptor/ensemble" 2>/dev/null | wc -l) confs  (data/receptor/ensemble/)
EOF

if [ "$INCLUDE_POSES" = "true" ]; then
  echo "  Pre-scored    : data/scored/{boltz,gnina,mmgbsa}_outputs/"
fi

cat <<EOF

  Next:
    bash examples/reproduce_top4.sh    # runs the full pipeline end-to-end
================================================================================
EOF
