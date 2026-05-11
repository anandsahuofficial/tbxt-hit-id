#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anand Sahu and contributors
# Package the bulk data assets into the layout that `fetch_data.sh`
# expects, write CHECKSUMS.sha256, and upload everything to a
# HuggingFace dataset repo.
#
# Required source paths:
#   --src-pool      <file>    570-compound SMILES CSV
#   --src-naar      <file>    653 Naar SPR Kd CSV
#   --src-ensemble  <dir>     directory of relaxed receptor conformations
#                              (every file inside is included; the dir is
#                              packed to receptor/tbxt_pocket_ensemble.tar.gz)
#
# Optional (with --include-poses):
#   --src-boltz   <dir>       Boltz prediction outputs root
#   --src-gnina   <dir>       GNINA dock outputs root
#   --src-mmgbsa  <dir>       MMGBSA refinement outputs root
#
# HuggingFace target (default: anandsahuofficial/tbxt-hit-id-data):
#   --hf-user <user>   override owner    (env: HF_USER)
#   --hf-repo <repo>   override repo     (env: HF_REPO)
#   --hf-token <tok>   bearer token      (env: HF_TOKEN; or run `hf auth login` once)
#
# Misc:
#   --staging <dir>     where to stage files (default: a fresh mktemp dir)
#   --include-poses     also package + upload the optional scored/ tarballs
#   --dry-run           do everything except the actual upload
#   --no-create         skip the create-repo-if-missing step
#   --help
#
# Idempotent: HF deduplicates by content; re-running with the same
# inputs uploads nothing new. SHA-256 manifest is regenerated each run.
#
# Prereqs:
#   pip install -U huggingface_hub
#   hf auth login   # one-time, opens browser
# Or set HF_TOKEN=<your-token> in the environment.
#
# Example:
#   bash setup/upload_hf_data.sh \
#       --src-pool     ~/Hackathon/TBXT/data/pool/candidate_pool_570.csv \
#       --src-naar     ~/Hackathon/TBXT/data/qsar/naar_spr_kd.csv \
#       --src-ensemble ~/Hackathon/TBXT/data/receptor/6F59_relaxed_ensemble \
#       --include-poses \
#       --src-boltz  ~/Hackathon/TBXT/data/boltz/runs \
#       --src-gnina  ~/Hackathon/TBXT/data/dock/full_pool_F \
#       --src-mmgbsa ~/Hackathon/TBXT/data/mmgbsa/top30
#
set -euo pipefail

# ─── Defaults ──────────────────────────────────────────────────────
HF_USER="${HF_USER:-anandsahuofficial}"
HF_REPO="${HF_REPO:-tbxt-hit-id-data}"
HF_TOKEN="${HF_TOKEN:-}"
STAGING=""
DRY_RUN=false
INCLUDE_POSES=false
SKIP_CREATE=false

SRC_POOL=""
SRC_NAAR=""
SRC_ENSEMBLE=""
SRC_BOLTZ=""
SRC_GNINA=""
SRC_MMGBSA=""

# ─── Args ──────────────────────────────────────────────────────────
while [ $# -gt 0 ]; do
  case "$1" in
    --src-pool)     SRC_POOL="$2"; shift 2 ;;
    --src-naar)     SRC_NAAR="$2"; shift 2 ;;
    --src-ensemble) SRC_ENSEMBLE="$2"; shift 2 ;;
    --src-boltz)    SRC_BOLTZ="$2"; shift 2 ;;
    --src-gnina)    SRC_GNINA="$2"; shift 2 ;;
    --src-mmgbsa)   SRC_MMGBSA="$2"; shift 2 ;;
    --hf-user)      HF_USER="$2"; shift 2 ;;
    --hf-repo)      HF_REPO="$2"; shift 2 ;;
    --hf-token)     HF_TOKEN="$2"; shift 2 ;;
    --staging)      STAGING="$2"; shift 2 ;;
    --include-poses) INCLUDE_POSES=true; shift ;;
    --dry-run)      DRY_RUN=true; shift ;;
    --no-create)    SKIP_CREATE=true; shift ;;
    --help|-h)      awk 'NR==1{next} /^[^#]/{exit} {print}' "$0"; exit 0 ;;
    *) echo "unknown flag: $1" >&2; exit 1 ;;
  esac
done

log() { printf "[upload_hf_data] %s\n" "$*"; }
err() { printf "[upload_hf_data] ERROR: %s\n" "$*" >&2; exit 1; }

# ─── Locate the HF CLI (only required for the real upload) ────────
HF_CLI=""
if command -v hf >/dev/null;                   then HF_CLI="hf"
elif command -v huggingface-cli >/dev/null;    then HF_CLI="huggingface-cli"
fi
if [ -n "$HF_CLI" ]; then
  log "using CLI: $HF_CLI"
elif [ "$DRY_RUN" = "true" ]; then
  log "no HF CLI on PATH (ok for --dry-run; install with 'pip install -U huggingface_hub' before real upload)"
else
  err "neither 'hf' nor 'huggingface-cli' on PATH; install with: pip install -U huggingface_hub"
fi

# Token: prefer env var, fall back to whatever 'hf auth login' cached.
if [ -n "$HF_TOKEN" ]; then
  export HF_TOKEN
  export HUGGING_FACE_HUB_TOKEN="$HF_TOKEN"
fi

# ─── Validate sources ──────────────────────────────────────────────
[ -f "$SRC_POOL" ]     || err "--src-pool not found: ${SRC_POOL:-<unset>}"
[ -f "$SRC_NAAR" ]     || err "--src-naar not found: ${SRC_NAAR:-<unset>}"
[ -d "$SRC_ENSEMBLE" ] || err "--src-ensemble not found or not a dir: ${SRC_ENSEMBLE:-<unset>}"

if [ "$INCLUDE_POSES" = "true" ]; then
  [ -d "$SRC_BOLTZ" ]  || err "--include-poses needs --src-boltz <dir>"
  [ -d "$SRC_GNINA" ]  || err "--include-poses needs --src-gnina <dir>"
  [ -d "$SRC_MMGBSA" ] || err "--include-poses needs --src-mmgbsa <dir>"
fi

# ─── Stage ─────────────────────────────────────────────────────────
[ -z "$STAGING" ] && STAGING=$(mktemp -d -t tbxt-hf-stage-XXXXXX)
mkdir -p "$STAGING"/{pool,naar,receptor}
[ "$INCLUDE_POSES" = "true" ] && mkdir -p "$STAGING/scored"

log "staging dir: $STAGING"

cp "$SRC_POOL" "$STAGING/pool/candidate_pool_570.csv"
cp "$SRC_NAAR" "$STAGING/naar/naar_spr_kd.csv"

log "packing receptor ensemble -> receptor/tbxt_pocket_ensemble.tar.gz"
ENS_PARENT=$(cd "$(dirname "$SRC_ENSEMBLE")" && pwd)
ENS_NAME=$(basename "$SRC_ENSEMBLE")
tar -czf "$STAGING/receptor/tbxt_pocket_ensemble.tar.gz" -C "$ENS_PARENT" "$ENS_NAME"

if [ "$INCLUDE_POSES" = "true" ]; then
  for pair in "boltz_outputs:$SRC_BOLTZ" "gnina_outputs:$SRC_GNINA" "mmgbsa_outputs:$SRC_MMGBSA"; do
    name="${pair%%:*}"
    src="${pair#*:}"
    log "packing $name -> scored/${name}.tar.gz"
    tar -czf "$STAGING/scored/${name}.tar.gz" -C "$(cd "$(dirname "$src")" && pwd)" "$(basename "$src")"
  done
fi

# ─── Checksums ─────────────────────────────────────────────────────
log "computing CHECKSUMS.sha256"
( cd "$STAGING" && find . -type f ! -name 'CHECKSUMS.sha256' -printf '%P\n' \
    | sort | xargs sha256sum > CHECKSUMS.sha256 )

log "staged contents:"
( cd "$STAGING" && find . -type f -printf '  %P\t' -exec stat -c '%s bytes' {} \; )

# ─── Upload ────────────────────────────────────────────────────────
if [ "$DRY_RUN" = "true" ]; then
  cat <<EOF

[dry-run] would upload contents of $STAGING to ${HF_USER}/${HF_REPO}.
Inspect the staged dir; re-run without --dry-run to upload.

EOF
  exit 0
fi

REPO_ID="${HF_USER}/${HF_REPO}"

if [ "$SKIP_CREATE" != "true" ]; then
  log "ensuring dataset repo exists: $REPO_ID"
  # `hf repo create --type dataset --exist-ok` is the modern path;
  # fall back to the older `huggingface-cli` syntax if needed.
  if [ "$HF_CLI" = "hf" ]; then
    $HF_CLI repo create "$REPO_ID" --type dataset --exist-ok 2>/dev/null \
      || $HF_CLI repo create "$REPO_ID" --type dataset 2>/dev/null \
      || log "  (repo may already exist; continuing)"
  else
    $HF_CLI repo create "$REPO_ID" --type dataset --yes 2>/dev/null \
      || log "  (repo may already exist; continuing)"
  fi
fi

log "uploading staged tree to dataset $REPO_ID"
if [ "$HF_CLI" = "hf" ]; then
  $HF_CLI upload "$REPO_ID" "$STAGING" . --repo-type dataset
else
  $HF_CLI upload "$REPO_ID" "$STAGING" . --repo-type dataset
fi

cat <<EOF

================================================================================
  ✅ Upload complete.

  Dataset: https://huggingface.co/datasets/${REPO_ID}
  Manifest: $STAGING/CHECKSUMS.sha256

  Verify the round-trip from a fresh clone:
    bash setup/fetch_data.sh
    # (or fetch_data.sh --include-poses if you uploaded the optional bundles)
================================================================================
EOF
