#!/usr/bin/env bash
# Pull the all-batteries-included tbxt-hit-id container from GHCR.
#
# Output: ./tbxt-hit-id.sif (Apptainer SIF file, ~6-12 GB)
#
# Usage:
#   bash setup/pull_container.sh                 # default: latest
#   bash setup/pull_container.sh --tag v1.0.0    # pinned version
#   bash setup/pull_container.sh --image-ref docker://ghcr.io/<other>/<repo>:tag
#   bash setup/pull_container.sh --out custom.sif

set -euo pipefail

OWNER="${GHCR_OWNER:-anandsahuofficial}"
REPO="${GHCR_REPO:-tbxt-hit-id}"
TAG="latest"
OUT="tbxt-hit-id.sif"
IMAGE_REF=""

# ─── Args ────────────────────────────────────────────────────────
while [ $# -gt 0 ]; do
  case "$1" in
    --tag)        TAG="$2"; shift 2 ;;
    --out)        OUT="$2"; shift 2 ;;
    --image-ref)  IMAGE_REF="$2"; shift 2 ;;
    --help|-h)    sed -n '2,12p' "$0"; exit 0 ;;
    *) echo "Unknown flag: $1" >&2; exit 1 ;;
  esac
done

[ -z "$IMAGE_REF" ] && IMAGE_REF="docker://ghcr.io/${OWNER}/${REPO}:${TAG}"

# ─── Apptainer / Singularity discovery ───────────────────────────
RUNNER=""
for cmd in apptainer singularity; do
  if command -v "$cmd" >/dev/null; then RUNNER="$cmd"; break; fi
done
if [ -z "$RUNNER" ]; then
  cat >&2 <<EOF
ERROR: neither apptainer nor singularity is on PATH.

Install one of:
  - Apptainer (recommended):  https://apptainer.org/docs/admin/main/installation.html
  - SingularityCE:            https://docs.sylabs.io/guides/latest/admin-guide/installation.html
EOF
  exit 1
fi

# ─── Pull ────────────────────────────────────────────────────────
echo "[pull_container] runner:    $RUNNER"
echo "[pull_container] source:    $IMAGE_REF"
echo "[pull_container] target:    $OUT"
echo "[pull_container] (this is a multi-GB pull; first run takes 5-15 min)"

$RUNNER pull --force "$OUT" "$IMAGE_REF"

# ─── Smoke test ──────────────────────────────────────────────────
echo
echo "[pull_container] smoke test: print env name + python version"
$RUNNER exec "$OUT" bash -lc 'echo "  env:    ${CONDA_DEFAULT_ENV:-?}"; echo "  python: $(python --version 2>&1)"; echo "  gnina:  $(command -v gnina || echo MISSING)"; echo "  rdkit:  $(python -c "import rdkit; print(rdkit.__version__)" 2>&1)"'

cat <<EOF

================================================================================
  ✅ Container ready at:  $OUT

  Run the demo path (no GPU needed if --include-poses cached):
    $RUNNER exec --bind \$PWD $OUT bash examples/reproduce_top4.sh --demo

  Run the full pipeline (needs GPU):
    $RUNNER exec --nv --bind \$PWD $OUT bash examples/reproduce_top4.sh

  Drop into the container shell for ad-hoc work:
    $RUNNER shell --nv --bind \$PWD $OUT
================================================================================
EOF
