#!/usr/bin/env bash
# Smoke test - confirms the local setup is ready for a real pipeline run.
#
# Run AFTER:    setup/fetch_receptor.sh + setup/fetch_data.sh
# Run BEFORE:   examples/reproduce_top4.sh
#
# What it checks:
#   1. Repo layout - critical files + directories present
#   2. Data assets - receptor PDB, pool CSV, Naar Kd CSV, ensemble dir
#   3. Python imports - all five src/ subpackages load cleanly
#   4. tier_classify round-trip - reclassifying the committed
#      results/all_candidates_tiered.csv reproduces it exactly
#      (137/137 tier matches). This is the strongest "the code works
#      against the documented data" signal we can give you without
#      running the full ~6 h scoring pipeline.
#   5. Optional environment - NVIDIA GPU, GNINA, Vina, container .sif
#
# Exit codes:
#   0  = ready - safe to run examples/reproduce_top4.sh
#   1  = something required is missing - fix before running the pipeline
#
# Run inside the container for the most reliable results:
#   apptainer exec --bind $PWD tbxt-hit-id.sif bash setup/smoke_test.sh

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

PASS=0; FAIL=0; WARN=0
ok()   { printf "    \033[32m✓\033[0m %s\n" "$*"; PASS=$((PASS+1)); }
warn() { printf "    \033[33m!\033[0m %s\n" "$*"; WARN=$((WARN+1)); }
bad()  { printf "    \033[31m✗\033[0m %s\n" "$*"; FAIL=$((FAIL+1)); }

section() { printf "\n\033[1m[%s]\033[0m %s\n" "$1" "$2"; }

cat <<EOF
╭─────────────────────────────────────────────────────────────────╮
│  tbxt-hit-id - smoke test                                       │
│  Verifies the local setup is ready for a real pipeline run.     │
╰─────────────────────────────────────────────────────────────────╯
EOF

# ─── 1. Repo layout ────────────────────────────────────────────────
section "1/5" "Repo layout"
for f in src tools setup examples docs results slides \
         README.md LICENSE AUTHORS.md CITATION.cff environment.yml; do
  [ -e "$f" ] && ok "$f" || bad "$f missing"
done

# ─── 2. Data assets ────────────────────────────────────────────────
section "2/5" "Data assets (run setup/fetch_*.sh first if any are missing)"
for f in data/receptor/6F59_chainA.pdb \
         data/pool/candidate_pool_570.csv \
         data/naar/naar_spr_kd.csv \
         data/receptor/ensemble; do
  [ -e "$f" ] && ok "$f" || bad "$f missing"
done

if [ -d data/receptor/ensemble ]; then
  n=$(find data/receptor/ensemble -maxdepth 1 -type f | wc -l)
  [ "$n" -ge 6 ] && ok "ensemble has $n files (>= 6 expected)" \
                 || warn "ensemble has $n files (expected >= 6)"
fi

[ -f data/CHECKSUMS.sha256 ] && {
  if (cd data && sha256sum -c CHECKSUMS.sha256 >/dev/null 2>&1); then
    ok "data/CHECKSUMS.sha256 verifies"
  else
    warn "data/CHECKSUMS.sha256 mismatch (re-run setup/fetch_data.sh)"
  fi
}

# ─── 3. Python imports ─────────────────────────────────────────────
section "3/5" "Python imports"
if ! command -v python >/dev/null 2>&1; then
  bad "python not on PATH (activate the env: 'conda activate tbxt-hit-id' or run inside the container)"
else
  python - <<'PY' 2>&1
import sys, importlib
mods = ['src',
        'src.ranking.tier_classify',
        'src.filters.strict_gate',
        'src.pipeline',
        'src.enumeration',
        'src.viz']
bad = []
for m in mods:
    try:
        importlib.import_module(m)
        print(f"    \033[32m✓\033[0m {m}")
    except Exception as e:
        print(f"    \033[31m✗\033[0m {m}: {type(e).__name__}: {e}")
        bad.append(m)
sys.exit(0 if not bad else 1)
PY
  if [ $? -eq 0 ]; then PASS=$((PASS + 6)); else FAIL=$((FAIL + 1)); fi
fi

# ─── 4. tier_classify round-trip ───────────────────────────────────
section "4/5" "tier_classify round-trip on results/all_candidates_tiered.csv"
if [ ! -f results/all_candidates_tiered.csv ]; then
  bad "results/all_candidates_tiered.csv missing"
elif ! command -v python >/dev/null 2>&1; then
  warn "skipped (python not on PATH)"
else
  TMP_OUT=$(mktemp -t _smoke_tier.XXXXXX.csv)
  if python -m src.ranking.tier_classify \
        --input  results/all_candidates_tiered.csv \
        --output "$TMP_OUT" >/dev/null 2>&1; then
    python - "$TMP_OUT" <<'PY'
import sys, pandas as pd
o = pd.read_csv('results/all_candidates_tiered.csv')
n = pd.read_csv(sys.argv[1])
m = o[['id','tier']].merge(n[['id','tier']], on='id', suffixes=('_orig','_new'))
match = (m.tier_orig == m.tier_new).sum()
total = len(m)
counts = n['tier'].value_counts().to_dict()
expected = {'T1_GOLD': 0, 'T2_SILVER': 16, 'T3_BRONZE': 89, 'T4_RELAXED': 32}
all_tiers_ok = True
for tier, want in expected.items():
    got = counts.get(tier, 0)
    mark = '\033[32m✓\033[0m' if got == want else '\033[31m✗\033[0m'
    print(f"    {mark} {tier:<11} expected {want}, got {got}")
    if got != want:
        all_tiers_ok = False
mark = '\033[32m✓\033[0m' if match == total else '\033[31m✗\033[0m'
print(f"    {mark} per-row exact match: {match}/{total}")
sys.exit(0 if (match == total and all_tiers_ok) else 1)
PY
    if [ $? -eq 0 ]; then
      PASS=$((PASS + 5))
      ok "round-trip 137/137 - tier rules match the committed data"
    else
      FAIL=$((FAIL + 1))
      bad "round-trip mismatch (see above)"
    fi
  else
    bad "tier_classify failed to run - check 'python -m src.ranking.tier_classify --help'"
  fi
  rm -f "$TMP_OUT"
fi

# ─── 5. Optional environment ───────────────────────────────────────
section "5/5" "Optional environment (warnings are not blockers)"
if command -v nvidia-smi >/dev/null && nvidia-smi -L 2>/dev/null | head -1 >/dev/null; then
  GPU=$(nvidia-smi -L 2>/dev/null | head -1)
  ok "GPU: $GPU"
else
  warn "no NVIDIA GPU detected (full pipeline needs GPU; --demo mode does not)"
fi
command -v gnina >/dev/null && ok "gnina on PATH ($(command -v gnina))" || warn "gnina not on PATH (needed for the full pipeline)"
command -v vina  >/dev/null && ok "vina on PATH"  || warn "vina not on PATH"
command -v boltz >/dev/null && ok "boltz on PATH" || warn "boltz not on PATH"
ls *.sif >/dev/null 2>&1 && ok "container .sif present in $PWD" \
                        || warn "no container .sif in $PWD (use setup/pull_container.sh if you want it)"

# ─── Summary ───────────────────────────────────────────────────────
GREEN=$'\033[32m'; RED=$'\033[31m'; YELLOW=$'\033[33m'; RESET=$'\033[0m'

printf "\n"
printf "  Summary\n"
printf "  -------\n"
printf "    passes:   %d\n" "$PASS"
printf "    warnings: %d\n" "$WARN"
printf "    failures: %d\n" "$FAIL"
printf "\n"

if [ $FAIL -eq 0 ]; then
  printf "  ${GREEN}Setup looks good - safe to run:${RESET}\n"
  printf "      bash examples/reproduce_top4.sh         (full pipeline)\n"
  printf "      bash examples/reproduce_top4.sh --demo  (no GPU needed)\n\n"
  exit 0
else
  printf "  ${RED}Setup not ready${RESET} - fix the failures above before running the full pipeline.\n"
  printf "  Most-common fixes:\n"
  printf "      - missing data:    rerun setup/fetch_receptor.sh + setup/fetch_data.sh\n"
  printf "      - missing imports: activate the env or run inside the container\n\n"
  exit 1
fi
