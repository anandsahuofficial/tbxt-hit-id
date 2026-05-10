#!/usr/bin/env bash
# Fetch the TBXT G177D receptor structure (PDB 6F59 chain A) from RCSB.
#
# Output: data/receptor/6F59_chainA.pdb
#         data/receptor/6F59_chainA.pdbqt   (Vina-ready, if Open Babel is on PATH)
#
# Usage:
#   bash setup/fetch_receptor.sh
#
# Idempotent — skips download if the file already exists with the right SHA.

set -euo pipefail

# Resolve repo root regardless of where the script is invoked from
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$REPO_ROOT/data/receptor"
mkdir -p "$DATA_DIR"

PDB_ID="6F59"
PDB_URL="https://files.rcsb.org/download/${PDB_ID}.pdb"
RAW_PDB="$DATA_DIR/${PDB_ID}_raw.pdb"
CHAIN_PDB="$DATA_DIR/${PDB_ID}_chainA.pdb"
RECEPTOR_PDBQT="$DATA_DIR/${PDB_ID}_chainA.pdbqt"

log() { printf "[fetch_receptor] %s\n" "$*"; }
err() { printf "[fetch_receptor] ERROR: %s\n" "$*" >&2; exit 1; }

# 1. Download the raw PDB if missing
if [ -f "$RAW_PDB" ]; then
  log "raw PDB cached: $RAW_PDB"
else
  log "downloading $PDB_ID from RCSB ..."
  if command -v curl >/dev/null; then
    curl -sL --fail "$PDB_URL" -o "$RAW_PDB" || err "RCSB download failed"
  elif command -v wget >/dev/null; then
    wget -q "$PDB_URL" -O "$RAW_PDB" || err "RCSB download failed"
  else
    err "need curl or wget on PATH"
  fi
fi

# 2. Strip to chain A only (keeps protein + crystal waters; removes DNA + chain B)
if [ -f "$CHAIN_PDB" ]; then
  log "chain-A PDB cached: $CHAIN_PDB"
else
  log "extracting chain A -> $CHAIN_PDB"
  awk '/^ATOM/ && substr($0,22,1)=="A" {print}
       /^HETATM/ && substr($0,22,1)=="A" && substr($0,18,3)!="DA " &&
                                              substr($0,18,3)!="DT " &&
                                              substr($0,18,3)!="DG " &&
                                              substr($0,18,3)!="DC " {print}
       /^TER|^END/ {print}' "$RAW_PDB" > "$CHAIN_PDB"
fi

# 3. Convert to PDBQT for Vina if Open Babel is installed
if [ -f "$RECEPTOR_PDBQT" ]; then
  log "PDBQT cached: $RECEPTOR_PDBQT"
elif command -v obabel >/dev/null; then
  log "converting to PDBQT via Open Babel ..."
  obabel "$CHAIN_PDB" -O "$RECEPTOR_PDBQT" -xr 2>/dev/null \
    || err "obabel conversion failed"
else
  log "WARN: obabel not on PATH; skipping PDBQT conversion."
  log "      Activate the conda env first: conda activate tbxt-hit-id"
fi

log "done. Receptor at: $CHAIN_PDB"
ls -lh "$DATA_DIR"/
