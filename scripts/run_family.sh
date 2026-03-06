#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_family.sh FAMILY

Examples:
  bash scripts/run_family.sh ACT
  bash scripts/run_family.sh UBQ
  bash scripts/run_family.sh EF1A
EOF
}

if [ $# -ne 1 ]; then
  usage
  exit 1
fi

FAMILY="$1"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

RUN_SCRIPT="scripts/run_${FAMILY}.sh"
POST_V2_SCRIPT="scripts/post_${FAMILY}_standardize.v2.sh"
POST_SCRIPT="scripts/post_${FAMILY}_standardize.sh"

echo "[INFO] Repo root: $REPO_ROOT"
echo "[INFO] Family: $FAMILY"

# 0) sanitize inputs first
if [ -f scripts/00_sanitize_inputs.sh ]; then
  echo "[INFO] Running sanitize step"
  bash scripts/00_sanitize_inputs.sh
else
  echo "[WARN] scripts/00_sanitize_inputs.sh not found, skip sanitize"
fi

# 1) run family pipeline
if [ ! -f "$RUN_SCRIPT" ]; then
  echo "[ERROR] Missing run script: $RUN_SCRIPT"
  exit 1
fi

echo "[INFO] Running: $RUN_SCRIPT"
bash "$RUN_SCRIPT"

# 2) post standardize
if [ -f "$POST_V2_SCRIPT" ]; then
  echo "[INFO] Running: $POST_V2_SCRIPT"
  bash "$POST_V2_SCRIPT"
elif [ -f "$POST_SCRIPT" ]; then
  echo "[INFO] Running: $POST_SCRIPT"
  bash "$POST_SCRIPT"
else
  echo "[WARN] No post standardize script found for $FAMILY"
fi

echo "[DONE] Family completed: $FAMILY"
