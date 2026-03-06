#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

FAMILIES=(
  ACT
  UBQ
  EF1A
  PRK
)

echo "[INFO] Repo root: $REPO_ROOT"
echo "[INFO] Families to run: ${FAMILIES[*]}"

for fam in "${FAMILIES[@]}"; do
  echo
  echo "=============================="
  echo "[INFO] Running family: $fam"
  echo "=============================="
  bash scripts/run_family.sh "$fam"
done

echo
echo "[DONE] All families completed."
