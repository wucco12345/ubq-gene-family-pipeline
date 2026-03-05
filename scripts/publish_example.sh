#!/usr/bin/env bash
set -euo pipefail

# Publish small, standardized outputs from results/<FAMILY>/ to example/<FAMILY>/,
# then git add/commit (and optionally push).
#
# Usage:
#   bash scripts/publish_example.sh EF1A
#   bash scripts/publish_example.sh UBQ --push
#
# Assumptions:
# - Your repo root contains: results/<FAMILY>/ and example/
# - "Small files" include: *_high_conf.{fa,id}, *_rescue.{fa,id}, *_qc.txt
# - Evidence: *_blast_top1.tsv (if exists), plus some UBQ evidence IDs if present.

FAMILY="${1:-}"
DO_PUSH="${2:-}"

if [[ -z "${FAMILY}" ]]; then
  echo "Usage: bash scripts/publish_example.sh <FAMILY> [--push]"
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RES_DIR="${REPO_ROOT}/results/${FAMILY}"
EX_DIR="${REPO_ROOT}/example/${FAMILY}"

if [[ ! -d "${RES_DIR}" ]]; then
  echo "[ERROR] results directory not found: ${RES_DIR}"
  exit 2
fi

mkdir -p "${EX_DIR}"

echo "[INFO] Publishing ${FAMILY}: ${RES_DIR} -> ${EX_DIR}"

# Copy standardized "small outputs" (if present)
shopt -s nullglob

# common
for f in "${RES_DIR}/${FAMILY}_high_conf.fa" "${RES_DIR}/${FAMILY}_high_conf.id" \
         "${RES_DIR}/${FAMILY}_rescue.fa"    "${RES_DIR}/${FAMILY}_rescue.id" \
         "${RES_DIR}/${FAMILY}_qc.txt"       "${RES_DIR}/${FAMILY}_excluded.tsv" \
         "${RES_DIR}/${FAMILY}_nohit.id"
do
  if [[ -f "$f" ]]; then
    cp -f "$f" "${EX_DIR}/"
  fi
done

# evidence (optional)
if [[ -d "${RES_DIR}/evidence" ]]; then
  # copy top1 if exists
  for f in "${RES_DIR}/evidence/${FAMILY}_blast_top1.tsv" \
           "${RES_DIR}/evidence/${FAMILY}_vs_sprot.tsv" \
           "${RES_DIR}/evidence/${FAMILY}_pfam_hits.tsv"
  do
    if [[ -f "$f" ]]; then
      cp -f "$f" "${EX_DIR}/"
    fi
  done

  # UBQ-specific small evidence IDs (if exist)
  for f in "${RES_DIR}/evidence/${FAMILY}_keep_strict.id" \
           "${RES_DIR}/evidence/${FAMILY}_poly_merged.id"
  do
    if [[ -f "$f" ]]; then
      cp -f "$f" "${EX_DIR}/"
    fi
  done
fi

shopt -u nullglob

# Show what was published
echo "[INFO] Published files:"
ls -lh "${EX_DIR}"

# Git add + commit
cd "${REPO_ROOT}"
git add "example/${FAMILY}"

# If nothing changed, exit cleanly
if git diff --cached --quiet; then
  echo "[INFO] No changes to commit for example/${FAMILY}."
  exit 0
fi

git commit -m "Publish ${FAMILY} example outputs"
echo "[OK] Committed."

if [[ "${DO_PUSH}" == "--push" ]]; then
  git push
  echo "[OK] Pushed to origin."
else
  echo "[INFO] Not pushed. To push: git push"
fi
