#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib.sh"

D="results/UBQ"
test -s "$D/UBQ_final_all.id" || die "missing $D/UBQ_final_all.id"
test -s "$D/UBQ_final_all.fa" || die "missing $D/UBQ_final_all.fa"
test -s "$D/UBQ_qc.txt" || die "missing $D/UBQ_qc.txt"

EX="example/UBQ"
mkdir -p "$EX"

cp -f "$D/UBQ_final_all.id" "$EX/UBQ_final_all.id"
cp -f "$D/UBQ_final_all.fa" "$EX/UBQ_final_all.fa"
cp -f "$D/UBQ_high_conf.id" "$EX/UBQ_high_conf.id"
cp -f "$D/UBQ_qc.txt" "$EX/UBQ_qc.txt"

# evidence: top1 table is useful and small-ish; include it if exists
if test -s "$D/evidence/UBQ_blast_top1.tsv"; then
  cp -f "$D/evidence/UBQ_blast_top1.tsv" "$EX/UBQ_blast_top1.tsv"
fi

echo "[DONE] wrote $EX"
