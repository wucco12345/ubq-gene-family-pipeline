#!/usr/bin/env bash
set -euo pipefail

F="UBQ"
D="results/${F}"
PROTEOME="inputs/proteome.fa"

# from run_UBQ.sh
CAND_RAW="${D}/${F}_candidates.id"          # PF00240 hits
FINAL_ID="${D}/${F}_final_all.id"           # poly_merged + ribofusion(strict)
POLY_ID="${D}/${F}_poly_merged.id"          # repeat>=2
RIBO_ID="${D}/${F}_keep_ribofusion.id"      # strict keep from BLAST
NOHIT_ID="${D}/${F}_nohit.id"               # repeat1 nohit
TOP1="${D}/evidence/${F}_blast_top1.tsv"    # repeat1 top1 evidence (if created by run_UBQ)

for f in "$PROTEOME" "$CAND_RAW" "$FINAL_ID" "$POLY_ID"; do
  test -s "$f" || { echo "[ERROR] missing/empty: $f"; exit 1; }
done
# optional files:
test -s "$RIBO_ID" || echo "[WARN] missing/empty (ok): $RIBO_ID"
test -s "$NOHIT_ID" || echo "[WARN] missing/empty (ok): $NOHIT_ID"
test -s "$TOP1" || echo "[WARN] missing/empty (ok but less traceability): $TOP1"

# ---------- 0) universe ----------
LC_ALL=C sort -u "$CAND_RAW" > "${D}/${F}_candidates.sorted.id"
LC_ALL=C sort -u "$FINAL_ID" > "${D}/${F}_final_all.sorted.id"

# ---------- 1) excluded_final = candidates - final_all ----------
comm -23 "${D}/${F}_candidates.sorted.id" "${D}/${F}_final_all.sorted.id" > "${D}/${F}_excluded_final.id"

# ---------- 2) subset ids (sorted) ----------
LC_ALL=C sort -u "$POLY_ID" > "${D}/${F}_poly_merged.sorted.id"
if test -s "$RIBO_ID"; then LC_ALL=C sort -u "$RIBO_ID" > "${D}/${F}_keep_ribofusion.sorted.id"; else : > "${D}/${F}_keep_ribofusion.sorted.id"; fi
if test -s "$NOHIT_ID"; then LC_ALL=C sort -u "$NOHIT_ID" > "${D}/${F}_nohit.sorted.id"; else : > "${D}/${F}_nohit.sorted.id"; fi

# ---------- 3) extract FASTA ----------
cat > "${D}/_extract_by_id.awk" <<'AWK'
BEGIN{
  while((getline line < IDS) > 0){
    gsub(/\r/, "", line)
    sub(/[ \t]+$/, "", line)
    if(line!="") want[line]=1
  }
  close(IDS)
  keep=0
}
/^>/{
  hdr=$0
  gsub(/\r/, "", hdr)
  sub(/^>/, "", hdr)
  split(hdr,a,/[ \t]/)
  id=a[1]
  keep = (id in want)
}
keep{
  gsub(/\r/, "", $0)
  print
}
AWK

awk -v IDS="${D}/${F}_final_all.sorted.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_final_all.fa"
awk -v IDS="${D}/${F}_poly_merged.sorted.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_poly_merged.fa"
awk -v IDS="${D}/${F}_keep_ribofusion.sorted.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_keep_ribofusion.fa"
awk -v IDS="${D}/${F}_excluded_final.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_excluded_final.fa"

# ---------- 4) evidence (top1) ----------
if test -s "$TOP1"; then
  grep -Ff "${D}/${F}_keep_ribofusion.sorted.id" "$TOP1" > "${D}/${F}_keep_ribofusion.top1.tsv" || true
  grep -Ff "${D}/${F}_excluded_final.id" "$TOP1" > "${D}/${F}_excluded_final.top1.tsv" || true
fi

# ---------- 5) QC ----------
{
  echo "UBQ QC (standardized v2)"
  echo "Candidates (PF00240): $(wc -l < "${D}/${F}_candidates.sorted.id")"
  echo "poly_merged (repeat>=2): $(wc -l < "${D}/${F}_poly_merged.sorted.id")"
  echo "keep_ribofusion (strict): $(wc -l < "${D}/${F}_keep_ribofusion.sorted.id")"
  echo "nohit (repeat1): $(wc -l < "${D}/${F}_nohit.sorted.id")"
  echo "Final ALL: $(wc -l < "${D}/${F}_final_all.sorted.id")"
  echo "Excluded final: $(wc -l < "${D}/${F}_excluded_final.id")"
  echo ""

  echo "[QC] final_all ∩ excluded_final (must be 0):"
  comm -12 "${D}/${F}_final_all.sorted.id" <(LC_ALL=C sort -u "${D}/${F}_excluded_final.id") | wc -l
  echo ""

  # sanity: poly should be subset of final_all
  echo "[QC] poly_merged - final_all (should be 0):"
  comm -23 "${D}/${F}_poly_merged.sorted.id" "${D}/${F}_final_all.sorted.id" | wc -l
  echo ""

  if test -s "$TOP1"; then
    echo "Top1 titles (keep_ribofusion; top 10):"
    cut -f7 "${D}/${F}_keep_ribofusion.top1.tsv" 2>/dev/null | sed 's/ OS=.*//' | LC_ALL=C sort | uniq -c | sort -nr | head -n 10 || true
    echo ""
    echo "Top1 titles (excluded_final; top 10):"
    cut -f7 "${D}/${F}_excluded_final.top1.tsv" 2>/dev/null | sed 's/ OS=.*//' | LC_ALL=C sort | uniq -c | sort -nr | head -n 10 || true
  else
    echo "[WARN] no top1 evidence file; only set arithmetic QC available."
  fi
} > "${D}/${F}_qc.txt"

echo "[DONE] UBQ standardized v2 outputs:"
ls -lh "${D}/${F}_final_all.id" "${D}/${F}_final_all.fa" \
      "${D}/${F}_excluded_final.id" "${D}/${F}_excluded_final.fa" \
      "${D}/${F}_poly_merged.fa" "${D}/${F}_keep_ribofusion.fa" \
      "${D}/${F}_qc.txt"
