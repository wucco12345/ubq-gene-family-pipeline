#!/usr/bin/env bash
set -euo pipefail

F="EF1A"
D="results/${F}"
PROTEOME="inputs/proteome.fa"

TOP1="${D}/evidence/${F}_blast_top1.tsv"
CAND_RAW="${D}/${F}_candidates.id"
HIGH_ID="${D}/${F}_high_conf.id"
RESC_ID="${D}/${F}_rescue.id"

for f in "$PROTEOME" "$CAND_RAW" "$HIGH_ID" "$RESC_ID"; do
  test -s "$f" || { echo "[ERROR] missing/empty: $f"; exit 1; }
done
test -s "$TOP1" || echo "[WARN] missing/empty (ok but less traceability): $TOP1"

# ---------- 0) universe ----------
LC_ALL=C sort -u "$CAND_RAW" > "${D}/${F}_candidates.sorted.id"

# ---------- 1) final_all = high_conf + rescue ----------
cat "$HIGH_ID" "$RESC_ID" | LC_ALL=C sort -u > "${D}/${F}_final_all.id"

# ---------- 2) excluded_final = candidates - final_all ----------
comm -23 "${D}/${F}_candidates.sorted.id" <(LC_ALL=C sort -u "${D}/${F}_final_all.id") \
  > "${D}/${F}_excluded_final.id"

# ---------- 3) evidence (top1) ----------
if test -s "$TOP1"; then
  grep -Ff "${D}/${F}_final_all.id" "$TOP1" > "${D}/${F}_final_all.top1.tsv" || true
  grep -Ff "${D}/${F}_excluded_final.id" "$TOP1" > "${D}/${F}_excluded_final.top1.tsv" || true
fi

# ---------- 4) extract FASTA ----------
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

awk -v IDS="${D}/${F}_final_all.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_final_all.fa"
awk -v IDS="${D}/${F}_excluded_final.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_excluded_final.fa"

# ---------- 5) QC ----------
{
  echo "EF1A QC (standardized v2)"
  echo "Candidates (PF00009): $(wc -l < "${D}/${F}_candidates.sorted.id")"
  echo "High_conf: $(wc -l < "$HIGH_ID")"
  echo "Rescue: $(wc -l < "$RESC_ID")"
  echo "Final ALL: $(wc -l < "${D}/${F}_final_all.id")"
  echo "Excluded final: $(wc -l < "${D}/${F}_excluded_final.id")"
  echo ""
  echo "[QC] final_all ∩ excluded_final (must be 0):"
  comm -12 <(LC_ALL=C sort -u "${D}/${F}_final_all.id") <(LC_ALL=C sort -u "${D}/${F}_excluded_final.id") | wc -l
  echo ""
  if test -s "$TOP1"; then
    echo "Top1 titles (final_all; top 10):"
    cut -f7 "${D}/${F}_final_all.top1.tsv" 2>/dev/null | sed 's/ OS=.*//' | LC_ALL=C sort | uniq -c | sort -nr | head -n 10 || true
    echo ""
    echo "Top1 titles (excluded_final; top 10):"
    cut -f7 "${D}/${F}_excluded_final.top1.tsv" 2>/dev/null | sed 's/ OS=.*//' | LC_ALL=C sort | uniq -c | sort -nr | head -n 10 || true
  else
    echo "[WARN] no top1 evidence file; only set arithmetic QC available."
  fi
} > "${D}/${F}_qc.txt"

echo "[DONE] EF1A standardized v2 outputs:"
ls -lh "${D}/${F}_final_all.id" "${D}/${F}_final_all.fa" "${D}/${F}_excluded_final.id" "${D}/${F}_excluded_final.fa" "${D}/${F}_qc.txt"
