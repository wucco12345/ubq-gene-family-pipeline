#!/usr/bin/env bash
set -euo pipefail

F="EF1A"
D="results/${F}"
PROTEOME="inputs/proteome.fa"

HIGH_ID="${D}/${F}_high_conf.id"
RESC_ID="${D}/${F}_rescue.id"
NOHIT_ID="${D}/${F}_nohit.id"
TOP1="${D}/evidence/${F}_blast_top1.tsv"

for f in "$HIGH_ID" "$RESC_ID" "$NOHIT_ID" "$PROTEOME"; do
  test -f "$f" || { echo "[ERROR] missing file: $f"; exit 1; }
done

if test -f "${D}/${F}_candidates.id"; then
  LC_ALL=C sort -u "${D}/${F}_candidates.id" > "${D}/${F}_candidates.sorted.id"
elif test -f "${D}/${F}_candidates.sorted.id"; then
  LC_ALL=C sort -u "${D}/${F}_candidates.sorted.id" > "${D}/${F}_candidates.sorted.id.tmp"
  mv "${D}/${F}_candidates.sorted.id.tmp" "${D}/${F}_candidates.sorted.id"
else
  echo "[ERROR] missing candidate id file"
  exit 1
fi

LC_ALL=C sort -u "$HIGH_ID" > "${D}/${F}_high_conf.sorted.id"
LC_ALL=C sort -u "$RESC_ID" > "${D}/${F}_rescue.sorted.id"
LC_ALL=C sort -u "$NOHIT_ID" > "${D}/${F}_nohit.sorted.id"

cat "${D}/${F}_high_conf.sorted.id" "${D}/${F}_rescue.sorted.id" | LC_ALL=C sort -u > "${D}/${F}_final_union.id"

comm -23 "${D}/${F}_candidates.sorted.id" "${D}/${F}_final_union.id" > "${D}/${F}_excluded.id"

awk -v IDS="${D}/${F}_high_conf.sorted.id" -f scripts/lib/extract_by_id.awk "$PROTEOME" > "${D}/${F}_high_conf.fa"
awk -v IDS="${D}/${F}_rescue.sorted.id" -f scripts/lib/extract_by_id.awk "$PROTEOME" > "${D}/${F}_rescue.fa"
awk -v IDS="${D}/${F}_excluded.id" -f scripts/lib/extract_by_id.awk "$PROTEOME" > "${D}/${F}_excluded.fa"

if test -s "$TOP1"; then
  grep -Ff "${D}/${F}_high_conf.sorted.id" "$TOP1" > "${D}/${F}_high_conf.top1.tsv" || true
  grep -Ff "${D}/${F}_rescue.sorted.id" "$TOP1" > "${D}/${F}_rescue.top1.tsv" || true
  grep -Ff "${D}/${F}_excluded.id" "$TOP1" > "${D}/${F}_excluded.top1.tsv" || true
fi

{
  echo -e "ID\treason\ttop1_stitle"
  while read -r id; do
    [ -z "$id" ] && continue
    if grep -Fxq "$id" "${D}/${F}_nohit.sorted.id"; then
      echo -e "${id}\tblast_nohit\tNA"
    else
      title=$(awk -F'\t' -v q="$id" '$1==q{print $7; exit}' "$TOP1" 2>/dev/null || true)
      [ -z "${title:-}" ] && title="NA"
      echo -e "${id}\tnot_in_final_set\t${title}"
    fi
  done < "${D}/${F}_excluded.id"
} > "${D}/${F}_excluded.tsv"

{
  echo "EF1A QC (standardized v2)"
  echo "Candidates: $(wc -l < "${D}/${F}_candidates.sorted.id")"
  echo "High_conf: $(wc -l < "${D}/${F}_high_conf.sorted.id")"
  echo "Rescue: $(wc -l < "${D}/${F}_rescue.sorted.id")"
  echo "Final total: $(wc -l < "${D}/${F}_final_union.id")"
  echo "NoHit: $(wc -l < "${D}/${F}_nohit.sorted.id")"
  echo "Excluded: $(wc -l < "${D}/${F}_excluded.id")"
  echo
  echo "[QC] high_conf ∩ rescue (must be 0):"
  comm -12 "${D}/${F}_high_conf.sorted.id" "${D}/${F}_rescue.sorted.id" | wc -l
  echo
  echo "[QC] final_union ∩ excluded (must be 0):"
  comm -12 "${D}/${F}_final_union.id" "${D}/${F}_excluded.id" | wc -l
} > "${D}/${F}_qc.txt"

echo "[DONE] standardized outputs in ${D}"
