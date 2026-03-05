#!/usr/bin/env bash
set -euo pipefail

F="ACT"
D="results/${F}"
PROTEOME="inputs/proteome.fa"

TOP1="${D}/evidence/${F}_blast_top1.tsv"
TOP5="${D}/evidence/${F}_vs_sprot.tsv"
CAND_RAW="${D}/${F}.id"   # run_ACT.sh 生成的 PF00022 候选 ID

for f in "$PROTEOME" "$TOP1" "$TOP5" "$CAND_RAW"; do
  test -s "$f" || { echo "[ERROR] missing/empty: $f"; exit 1; }
done

# ---------- 0) normalize candidates ----------
LC_ALL=C sort -u "$CAND_RAW" > "${D}/${F}_candidates.sorted.id"

# ---------- 1) top5 consensus per qid ----------
awk -F'\t' 'BEGIN{IGNORECASE=1; OFS="\t"}
{
  q=$1; t=$7
  # ACT-like: has actin but NOT actin-related/arp
  if (t ~ /actin/ && t !~ /actin-related|arp/) act[q]++
  # ARP-like: has actin-related OR arp
  if (t ~ /actin-related|arp/) arp[q]++
  n[q]++
}
END{
  for (q in n){
    if (!(q in act)) act[q]=0
    if (!(q in arp)) arp[q]=0
    print q, act[q], arp[q], n[q]
  }
}' "$TOP5" | LC_ALL=C sort -k1,1 > "${D}/${F}_top5_consensus.tsv"
# qid act_like arp_like topN

# ---------- 2) join with top1 title ----------
cut -f1,7 "$TOP1" | LC_ALL=C sort -k1,1 > "${D}/${F}_top1_title.tsv"

join -t $'\t' -a 1 -e 'NA' -o '0,1.2,2.2,2.3,2.4' \
  "${D}/${F}_top1_title.tsv" "${D}/${F}_top5_consensus.tsv" \
  > "${D}/${F}_classify.tsv"
# qid top1_title act_like arp_like topN

# ---------- 3) classify (NP-safe: top1 + top5 consensus) ----------
: > "${D}/${F}_final_ACT.id"
: > "${D}/${F}_final_ARP.id"
: > "${D}/${F}_excluded_final.id"

awk -F'\t' 'BEGIN{IGNORECASE=1}
{
  q=$1; title=$2; act=$3+0; arp=$4+0;

  is_act_top1 = (title ~ /actin/ && title !~ /actin-related|arp/);
  is_arp_top1 = (title ~ /actin-related|arp/);

  # strong consensus thresholds:
  #   strict: top1 consistent AND >=3/5 (or >=3/topN if topN<5)
  #   rescue: top1 odd but consensus very strong (>=4)
  if (is_act_top1 && act>=3) {print q >> "'"${D}/${F}_final_ACT.id"'"; next}
  if (is_arp_top1 && arp>=3) {print q >> "'"${D}/${F}_final_ARP.id"'"; next}

  if (!is_act_top1 && act>=4) {print q >> "'"${D}/${F}_final_ACT.id"'"; next}
  if (!is_arp_top1 && arp>=4) {print q >> "'"${D}/${F}_final_ARP.id"'"; next}

  print q >> "'"${D}/${F}_excluded_final.id"'"
}' "${D}/${F}_classify.tsv"

cat "${D}/${F}_final_ACT.id" "${D}/${F}_final_ARP.id" | LC_ALL=C sort -u > "${D}/${F}_final_all.id"
LC_ALL=C sort -u "${D}/${F}_final_all.id" > "${D}/${F}_final_all.sorted.id"
LC_ALL=C sort -u "${D}/${F}_excluded_final.id" > "${D}/${F}_excluded_final.sorted.id"

# QC: disjoint
echo -n "[QC] final_all ∩ excluded_final = "
comm -12 "${D}/${F}_final_all.sorted.id" "${D}/${F}_excluded_final.sorted.id" | wc -l

# ---------- 4) extract FASTA helper ----------
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

# final / excluded FASTA
awk -v IDS="${D}/${F}_final_all.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_final_all.fa"
awk -v IDS="${D}/${F}_excluded_final.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_excluded_final.fa"
awk -v IDS="${D}/${F}_final_ACT.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_final_ACT.fa"
awk -v IDS="${D}/${F}_final_ARP.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_final_ARP.fa"

# ---------- 5) evidence snapshots ----------
grep -Ff "${D}/${F}_final_ACT.id" "$TOP1" > "${D}/${F}_final_ACT.top1.tsv" || true
grep -Ff "${D}/${F}_final_ARP.id" "$TOP1" > "${D}/${F}_final_ARP.top1.tsv" || true
grep -Ff "${D}/${F}_excluded_final.id" "$TOP1" > "${D}/${F}_excluded_final.top1.tsv" || true

# ---------- 6) QC report ----------
{
  echo "ACT QC (standardized v2)"
  echo "Candidates (PF00022): $(wc -l < "${D}/${F}_candidates.sorted.id")"
  echo "Final ACT: $(wc -l < "${D}/${F}_final_ACT.id")"
  echo "Final ARP: $(wc -l < "${D}/${F}_final_ARP.id")"
  echo "Final ALL: $(wc -l < "${D}/${F}_final_all.id")"
  echo "Excluded final: $(wc -l < "${D}/${F}_excluded_final.id")"
  echo ""
  echo "[QC] excluded_final contains obvious ACT top1 (should be 0):"
  awk -F'\t' 'BEGIN{IGNORECASE=1} $7 ~ /actin/ && $7 !~ /actin-related|arp/ {print}' "${D}/${F}_excluded_final.top1.tsv" | wc -l
  echo ""
  echo "Top1 ACT examples:"
  head -n 5 "${D}/${F}_final_ACT.top1.tsv" 2>/dev/null || true
  echo ""
  echo "Top1 ARP examples:"
  head -n 5 "${D}/${F}_final_ARP.top1.tsv" 2>/dev/null || true
} > "${D}/${F}_qc.txt"

echo "[DONE] ACT standardized v2 outputs:"
ls -lh \
  "${D}/${F}_final_all.id" "${D}/${F}_final_all.fa" \
  "${D}/${F}_excluded_final.id" "${D}/${F}_excluded_final.fa" \
  "${D}/${F}_qc.txt" \
  "${D}/${F}_final_ACT.id" "${D}/${F}_final_ACT.fa" \
  "${D}/${F}_final_ARP.id" "${D}/${F}_final_ARP.fa"
