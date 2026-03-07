#!/usr/bin/env bash
set -euo pipefail

ACT_DIR="results/ACT"
TOP1="${ACT_DIR}/evidence/ACT_blast_top1.tsv"
TOP5="${ACT_DIR}/evidence/ACT_vs_sprot.tsv"
CAND_RAW="${ACT_DIR}/ACT.id"

for f in "$TOP1" "$TOP5" "$CAND_RAW"; do
  test -s "$f" || { echo "[ERROR] missing/empty: $f"; exit 1; }
done

# ---- Helper: normalize ids ----
sort -u "$CAND_RAW" > "${ACT_DIR}/ACT_candidates.id"

# ---- 1) Build top5 consensus counts per query ----
# Count in top5:
#   ACT-like: contains 'actin' but NOT 'actin-related' or 'arp'
#   ARP-like: contains 'actin-related' OR 'arp'
# We count among top5 lines for each query id.
awk -F'\t' 'BEGIN{IGNORECASE=1}
{
  q=$1; title=$7;
  # actin but not actin-related/arp
  if (title ~ /actin/ && title !~ /actin-related|arp/) act[q]++
  # actin-related / arp
  if (title ~ /actin-related|arp/) arp[q]++
  total[q]++
}
END{
  for (q in total){
    # ensure keys exist
    if (!(q in act)) act[q]=0;
    if (!(q in arp)) arp[q]=0;
    print q"\t"act[q]"\t"arp[q]"\t"total[q]
  }
}' "$TOP5" | sort -k1,1 > "${ACT_DIR}/ACT_top5_consensus.tsv"
# columns: qid  act_like_count  arp_like_count  topN(should be <=5)

# ---- 2) Classify by top1 + top5 consensus (>=3/5) ----
# Prepare top1 title table
cut -f1,7 "$TOP1" | sort -k1,1 > "${ACT_DIR}/ACT_top1_title.tsv"

join -t $'\t' -a 1 -e 'NA' -o '0,1.2,2.2,2.3,2.4' \
  "${ACT_DIR}/ACT_top1_title.tsv" "${ACT_DIR}/ACT_top5_consensus.tsv" \
  > "${ACT_DIR}/ACT_classify.tsv"
# columns: qid  top1_title  act_like  arp_like  topN

awk -F'\t' 'BEGIN{IGNORECASE=1}
{
  q=$1; title=$2; act=$3+0; arp=$4+0;
  # rule:
  # ACT: top1 has actin (not actin-related/arp) AND act>=3
  # ARP: top1 has actin-related/arp AND arp>=3
  is_act_top1 = (title ~ /actin/ && title !~ /actin-related|arp/)
  is_arp_top1 = (title ~ /actin-related|arp/)
  if (is_act_top1 && act>=3) {print q > "'"${ACT_DIR}/ACT_final_ACT.id"'"; next}
  if (is_arp_top1 && arp>=3) {print q > "'"${ACT_DIR}/ACT_final_ARP.id"'"; next}
  # if top1 is odd, but consensus very strong, allow rescue:
  if (!is_act_top1 && act>=4) {print q > "'"${ACT_DIR}/ACT_final_ACT.id"'"; next}
  if (!is_arp_top1 && arp>=4) {print q > "'"${ACT_DIR}/ACT_final_ARP.id"'"; next}
  print q > "'"${ACT_DIR}/ACT_excluded_final.id"'"
}' "${ACT_DIR}/ACT_classify.tsv"

# ---- 3) Final all ----
cat "${ACT_DIR}/ACT_final_ACT.id" "${ACT_DIR}/ACT_final_ARP.id" 2>/dev/null | sort -u > "${ACT_DIR}/ACT_final_all.id"

# ---- 4) Evidence extraction for NP traceability ----
grep -Ff "${ACT_DIR}/ACT_final_ACT.id" "$TOP1" > "${ACT_DIR}/ACT_final_ACT.top1.tsv" || true
grep -Ff "${ACT_DIR}/ACT_final_ARP.id" "$TOP1" > "${ACT_DIR}/ACT_final_ARP.top1.tsv" || true
grep -Ff "${ACT_DIR}/ACT_excluded_final.id" "$TOP1" > "${ACT_DIR}/ACT_excluded_final.top1.tsv" || true

# ---- 5) QC ----
{
  echo "ACT QC (standardized)"
  echo "Candidates (PF00022): $(wc -l < "${ACT_DIR}/ACT_candidates.id")"
  echo "Final ACT: $(wc -l < "${ACT_DIR}/ACT_final_ACT.id")"
  echo "Final ARP: $(wc -l < "${ACT_DIR}/ACT_final_ARP.id")"
  echo "Final ALL: $(wc -l < "${ACT_DIR}/ACT_final_all.id")"
  echo "Excluded final: $(wc -l < "${ACT_DIR}/ACT_excluded_final.id")"
  echo ""

  echo "[QC] final_all ∩ excluded_final (must be 0):"
  comm -12 <(sort -u "${ACT_DIR}/ACT_final_all.id") <(sort -u "${ACT_DIR}/ACT_excluded_final.id") | wc -l
  echo ""

  echo "[QC] excluded_final contains obvious ACT top1 (should be 0):"
  awk -F'\t' 'BEGIN{IGNORECASE=1} $7 ~ /actin/ && $7 !~ /actin-related|arp/ {print}' "${ACT_DIR}/ACT_excluded_final.top1.tsv" | wc -l
  echo ""
  echo "Top1 ACT examples:"
  head -n 5 "${ACT_DIR}/ACT_final_ACT.top1.tsv" 2>/dev/null || true
  echo ""
  echo "Top1 ARP examples:"
  head -n 5 "${ACT_DIR}/ACT_final_ARP.top1.tsv" 2>/dev/null || true
} > "${ACT_DIR}/ACT_qc.txt"

echo "[DONE] Standardized classification written:"
ls -lh "${ACT_DIR}/ACT_final_ACT.id" "${ACT_DIR}/ACT_final_ARP.id" "${ACT_DIR}/ACT_final_all.id" "${ACT_DIR}/ACT_qc.txt"
