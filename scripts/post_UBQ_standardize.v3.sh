#!/usr/bin/env bash
set -euo pipefail

F="UBQ"
D="results/${F}"
PROTEOME="inputs/proteome.fa"

# Inputs from run_UBQ.sh
CAND_RAW="${D}/${F}_candidates.id"                 # PF00240 candidates
POLY_ID="${D}/${F}_poly_merged.id"                 # repeat>=2 (true polyubq)
REP1_ID="${D}/${F}_repeat1.id"                     # repeat==1
TOP1="${D}/evidence/${F}_blast_top1.tsv"           # repeat1 top1 (from blast)
TOP5="${D}/evidence/${F}_vs_sprot.tsv"             # repeat1 top5 (from blast)

for f in "$PROTEOME" "$CAND_RAW" "$POLY_ID" "$REP1_ID" "$TOP1" "$TOP5"; do
  test -s "$f" || { echo "[ERROR] missing/empty: $f"; exit 1; }
done

mkdir -p "${D}/evidence/repeat1_desc"

# ---------- 0) normalize candidate universe ----------
LC_ALL=C sort -u "$CAND_RAW" > "${D}/${F}_candidates.sorted.id"
LC_ALL=C sort -u "$POLY_ID"  > "${D}/${F}_poly_merged.sorted.id"
LC_ALL=C sort -u "$REP1_ID"  > "${D}/${F}_repeat1.sorted.id"

# ---------- helper: extract fasta by id ----------
EXTR_AWK="${D}/_extract_by_id.awk"
cat > "$EXTR_AWK" <<'AWK'
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

extract_fa () {
  local idf="$1" outfa="$2"
  : > "$outfa"
  if test -s "$idf"; then
    awk -v IDS="$idf" -f "$EXTR_AWK" "$PROTEOME" > "$outfa"
  fi
}

# ---------- 1) repeat1: collect top1/top5 for traceability ----------
grep -Ff "${D}/${F}_repeat1.sorted.id" "$TOP1" > "${D}/evidence/repeat1_desc/${F}_repeat1_top1.tsv" || true
grep -Ff "${D}/${F}_repeat1.sorted.id" "$TOP5" > "${D}/evidence/repeat1_desc/${F}_repeat1_top5.tsv" || true

# Top1 title counts (human-readable)
cut -f7 "${D}/evidence/repeat1_desc/${F}_repeat1_top1.tsv" \
  | sed 's/ OS=.*//' | sort | uniq -c | sort -nr \
  > "${D}/evidence/repeat1_desc/${F}_repeat1_top1.stitle.counts.txt" || true

# ---------- 2) repeat1 classification by top5 consensus (NP-safe) ----------
# We classify repeat1 into:
#   keep_polyubq: strong positive signal for polyubiquitin/ubiquitin (>=3 pos in top5) AND no strong negatives
#   remove_nonubq: strong negative (RUB/NEDD8/SUMO or clearly non-UBQ protein)
#   ambiguous: everything else (for manual review)

# Keywords:
# POS (true UBQ-like): polyubiquitin / ubiquitin (but beware "ubiquitin-like")
# NEG (exclude): RUB1/RUB2/NEDD8/SUMO, "ubiquitin-like domain-containing", obvious non-ubq enzymes (PI4K)

awk -F'\t' 'BEGIN{IGNORECASE=1; OFS="\t"}
{
  q=$1; t=$7

  # positive: polyubiquitin or canonical ubiquitin title
  if (t ~ /polyubiquitin/) pos[q]++
  # "ubiquitin" as a word (avoid counting ubiquitin-like)
  if (t ~ /(^|[^-])ubiquitin([^a-z]|$)/ && t !~ /ubiquitin-like/) pos[q]++

  # negative: RUB/NEDD8/SUMO
  if (t ~ /(rub1|rub2|nedd8|sumo)/) neg[q]++
  # negative: UBL-domain containing (big proteins)
  if (t ~ /(ubiquitin-like domain-containing)/) neg[q]++
  # negative: obvious non-ubq proteins seen in your data
  if (t ~ /(phosphatidylinositol 4-kinase|pi4k|p4kg)/) neg[q]++

  n[q]++
}
END{
  for(q in n){
    if(!(q in pos)) pos[q]=0
    if(!(q in neg)) neg[q]=0
    print q, pos[q], neg[q], n[q]
  }
}' "${D}/evidence/repeat1_desc/${F}_repeat1_top5.tsv" \
  | sort -k1,1 > "${D}/evidence/repeat1_desc/${F}_repeat1_top5_keywords.tsv"
# columns: qid pos_hits neg_hits topN

# Merge with top1 title for classification traceability
cut -f1,7 "${D}/evidence/repeat1_desc/${F}_repeat1_top1.tsv" | sort -k1,1 > "${D}/evidence/repeat1_desc/${F}_repeat1_top1_title.tsv"

join -t $'\t' -a 1 -e 'NA' -o '0,1.2,2.2,2.3,2.4' \
  "${D}/evidence/repeat1_desc/${F}_repeat1_top1_title.tsv" \
  "${D}/evidence/repeat1_desc/${F}_repeat1_top5_keywords.tsv" \
  > "${D}/evidence/repeat1_desc/${F}_repeat1_classify.tsv"
# columns: qid top1_title pos neg topN

: > "${D}/${F}_repeat1_keep_polyubq.id"
: > "${D}/${F}_repeat1_remove_nonubq.id"
: > "${D}/${F}_repeat1_ambiguous.id"

awk -F'\t' 'BEGIN{OFS="\t"; IGNORECASE=1}
{
  q=$1; title=$2; pos=$3+0; neg=$4+0; topN=$5+0;

  # Hard negative if top1 clearly negative OR neg hits strong
  hard_neg = (title ~ /(rub1|rub2|nedd8|sumo|phosphatidylinositol 4-kinase|pi4k|p4kg|ubiquitin-like domain-containing)/) || (neg>=2)

  # Hard positive if top1 is polyubiquitin/ubiquitin AND pos strong
  top1_pos = (title ~ /polyubiquitin/) || ((title ~ /(^|[^-])ubiquitin([^a-z]|$)/) && title !~ /ubiquitin-like/)
  hard_pos = (top1_pos && pos>=3 && neg==0)

  if (hard_pos) { print q > "'"${D}/${F}_repeat1_keep_polyubq.id"'"; next }
  if (hard_neg) { print q > "'"${D}/${F}_repeat1_remove_nonubq.id"'"; next }

  # Rescue: very strong pos even if top1 label not perfect
  if (pos>=4 && neg==0) { print q > "'"${D}/${F}_repeat1_keep_polyubq.id"'"; next }

  print q > "'"${D}/${F}_repeat1_ambiguous.id"'"
}' "${D}/evidence/repeat1_desc/${F}_repeat1_classify.tsv"

LC_ALL=C sort -u "${D}/${F}_repeat1_keep_polyubq.id"    -o "${D}/${F}_repeat1_keep_polyubq.id"
LC_ALL=C sort -u "${D}/${F}_repeat1_remove_nonubq.id"   -o "${D}/${F}_repeat1_remove_nonubq.id"
LC_ALL=C sort -u "${D}/${F}_repeat1_ambiguous.id"       -o "${D}/${F}_repeat1_ambiguous.id"

# Evidence extraction (top1 lines)
grep -Ff "${D}/${F}_repeat1_keep_polyubq.id"  "${D}/evidence/repeat1_desc/${F}_repeat1_top1.tsv" > "${D}/${F}_repeat1_keep_polyubq.top1.tsv" || true
grep -Ff "${D}/${F}_repeat1_remove_nonubq.id" "${D}/evidence/repeat1_desc/${F}_repeat1_top1.tsv" > "${D}/${F}_repeat1_remove_nonubq.top1.tsv" || true
grep -Ff "${D}/${F}_repeat1_ambiguous.id"     "${D}/evidence/repeat1_desc/${F}_repeat1_top1.tsv" > "${D}/${F}_repeat1_ambiguous.top1.tsv" || true

# ---------- 3) final_all = poly_merged + repeat1_keep_polyubq ----------
cat "${D}/${F}_poly_merged.sorted.id" "${D}/${F}_repeat1_keep_polyubq.id" 2>/dev/null \
  | LC_ALL=C sort -u > "${D}/${F}_final_all.id"

LC_ALL=C sort -u "${D}/${F}_final_all.id" > "${D}/${F}_final_all.sorted.id"

# excluded_final = candidates - final_all
comm -23 "${D}/${F}_candidates.sorted.id" "${D}/${F}_final_all.sorted.id" > "${D}/${F}_excluded_final.id"

# ---------- 4) extract FASTA for final/excluded ----------
extract_fa "${D}/${F}_final_all.id"     "${D}/${F}_final_all.fa"
extract_fa "${D}/${F}_excluded_final.id" "${D}/${F}_excluded_final.fa"
extract_fa "${D}/${F}_poly_merged.sorted.id" "${D}/${F}_poly_merged.fa"
extract_fa "${D}/${F}_repeat1_keep_polyubq.id" "${D}/${F}_repeat1_keep_polyubq.fa"

# ---------- 5) QC ----------
{
  echo "UBQ QC (standardized v3)"
  echo "Candidates (PF00240): $(wc -l < "${D}/${F}_candidates.sorted.id")"
  echo "Polyubiquitin repeat>=2: $(wc -l < "${D}/${F}_poly_merged.sorted.id")"
  echo "Repeat1 total: $(wc -l < "${D}/${F}_repeat1.sorted.id")"
  echo "Repeat1 keep_polyubq: $(wc -l < "${D}/${F}_repeat1_keep_polyubq.id")"
  echo "Repeat1 remove_nonubq: $(wc -l < "${D}/${F}_repeat1_remove_nonubq.id")"
  echo "Repeat1 ambiguous: $(wc -l < "${D}/${F}_repeat1_ambiguous.id")"
  echo "FINAL all: $(wc -l < "${D}/${F}_final_all.id")"
  echo "Excluded_final: $(wc -l < "${D}/${F}_excluded_final.id")"
  echo ""

  echo "[QC] final_all ∩ excluded_final (must be 0):"
  comm -12 <(LC_ALL=C sort -u "${D}/${F}_final_all.id") <(LC_ALL=C sort -u "${D}/${F}_excluded_final.id") | wc -l
  echo ""

  echo "[QC] repeat1_keep top1 title examples:"
  head -n 10 "${D}/${F}_repeat1_keep_polyubq.top1.tsv" 2>/dev/null || true
  echo ""

  echo "[QC] repeat1_remove top1 title examples:"
  head -n 10 "${D}/${F}_repeat1_remove_nonubq.top1.tsv" 2>/dev/null || true
  echo ""

  echo "[QC] repeat1 ambiguous top1 title examples:"
  head -n 10 "${D}/${F}_repeat1_ambiguous.top1.tsv" 2>/dev/null || true
} > "${D}/${F}_qc.txt"

echo "[DONE] UBQ standardized v3 outputs:"
ls -lh \
  "${D}/${F}_final_all.id" "${D}/${F}_final_all.fa" \
  "${D}/${F}_excluded_final.id" "${D}/${F}_excluded_final.fa" \
  "${D}/${F}_qc.txt" \
  "${D}/${F}_repeat1_keep_polyubq.id" "${D}/${F}_repeat1_remove_nonubq.id" "${D}/${F}_repeat1_ambiguous.id" \
  "${D}/${F}_poly_merged.fa" "${D}/${F}_repeat1_keep_polyubq.fa"
