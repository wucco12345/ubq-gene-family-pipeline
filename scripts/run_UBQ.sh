#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
source scripts/lib/yaml.sh

F="UBQ"
CFG="config/UBQ.yml"
D="results/${F}"
EVD="${D}/evidence"
PROTEOME="inputs/proteome.fa"
PFAM="results/pfam.tsv"
DB="db/uniprot/sprot"

mkdir -p "$EVD"

# ---- read config ----
PF_RE="$(yget "$CFG" pfam.pf_re)"
IE_COL="$(yget "$CFG" pfam.ie_col 9)"
PROT_COL="$(yget "$CFG" pfam.protein_col 3)"
ALI_FROM_COL="$(yget "$CFG" pfam.ali_from_col 12)"
ALI_TO_COL="$(yget "$CFG" pfam.ali_to_col 13)"

BLAST_TASK="$(yget "$CFG" blast.task blastp)"
BLAST_TASK_SHORT="$(yget "$CFG" blast.task_short blastp-short)"
BLAST_EVALUE="$(yget "$CFG" blast.evalue 1e-10)"

POS_RE="$(yget "$CFG" rules.pos_re)"
NEG_RE="$(yget "$CFG" rules.neg_re)"

log "[1/7] Pfam candidates: PF00240 (iE<=1e-5)"
pfam_pick_ids "$PFAM" "$PF_RE" "$IE_COL" "$PROT_COL" "${D}/${F}_candidates.id"
cp "${D}/${F}_candidates.id" "${D}/${F}_candidates.raw.id"
log "Candidates: $(wc -l < ${D}/${F}_candidates.id)"

log "[2/7] Extract candidate FASTA"
extract_fa_by_id "${D}/${F}_candidates.id" "$PROTEOME" "${D}/${F}_candidates.fa"
log "Candidate FASTA seqs: $(grep -c '^>' ${D}/${F}_candidates.fa)"

log "[3/7] Save PFam hit intervals + merged repeat counts"
# evidence: qid ali_from ali_to ie pfacc
awk -F'\t' -v re="$PF_RE" -v IEC="$IE_COL" -v PC="$PROT_COL" -v AF="$ALI_FROM_COL" -v AT="$ALI_TO_COL" '
  BEGIN{OFS="\t"}
  {gsub(/\r/,"",$0)}
  ($2 ~ re) {
    ie=$(IEC)+0
    if(ie<=1e-5){
      id=$(PC); gsub(/\r/,"",id); sub(/[ \t]+$/,"",id)
      print id, $(AF)+0, $(AT)+0, ie, $2
    }
  }' "$PFAM" | LC_ALL=C sort -k1,1 -k2,2n > "${EVD}/${F}_pfam_hits.tsv"

# merge overlapping/adjacent intervals per protein; count merged segments
awk '
  BEGIN{OFS="\t"}
  {
    id=$1; s=$2+0; e=$3+0
    if(id!=cur && cur!=""){
      print cur, seg
    }
    if(id!=cur){
      cur=id; seg=1; ms=s; me=e
    } else {
      if(s <= me+1){
        if(e>me) me=e
      } else {
        seg++
        ms=s; me=e
      }
    }
  }
  END{ if(cur!="") print cur, seg }
' "${EVD}/${F}_pfam_hits.tsv" > "${D}/${F}_repeat_merged_counts.tsv"

awk '$2>=2{print $1}' "${D}/${F}_repeat_merged_counts.tsv" | LC_ALL=C sort -u > "${D}/${F}_poly_merged.id"
awk '$2==1{print $1}' "${D}/${F}_repeat_merged_counts.tsv" | LC_ALL=C sort -u > "${D}/${F}_repeat1.id"
log "poly_merged (repeat>=2): $(wc -l < ${D}/${F}_poly_merged.id)"
log "repeat1 (repeat==1): $(wc -l < ${D}/${F}_repeat1.id)"

log "[4/7] Extract repeat1 FASTA"
extract_fa_by_id "${D}/${F}_repeat1.id" "$PROTEOME" "${D}/${F}_repeat1.fa"
log "repeat1 FASTA seqs: $(grep -c '^>' ${D}/${F}_repeat1.fa)"

log "[5/7] BLASTp repeat1 vs Swiss-Prot (top5 -> top1)"
# repeat1 可能很短：默认用 blastp-short 更稳；但你也可以在 config 里把 task_short 改成 blastp
blast_top5 "${D}/${F}_repeat1.fa" "$DB" "${EVD}/${F}_vs_sprot.tsv" "$THREADS_DEFAULT" "$BLAST_TASK_SHORT" "$BLAST_EVALUE"
blast_top1_from_top5 "${EVD}/${F}_vs_sprot.tsv" "${EVD}/${F}_blast_top1.tsv" "${D}/${F}_top1.id"

log "[6/7] Classify repeat1 by top1 (configurable pos/neg); nohit handled"
# nohit = repeat1 - top1.id
comm -23 <(LC_ALL=C sort -u "${D}/${F}_repeat1.id") <(LC_ALL=C sort -u "${D}/${F}_top1.id") > "${D}/${F}_nohit.id"

# keep/remove (based on top1 stitle)
awk -F'\t' -v POS_RE="$POS_RE" -v NEG_RE="$NEG_RE" 'BEGIN{IGNORECASE=1}
{
  q=$1; t=$7
  if(t ~ NEG_RE){print q > "'"${D}/${F}_remove_by_blast.id"'"; next}
  if(t ~ POS_RE){print q > "'"${D}/${F}_keep_ribofusion.id"'"; next}  # 文件名沿用旧名（即使现在不是 ribofusion）
}' "${EVD}/${F}_blast_top1.tsv"

# ensure files exist even if empty
test -f "${D}/${F}_keep_ribofusion.id" || : > "${D}/${F}_keep_ribofusion.id"
test -f "${D}/${F}_remove_by_blast.id" || : > "${D}/${F}_remove_by_blast.id"

# strict UBQ final: poly_merged + keep_from_repeat1
cat "${D}/${F}_poly_merged.id" "${D}/${F}_keep_ribofusion.id" 2>/dev/null | LC_ALL=C sort -u > "${D}/${F}_final_all.id"
cp "${D}/${F}_final_all.id" "${D}/${F}_high_conf.id"
: > "${D}/${F}_rescue.id"

log "[7/7] Extract final FASTA + excluded_final set"
extract_fa_by_id "${D}/${F}_final_all.id" "$PROTEOME" "${D}/${F}_final_all.fa"
extract_fa_by_id "${D}/${F}_poly_merged.id" "$PROTEOME" "${D}/${F}_poly_merged.fa" || true
extract_fa_by_id "${D}/${F}_keep_ribofusion.id" "$PROTEOME" "${D}/${F}_keep_ribofusion.fa" || true

# excluded_final = candidates - final_all
comm -23 <(LC_ALL=C sort -u "${D}/${F}_candidates.id") <(LC_ALL=C sort -u "${D}/${F}_final_all.id") > "${D}/${F}_excluded_final.id"

# evidence for excluded_final
grep -Ff "${D}/${F}_excluded_final.id" "${EVD}/${F}_blast_top1.tsv" > "${D}/${F}_excluded_final.top1.tsv" || true

# excluded.tsv（粗标签，NP-safe）
awk '{print $1"\tEXCLUDED_BY_RULES_OR_NOHIT"}' "${D}/${F}_excluded_final.id" > "${D}/${F}_excluded.tsv"

{
  echo "UBQ QC (run_UBQ)"
  echo "Candidates (PF00240): $(wc -l < ${D}/${F}_candidates.id)"
  echo "poly_merged (repeat>=2): $(wc -l < ${D}/${F}_poly_merged.id)"
  echo "repeat1: $(wc -l < ${D}/${F}_repeat1.id)"
  echo "repeat1 keep (POS_RE): $(wc -l < ${D}/${F}_keep_ribofusion.id)"
  echo "repeat1 remove (NEG_RE): $(wc -l < ${D}/${F}_remove_by_blast.id)"
  echo "NoHit: $(wc -l < ${D}/${F}_nohit.id)"
  echo "FINAL strict: $(wc -l < ${D}/${F}_final_all.id)"
} > "${D}/${F}_qc.txt"

echo "[DONE] outputs in ${D}"
