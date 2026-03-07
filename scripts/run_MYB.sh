#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
source scripts/lib/yaml.sh

F="MYB"
CFG="config/MYB.yml"
D="results/${F}"
EVD="${D}/evidence"
PROTEOME="inputs/proteome.fa"
PFAM="results/pfam.tsv"
DB="db/uniprot/sprot"

mkdir -p "$D" "$EVD"

PF_RE="$(yget "$CFG" pfam.pf_re)"
IE_COL="$(yget "$CFG" pfam.ie_col 9)"
PROT_COL="$(yget "$CFG" pfam.protein_col 3)"
ALI_FROM_COL="$(yget "$CFG" pfam.ali_from_col 12)"
ALI_TO_COL="$(yget "$CFG" pfam.ali_to_col 13)"

HIGH_MIN="$(yget "$CFG" length.high_min 250)"
HIGH_MAX="$(yget "$CFG" length.high_max 450)"
RESCUE_MIN="$(yget "$CFG" length.rescue_min 200)"
RESCUE_MAX="$(yget "$CFG" length.rescue_max 600)"

BLAST_TASK="$(yget "$CFG" blast.task blastp)"
BLAST_EVALUE="$(yget "$CFG" blast.evalue 1e-10)"

KEEP_RE="$(yget "$CFG" rules.keep_top1_re)"
REMOVE_RE="$(yget "$CFG" rules.hard_remove_top1_re)"

log "[1/10] Pfam candidates: PF00249 (iE<=1e-5)"
pfam_pick_ids "$PFAM" "$PF_RE" "$IE_COL" "$PROT_COL" "${D}/${F}_candidates.id"
cp "${D}/${F}_candidates.id" "${D}/${F}_candidates.raw.id"
log "Candidates: $(wc -l < "${D}/${F}_candidates.id")"

log "[2/10] Extract candidate FASTA"
extract_fa_by_id "${D}/${F}_candidates.id" "$PROTEOME" "${D}/${F}_candidates.fa"
log "Candidate FASTA seqs: $(grep -c '^>' "${D}/${F}_candidates.fa")"

log "[3/10] Save PFam hit intervals"
awk -F'\t' -v re="$PF_RE" -v IEC="$IE_COL" -v PC="$PROT_COL" -v AF="$ALI_FROM_COL" -v AT="$ALI_TO_COL" '
BEGIN{OFS="\t"}
{
  gsub(/\r/,"",$0)
  ie=$(IEC)+0
  if(ie<=1e-5 && $2 ~ re){
    id=$(PC)
    gsub(/\r/,"",id)
    sub(/[ \t]+$/,"",id)
    s=$(AF)+0; e=$(AT)+0
    if(s>e){t=s;s=e;e=t}
    print id, s, e, ie, $2
  }
}' "$PFAM" | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "${EVD}/${F}_pfam_hits.tsv"

log "[4/10] Merge intervals and count repeats"
awk '
BEGIN{OFS="\t"}
{
  id=$1; s=$2+0; e=$3+0
  if(id!=cur && cur!=""){
    print cur, seg, min_from, max_to
  }
  if(id!=cur){
    cur=id
    seg=1
    ms=s; me=e
    min_from=s
    max_to=e
  } else {
    if(s <= me+1){
      if(e>me) me=e
    } else {
      seg++
      ms=s; me=e
    }
    if(s<min_from) min_from=s
    if(e>max_to) max_to=e
  }
}
END{
  if(cur!="") print cur, seg, min_from, max_to
}' "${EVD}/${F}_pfam_hits.tsv" > "${D}/${F}_repeat_merged_counts.tsv"

awk '$2==2{print $1}' "${D}/${F}_repeat_merged_counts.tsv" | LC_ALL=C sort -u > "${D}/${F}_repeat2.id"
awk '$2==1{print $1}' "${D}/${F}_repeat_merged_counts.tsv" | LC_ALL=C sort -u > "${D}/${F}_repeat1.id"
awk '$2>=3{print $1}' "${D}/${F}_repeat_merged_counts.tsv" | LC_ALL=C sort -u > "${D}/${F}_repeat3plus.id"

log "repeat=2: $(wc -l < "${D}/${F}_repeat2.id")"
log "repeat=1: $(wc -l < "${D}/${F}_repeat1.id")"
log "repeat>=3: $(wc -l < "${D}/${F}_repeat3plus.id")"

log "[5/10] Length table + QC bins for repeat=2"
extract_fa_by_id "${D}/${F}_repeat2.id" "$PROTEOME" "${D}/${F}_repeat2.fa"
fa_len_tsv "${D}/${F}_repeat2.fa" "${EVD}/${F}_len.tsv"
len_filter_ids "${EVD}/${F}_len.tsv" "$HIGH_MIN" "$HIGH_MAX" "${D}/${F}_high_len.id"
len_filter_ids "${EVD}/${F}_len.tsv" "$RESCUE_MIN" "$RESCUE_MAX" "${D}/${F}_rescue_len.id"

log "[6/10] BLASTp repeat=2 candidates vs Swiss-Prot (top5 -> top1)"
blast_top5 "${D}/${F}_repeat2.fa" "$DB" "${EVD}/${F}_vs_sprot.tsv" "$THREADS_DEFAULT" "$BLAST_TASK" "$BLAST_EVALUE"
blast_top1_from_top5 "${EVD}/${F}_vs_sprot.tsv" "${EVD}/${F}_blast_top1.tsv" "${D}/${F}_top1.id"

log "[7/10] nohit = repeat2 - top1.id"
comm -23 <(LC_ALL=C sort -u "${D}/${F}_repeat2.id") <(LC_ALL=C sort -u "${D}/${F}_top1.id") > "${D}/${F}_nohit.id"

log "[8/10] Classify by top1 annotation"
: > "${D}/${F}_keep_by_blast.id"
: > "${D}/${F}_remove_by_blast.id"
: > "${D}/${F}_ambiguous_by_blast.id"

awk -F'\t' -v KEEP_RE="$KEEP_RE" -v REMOVE_RE="$REMOVE_RE" 'BEGIN{IGNORECASE=1}
{
  q=$1; t=$7
  if(t ~ REMOVE_RE){print q >> "'"${D}/${F}_remove_by_blast.id"'"; next}
  if(t ~ KEEP_RE){print q >> "'"${D}/${F}_keep_by_blast.id"'"; next}
  print q >> "'"${D}/${F}_ambiguous_by_blast.id"'"
}' "${EVD}/${F}_blast_top1.tsv"

LC_ALL=C sort -u "${D}/${F}_keep_by_blast.id" > "${D}/${F}_keep_by_blast.sorted.id"
LC_ALL=C sort -u "${D}/${F}_high_len.id" > "${D}/${F}_high_len.sorted.id"
LC_ALL=C sort -u "${D}/${F}_rescue_len.id" > "${D}/${F}_rescue_len.sorted.id"

log "[9/10] Build classic R2R3-MYB sets"
# 主集：keep_by_blast ∩ high_len
comm -12 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_len.sorted.id" > "${D}/${F}_high_conf.id"

# 边界集 seed：
#   1) keep_by_blast 但不在 high_conf
#   2) ambiguous_by_blast 且落在 rescue_len
comm -23 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_conf.id" > "${D}/${F}_strong_not_high.id"
comm -12 <(LC_ALL=C sort -u "${D}/${F}_ambiguous_by_blast.id") "${D}/${F}_rescue_len.sorted.id" > "${D}/${F}_ambiguous_rescue.id"

cat "${D}/${F}_strong_not_high.id" "${D}/${F}_ambiguous_rescue.id" 2>/dev/null | LC_ALL=C sort -u > "${D}/${F}_rescue.id"

log "[10/10] QC"
{
  echo "MYB QC (run_MYB classic R2R3)"
  echo "Candidates (PF00249): $(wc -l < "${D}/${F}_candidates.id")"
  echo "repeat=2: $(wc -l < "${D}/${F}_repeat2.id")"
  echo "repeat=1: $(wc -l < "${D}/${F}_repeat1.id")"
  echo "repeat>=3: $(wc -l < "${D}/${F}_repeat3plus.id")"
  echo "High_len (${HIGH_MIN}-${HIGH_MAX}): $(wc -l < "${D}/${F}_high_len.id")"
  echo "Rescue_len (${RESCUE_MIN}-${RESCUE_MAX}): $(wc -l < "${D}/${F}_rescue_len.id")"
  echo "Keep_by_blast: $(wc -l < "${D}/${F}_keep_by_blast.id")"
  echo "Remove_by_blast: $(wc -l < "${D}/${F}_remove_by_blast.id")"
  echo "Ambiguous_by_blast: $(wc -l < "${D}/${F}_ambiguous_by_blast.id")"
  echo "NoHit: $(wc -l < "${D}/${F}_nohit.id")"
  echo "High_conf: $(wc -l < "${D}/${F}_high_conf.id")"
  echo "Rescue: $(wc -l < "${D}/${F}_rescue.id")"
} > "${D}/${F}_qc.txt"

echo "[DONE] outputs in ${D}"
