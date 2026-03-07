#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
source scripts/lib/yaml.sh

F="bZIP"
CFG="config/bZIP.yml"
D="results/${F}"
EVD="${D}/evidence"
PROTEOME="inputs/proteome.fa"
PFAM="results/pfam.tsv"
DB="db/uniprot/sprot"

mkdir -p "$D" "$EVD"

PF00170_RE="$(yget "$CFG" pfam.pf00170_re)"
PF07716_RE="$(yget "$CFG" pfam.pf07716_re)"
IE_COL="$(yget "$CFG" pfam.ie_col 9)"
PROT_COL="$(yget "$CFG" pfam.protein_col 3)"
ALI_FROM_COL="$(yget "$CFG" pfam.ali_from_col 12)"
ALI_TO_COL="$(yget "$CFG" pfam.ali_to_col 13)"

HIGH_MIN="$(yget "$CFG" length.high_min 150)"
HIGH_MAX="$(yget "$CFG" length.high_max 800)"
RESCUE_MIN="$(yget "$CFG" length.rescue_min 1)"
RESCUE_MAX="$(yget "$CFG" length.rescue_max 2000)"

BLAST_TASK="$(yget "$CFG" blast.task $BLAST_TASK_DEFAULT)"
BLAST_EVALUE="$(yget "$CFG" blast.evalue $BLAST_EVALUE_DEFAULT)"

KEEP_RE="$(yget "$CFG" rules.keep_top1_re)"
REMOVE_RE="$(yget "$CFG" rules.remove_top1_re)"

log "[1/8] Pfam candidates: PF00170 OR PF07716"
awk -F'\t' -v re1="$PF00170_RE" -v re2="$PF07716_RE" -v IEC="$IE_COL" -v PC="$PROT_COL" '
{
  gsub(/\r/,"",$0)
  ie=$(IEC)+0
  if(ie<=1e-5 && ($2 ~ re1 || $2 ~ re2)){
    id=$(PC)
    gsub(/\r/,"",id)
    sub(/[ \t]+$/,"",id)
    print id
  }
}' "$PFAM" | LC_ALL=C sort -u > "${D}/${F}_candidates.id"
log "Candidates: $(wc -l < "${D}/${F}_candidates.id")"

log "[2/8] Extract candidate FASTA"
extract_fa_by_id "${D}/${F}_candidates.id" "$PROTEOME" "${D}/${F}_candidates.fa"
log "Candidate FASTA seqs: $(grep -c '^>' "${D}/${F}_candidates.fa")"

log "[3/8] Domain status table"
awk -F'\t' \
  -v re1="$PF00170_RE" \
  -v re2="$PF07716_RE" \
  -v IEC="$IE_COL" \
  -v PC="$PROT_COL" \
  -v AFC="$ALI_FROM_COL" \
  -v ATC="$ALI_TO_COL" '
BEGIN{OFS="\t"}
{
  gsub(/\r/,"",$0)
  ie=$(IEC)+0
  if(ie<=1e-5){
    id=$(PC)
    gsub(/\r/,"",id)
    sub(/[ \t]+$/,"",id)
    if($2 ~ re1){
      pf00170[id]=1
      s=$(AFC)+0; e=$(ATC)+0
      if(!(id in min_from) || s<min_from[id]) min_from[id]=s
      if(!(id in max_to) || e>max_to[id]) max_to[id]=e
    }
    if($2 ~ re2){
      pf07716[id]=1
      s=$(AFC)+0; e=$(ATC)+0
      if(!(id in min_from) || s<min_from[id]) min_from[id]=s
      if(!(id in max_to) || e>max_to[id]) max_to[id]=e
    }
  }
}
END{
  for(id in min_from){
    has1=(id in pf00170 ? 1 : 0)
    has2=(id in pf07716 ? 1 : 0)
    cls="other"
    if(has1 && has2) cls="both"
    else if(has1) cls="pf00170_only"
    else if(has2) cls="pf07716_only"
    print id, has1, has2, cls, min_from[id], max_to[id]
  }
}' "$PFAM" | LC_ALL=C sort -k1,1 > "${D}/${F}_domain_status.tsv"

awk -F'\t' '$4=="both"{print $1}' "${D}/${F}_domain_status.tsv" | LC_ALL=C sort -u > "${D}/${F}_both.id"
awk -F'\t' '$4=="pf00170_only"{print $1}' "${D}/${F}_domain_status.tsv" | LC_ALL=C sort -u > "${D}/${F}_pf00170_only.id"
awk -F'\t' '$4=="pf07716_only"{print $1}' "${D}/${F}_domain_status.tsv" | LC_ALL=C sort -u > "${D}/${F}_pf07716_only.id"
awk -F'\t' '$5+0<=150{print $1}' "${D}/${F}_domain_status.tsv" | LC_ALL=C sort -u > "${D}/${F}_domain_start_le150.id"

log "both: $(wc -l < "${D}/${F}_both.id")"
log "pf00170_only: $(wc -l < "${D}/${F}_pf00170_only.id")"
log "pf07716_only: $(wc -l < "${D}/${F}_pf07716_only.id")"

log "[4/8] Length table + QC bins"
fa_len_tsv "${D}/${F}_candidates.fa" "${EVD}/${F}_len.tsv"
len_filter_ids "${EVD}/${F}_len.tsv" "$HIGH_MIN" "$HIGH_MAX" "${D}/${F}_high_len.id"
len_filter_ids "${EVD}/${F}_len.tsv" "$RESCUE_MIN" "$RESCUE_MAX" "${D}/${F}_rescue_len.id"

log "[5/8] BLASTp ALL candidates vs Swiss-Prot (top5 -> top1)"
blast_top5 "${D}/${F}_candidates.fa" "$DB" "${EVD}/${F}_vs_sprot.tsv" "$THREADS_DEFAULT" "$BLAST_TASK"
blast_top1_from_top5 "${EVD}/${F}_vs_sprot.tsv" "${EVD}/${F}_blast_top1.tsv" "${D}/${F}_top1.id"

log "[6/8] nohit = candidates - top1.id"
comm -23 <(LC_ALL=C sort -u "${D}/${F}_candidates.id") <(LC_ALL=C sort -u "${D}/${F}_top1.id") > "${D}/${F}_nohit.id"

log "[7/8] Classify by top1 annotation"
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
LC_ALL=C sort -u "${D}/${F}_domain_start_le150.id" > "${D}/${F}_domain_start_le150.sorted.id"

# high_conf = keep_by_blast ∩ high_len ∩ domain_start<=150
comm -12 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_len.sorted.id"   | comm -12 - "${D}/${F}_domain_start_le150.sorted.id"   > "${D}/${F}_high_conf.id"

# rescue = keep_by_blast - high_conf
comm -23 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_conf.id" > "${D}/${F}_rescue.id"

log "[8/8] QC"
{
  echo "bZIP QC (run_bZIP first-pass)"
  echo "Candidates (PF00170 or PF07716): $(wc -l < "${D}/${F}_candidates.id")"
  echo "both: $(wc -l < "${D}/${F}_both.id")"
  echo "pf00170_only: $(wc -l < "${D}/${F}_pf00170_only.id")"
  echo "pf07716_only: $(wc -l < "${D}/${F}_pf07716_only.id")"
  echo "domain_start<=150: $(wc -l < "${D}/${F}_domain_start_le150.id")"
  echo "High_len (${HIGH_MIN}-${HIGH_MAX}): $(wc -l < "${D}/${F}_high_len.id")"
  echo "Keep_by_blast: $(wc -l < "${D}/${F}_keep_by_blast.id")"
  echo "Remove_by_blast: $(wc -l < "${D}/${F}_remove_by_blast.id")"
  echo "Ambiguous_by_blast: $(wc -l < "${D}/${F}_ambiguous_by_blast.id")"
  echo "NoHit: $(wc -l < "${D}/${F}_nohit.id")"
  echo "High_conf: $(wc -l < "${D}/${F}_high_conf.id")"
  echo "Rescue: $(wc -l < "${D}/${F}_rescue.id")"
} > "${D}/${F}_qc.txt"

echo "[DONE] outputs in ${D}"
