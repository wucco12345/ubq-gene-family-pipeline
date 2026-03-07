#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
source scripts/lib/yaml.sh

F="PRR"
CFG="config/PRR.yml"
D="results/${F}"
EVD="${D}/evidence"
PROTEOME="inputs/proteome.fa"
PFAM="results/pfam.tsv"
DB="db/uniprot/sprot"

mkdir -p "$D" "$EVD"

PF00072_RE="$(yget "$CFG" pfam.pf00072_re)"
PF06203_RE="$(yget "$CFG" pfam.pf06203_re)"
IE_COL="$(yget "$CFG" pfam.ie_col 9)"
PROT_COL="$(yget "$CFG" pfam.protein_col 3)"
ALI_FROM_COL="$(yget "$CFG" pfam.ali_from_col 12)"
ALI_TO_COL="$(yget "$CFG" pfam.ali_to_col 13)"

HIGH_MIN="$(yget "$CFG" length.high_min 200)"
HIGH_MAX="$(yget "$CFG" length.high_max 500)"
RESCUE_MIN="$(yget "$CFG" length.rescue_min 150)"
RESCUE_MAX="$(yget "$CFG" length.rescue_max 700)"

BLAST_TASK="$(yget "$CFG" blast.task $BLAST_TASK_DEFAULT)"
BLAST_EVALUE="$(yget "$CFG" blast.evalue $BLAST_EVALUE_DEFAULT)"

KEEP_RE="$(yget "$CFG" rules.keep_top1_re)"
REMOVE_RE="$(yget "$CFG" rules.remove_top1_re)"

log "[1/8] Pfam candidates: PF00072 + PF06203 (dual-domain required)"
awk -F'\t' \
  -v re1="$PF00072_RE" \
  -v re2="$PF06203_RE" \
  -v IEC="$IE_COL" \
  -v PC="$PROT_COL" '
{
  gsub(/\r/,"",$0)
  ie=$(IEC)+0
  if(ie<=1e-5){
    id=$(PC)
    gsub(/\r/,"",id)
    sub(/[ \t]+$/,"",id)
    if($2 ~ re1) a[id]=1
    if($2 ~ re2) b[id]=1
  }
}
END{
  for(id in a){
    if(id in b) print id
  }
}' "$PFAM" | LC_ALL=C sort -u > "${D}/${F}_candidates.id"
log "Candidates (dual-domain): $(wc -l < "${D}/${F}_candidates.id")"

log "[2/8] Extract candidate FASTA"
extract_fa_by_id "${D}/${F}_candidates.id" "$PROTEOME" "${D}/${F}_candidates.fa"
log "Candidate FASTA seqs: $(grep -c '^>' "${D}/${F}_candidates.fa")"

log "[3/8] Domain status table"
awk -F'\t' \
  -v re1="$PF00072_RE" \
  -v re2="$PF06203_RE" \
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
      pf00072[id]=1
      s=$(AFC)+0; e=$(ATC)+0
      if(!(id in min_from) || s<min_from[id]) min_from[id]=s
      if(!(id in max_to) || e>max_to[id]) max_to[id]=e
    }
    if($2 ~ re2){
      pf06203[id]=1
      s=$(AFC)+0; e=$(ATC)+0
      if(!(id in min_from) || s<min_from[id]) min_from[id]=s
      if(!(id in max_to) || e>max_to[id]) max_to[id]=e
    }
  }
}
END{
  for(id in min_from){
    has1=(id in pf00072 ? 1 : 0)
    has2=(id in pf06203 ? 1 : 0)
    cls="other"
    if(has1 && has2) cls="dual_domain"
    else if(has1) cls="pf00072_only"
    else if(has2) cls="pf06203_only"
    print id, has1, has2, cls, min_from[id], max_to[id]
  }
}' "$PFAM" | LC_ALL=C sort -k1,1 > "${D}/${F}_domain_status.tsv"

awk -F'\t' '$4=="dual_domain"{print $1}' "${D}/${F}_domain_status.tsv" | LC_ALL=C sort -u > "${D}/${F}_dual_domain.id"

log "dual_domain: $(wc -l < "${D}/${F}_dual_domain.id")"

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
  if(t ~ KEEP_RE){print q >> "'"${D}/${F}_keep_by_blast.id"'"; next}
  if(t ~ REMOVE_RE){print q >> "'"${D}/${F}_remove_by_blast.id"'"; next}
  print q >> "'"${D}/${F}_ambiguous_by_blast.id"'"
}' "${EVD}/${F}_blast_top1.tsv"

LC_ALL=C sort -u "${D}/${F}_keep_by_blast.id" > "${D}/${F}_keep_by_blast.sorted.id"
LC_ALL=C sort -u "${D}/${F}_high_len.id" > "${D}/${F}_high_len.sorted.id"

comm -12 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_len.sorted.id" > "${D}/${F}_high_conf.id"
comm -23 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_conf.id" > "${D}/${F}_rescue.id"

log "[8/8] QC"
{
  echo "PRR QC (run_PRR first-pass)"
  echo "Candidates (PF00072 + PF06203): $(wc -l < "${D}/${F}_candidates.id")"
  echo "Dual_domain: $(wc -l < "${D}/${F}_dual_domain.id")"
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
