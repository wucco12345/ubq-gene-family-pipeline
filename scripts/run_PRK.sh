#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
source scripts/lib/yaml.sh

F="PRK"
CFG="config/PRK.yml"
D="results/${F}"
EVD="${D}/evidence"
PROTEOME="inputs/proteome.fa"
PFAM="results/pfam.tsv"
DB="db/uniprot/sprot"

mkdir -p "$D" "$EVD"

PF_RE="$(yget "$CFG" pfam.pf_re)"
IE_COL="$(yget "$CFG" pfam.ie_col 9)"
PROT_COL="$(yget "$CFG" pfam.protein_col 3)"

HIGH_MIN="$(yget "$CFG" length.high_min 250)"
HIGH_MAX="$(yget "$CFG" length.high_max 450)"
RESCUE_MIN="$(yget "$CFG" length.rescue_min 200)"
RESCUE_MAX="$(yget "$CFG" length.rescue_max 500)"

BLAST_TASK="$(yget "$CFG" blast.task $BLAST_TASK_DEFAULT)"
BLAST_EVALUE="$(yget "$CFG" blast.evalue $BLAST_EVALUE_DEFAULT)"

KEEP_RE="$(yget "$CFG" rules.keep_top1_re)"
REMOVE_RE="$(yget "$CFG" rules.remove_top1_re)"

log "[1/7] Pfam candidates: PF00485 (iE<=1e-5)"
pfam_pick_ids "$PFAM" "$PF_RE" "$IE_COL" "$PROT_COL" "${D}/${F}_candidates.id"
log "Candidates: $(wc -l < "${D}/${F}_candidates.id")"

log "[2/7] Extract candidate FASTA"
extract_fa_by_id "${D}/${F}_candidates.id" "$PROTEOME" "${D}/${F}_candidates.fa"
log "Candidate FASTA seqs: $(grep -c '^>' "${D}/${F}_candidates.fa")"

log "[3/7] Length table + length bins (QC only, not BLAST gate)"
fa_len_tsv "${D}/${F}_candidates.fa" "${EVD}/${F}_len.tsv"
len_filter_ids "${EVD}/${F}_len.tsv" "$HIGH_MIN" "$HIGH_MAX" "${D}/${F}_high_len.id"
len_filter_ids "${EVD}/${F}_len.tsv" "$RESCUE_MIN" "$RESCUE_MAX" "${D}/${F}_rescue_len.id"

log "[4/7] BLASTp ALL candidates vs Swiss-Prot (top5 -> top1)"
blast_top5 "${D}/${F}_candidates.fa" "$DB" "${EVD}/${F}_vs_sprot.tsv" "$THREADS_DEFAULT" "$BLAST_TASK"
blast_top1_from_top5 "${EVD}/${F}_vs_sprot.tsv" "${EVD}/${F}_blast_top1.tsv" "${D}/${F}_top1.id"

log "[5/7] nohit = candidates - top1.id"
comm -23 <(LC_ALL=C sort -u "${D}/${F}_candidates.id") <(LC_ALL=C sort -u "${D}/${F}_top1.id") > "${D}/${F}_nohit.id"

log "[6/7] Classify by BLAST top1"
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
LC_ALL=C sort -u "${D}/${F}_candidates.id" > "${D}/${F}_candidates.sorted.id"

# high_conf = BLAST membership + high_len
comm -12 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_len.sorted.id" > "${D}/${F}_high_conf.id"

# rescue = BLAST membership but NOT high_conf
# 这里不再要求一定落在 rescue_len；只要 BLAST 支持但长度不在 high_len，就先放 rescue
comm -23 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_conf.id" > "${D}/${F}_rescue.id"

log "[7/7] Run-level QC"
{
  echo "PRK QC (run_PRK all-candidate-blast)"
  echo "Candidates (PF00485): $(wc -l < "${D}/${F}_candidates.id")"
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
