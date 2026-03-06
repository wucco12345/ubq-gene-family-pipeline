#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
source scripts/lib/yaml.sh

F="ACT"
CFG="config/ACT.yml"
D="results/${F}"
PROTEOME="inputs/proteome.fa"
PFAM="results/pfam.tsv"
DB="db/uniprot/sprot"

mkdir -p "${D}/evidence"

# ---- read config ----
PF_RE="$(yget "$CFG" pfam.pf_re)"
IE_COL="$(yget "$CFG" pfam.ie_col 9)"
PROT_COL="$(yget "$CFG" pfam.protein_col 3)"

HIGH_MIN="$(yget "$CFG" length.high_min 300)"
HIGH_MAX="$(yget "$CFG" length.high_max 450)"
RESCUE_MIN="$(yget "$CFG" length.rescue_min 250)"
RESCUE_MAX="$(yget "$CFG" length.rescue_max 500)"

BLAST_TASK="$(yget "$CFG" blast.task blastp)"
BLAST_EVALUE="$(yget "$CFG" blast.evalue 1e-10)"

log "[1/6] Pfam candidates: PF00022 (iE<=1e-5)"
pfam_pick_ids "$PFAM" "$PF_RE" "$IE_COL" "$PROT_COL" "${D}/${F}.id"
log "Unique proteins: $(wc -l < ${D}/${F}.id)"

log "[2/6] Extract candidate FASTA"
extract_fa_by_id "${D}/${F}.id" "$PROTEOME" "${D}/${F}_candidates.fa"
log "Candidate FASTA seqs: $(grep -c '^>' ${D}/${F}_candidates.fa)"

log "[3/6] Length table + length filters"
fa_len_tsv "${D}/${F}_candidates.fa" "${D}/evidence/${F}_len.tsv"
len_filter_ids "${D}/evidence/${F}_len.tsv" "$HIGH_MIN" "$HIGH_MAX" "${D}/${F}_high_len.id"
len_filter_ids "${D}/evidence/${F}_len.tsv" "$RESCUE_MIN" "$RESCUE_MAX" "${D}/${F}_rescue_len_all.id"

extract_fa_by_id "${D}/${F}_high_len.id" "${D}/${F}_candidates.fa" "${D}/${F}_high_len.fa" || true
extract_fa_by_id "${D}/${F}_rescue_len_all.id" "${D}/${F}_candidates.fa" "${D}/${F}_rescue_len.fa" || true

log "[4/6] BLASTp candidates vs Swiss-Prot (top5 -> top1)"
blast_top5 "${D}/${F}_candidates.fa" "$DB" "${D}/evidence/${F}_vs_sprot.tsv" "$THREADS_DEFAULT" "$BLAST_TASK" "$BLAST_EVALUE"
blast_top1_from_top5 "${D}/evidence/${F}_vs_sprot.tsv" "${D}/evidence/${F}_blast_top1.tsv" "${D}/evidence/${F}_top1.id"

log "[5/6] Keep/remove by BLAST top1 annotation (loose; final split in post script)"
awk -F'\t' 'BEGIN{IGNORECASE=1}
{
  q=$1; t=$7
  if(t ~ /actin/){print q > "'"${D}/${F}_keep_by_blast.id"'"}
  else{print q > "'"${D}/${F}_remove_by_blast.id"'"}
}' "${D}/evidence/${F}_blast_top1.tsv"

LC_ALL=C sort -u "${D}/${F}_keep_by_blast.id" > "${D}/${F}_keep_by_blast.sorted.id"
LC_ALL=C sort -u "${D}/${F}_high_len.id" > "${D}/${F}_high_len.sorted.id"
LC_ALL=C sort -u "${D}/${F}_rescue_len_all.id" > "${D}/${F}_rescue_len_all.sorted.id"

comm -12 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_len.sorted.id" > "${D}/${F}_high_conf.id"
comm -12 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_rescue_len_all.sorted.id" > "${D}/${F}_rescue.id"

extract_fa_by_id "${D}/${F}_high_conf.id" "${D}/${F}_candidates.fa" "${D}/${F}_high_conf.fa" || true
extract_fa_by_id "${D}/${F}_rescue.id" "${D}/${F}_candidates.fa" "${D}/${F}_rescue.fa" || true

log "[6/6] QC"
{
  echo "ACT QC (run_ACT)"
  echo "Candidates (PF00022): $(wc -l < ${D}/${F}.id)"
  echo "High_len (${HIGH_MIN}-${HIGH_MAX}): $(wc -l < ${D}/${F}_high_len.id)"
  echo "Rescue_len (${RESCUE_MIN}-${RESCUE_MAX}): $(wc -l < ${D}/${F}_rescue_len_all.id)"
  echo "Keep_by_blast(actin*): $(wc -l < ${D}/${F}_keep_by_blast.id)"
  echo "FINAL high_conf: $(wc -l < ${D}/${F}_high_conf.id)"
  echo "FINAL rescue: $(wc -l < ${D}/${F}_rescue.id)"
} > "${D}/${F}_qc.txt"

echo "[DONE] outputs in ${D}"
