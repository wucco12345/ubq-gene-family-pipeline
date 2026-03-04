#!/usr/bin/env bash
set -euo pipefail

FAMILY="ACT"
PF="PF00022"
IE="1e-5"

# Inputs (local big files; not tracked by git)
PROTEOME="inputs/proteome.fa"
PFAM_TSV="results/pfam.tsv"

# Swiss-Prot BLAST DB prefix (change if your db is elsewhere)
# Examples:
#   SPROT_DB="db/uniprot/sprot"
#   SPROT_DB="$HOME/gene_family/db/uniprot/sprot"
SPROT_DB="db/uniprot/sprot"

# Output dirs
OUTDIR="results/${FAMILY}"
EVDIR="${OUTDIR}/evidence"
mkdir -p "${OUTDIR}" "${EVDIR}"

# Sanity checks
if [[ ! -f "${PROTEOME}" ]]; then
  echo "[ERROR] Missing ${PROTEOME} (put your proteome at inputs/proteome.fa)"
  exit 2
fi
if [[ ! -f "${PFAM_TSV}" ]]; then
  echo "[ERROR] Missing ${PFAM_TSV} (put parsed pfam table at results/pfam.tsv)"
  exit 2
fi
if [[ ! -f "${SPROT_DB}.pin" && ! -f "${SPROT_DB}.psq" ]]; then
  echo "[ERROR] Swiss-Prot BLAST DB not found at prefix: ${SPROT_DB}"
  echo "        Expected files like ${SPROT_DB}.pin/.psq/.phr"
  exit 2
fi

echo "[1/6] Pfam candidates: ${PF} (iE<=${IE})"
bash scripts/pfam_pipeline.sh list "${PFAM_TSV}" "${PF}" "${IE}" "${OUTDIR}/${FAMILY}"

# Save pfam hits as evidence
cp -f "${OUTDIR}/${FAMILY}.all_hits.tsv" "${EVDIR}/${FAMILY}_pfam_hits.tsv" || true
cp -f "${OUTDIR}/${FAMILY}.id" "${OUTDIR}/${FAMILY}_candidates.id" || true

echo "[2/6] Extract candidate FASTA"
bash scripts/pfam_pipeline.sh fasta "${PROTEOME}" "${OUTDIR}/${FAMILY}.id" "${OUTDIR}/${FAMILY}_candidates.fa"

echo "[3/6] Length table + length filters"
# length table (id \t length)
seqkit fx2tab -n -l "${OUTDIR}/${FAMILY}_candidates.fa" > "${EVDIR}/${FAMILY}_len.tsv"

# High-confidence length: 300-450
awk -F'\t' '$2>=300 && $2<=450 {print $1}' "${EVDIR}/${FAMILY}_len.tsv" | sort -u > "${OUTDIR}/${FAMILY}_high_len.id"

# Rescue length: 250-500 (and not in high_len)
awk -F'\t' '$2>=250 && $2<=500 {print $1}' "${EVDIR}/${FAMILY}_len.tsv" | sort -u > "${OUTDIR}/${FAMILY}_rescue_len_all.id"
comm -23 <(sort "${OUTDIR}/${FAMILY}_rescue_len_all.id") <(sort "${OUTDIR}/${FAMILY}_high_len.id") > "${OUTDIR}/${FAMILY}_rescue_len.id"

# Excluded by length
comm -23 <(sort "${OUTDIR}/${FAMILY}.id") <(sort "${OUTDIR}/${FAMILY}_rescue_len_all.id") > "${OUTDIR}/${FAMILY}_excluded_len.id"

# Extract FASTA for high/rescue length sets
bash scripts/pfam_pipeline.sh fasta "${PROTEOME}" "${OUTDIR}/${FAMILY}_high_len.id"   "${OUTDIR}/${FAMILY}_high_len.fa"
bash scripts/pfam_pipeline.sh fasta "${PROTEOME}" "${OUTDIR}/${FAMILY}_rescue_len.id" "${OUTDIR}/${FAMILY}_rescue_len.fa"

echo "[4/6] BLASTp candidates vs Swiss-Prot (top5; classify by top1)"
blastp \
  -query "${OUTDIR}/${FAMILY}_candidates.fa" \
  -db "${SPROT_DB}" \
  -evalue 1e-10 \
  -max_target_seqs 5 \
  -num_threads 10 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
  -out "${EVDIR}/${FAMILY}_vs_sprot.tsv"

# top1 per query
awk '!seen[$1]++{print}' "${EVDIR}/${FAMILY}_vs_sprot.tsv" > "${EVDIR}/${FAMILY}_blast_top1.tsv"
cut -f1 "${EVDIR}/${FAMILY}_blast_top1.tsv" | sort -u > "${EVDIR}/${FAMILY}_top1.id"

# nohit = candidates - top1
comm -23 <(sort "${OUTDIR}/${FAMILY}.id") <(sort "${EVDIR}/${FAMILY}_top1.id") > "${OUTDIR}/${FAMILY}_nohit.id"

echo "[5/6] Keep/remove by BLAST top1 annotation"
# keep if top1 contains actin or actin-related (ARP)
awk -F'\t' 'BEGIN{IGNORECASE=1}
  $7 ~ /(^| )actin( |$)/ || $7 ~ /actin-related/ || $7 ~ /actin related/ {print $1}
' "${EVDIR}/${FAMILY}_blast_top1.tsv" | sort -u > "${OUTDIR}/${FAMILY}_keep_by_blast.id"

# remove = top1 - keep
comm -23 <(sort "${EVDIR}/${FAMILY}_top1.id") <(sort "${OUTDIR}/${FAMILY}_keep_by_blast.id") > "${OUTDIR}/${FAMILY}_remove_by_blast.id"

# FINAL sets:
# high_conf = high_len ∩ keep_by_blast
comm -12 <(sort "${OUTDIR}/${FAMILY}_high_len.id") <(sort "${OUTDIR}/${FAMILY}_keep_by_blast.id") > "${OUTDIR}/${FAMILY}_high_conf.id"

# rescue = rescue_len ∩ keep_by_blast
comm -12 <(sort "${OUTDIR}/${FAMILY}_rescue_len.id") <(sort "${OUTDIR}/${FAMILY}_keep_by_blast.id") > "${OUTDIR}/${FAMILY}_rescue.id"

# Extract final FASTA
bash scripts/pfam_pipeline.sh fasta "${PROTEOME}" "${OUTDIR}/${FAMILY}_high_conf.id" "${OUTDIR}/${FAMILY}_high_conf.fa"
bash scripts/pfam_pipeline.sh fasta "${PROTEOME}" "${OUTDIR}/${FAMILY}_rescue.id"    "${OUTDIR}/${FAMILY}_rescue.fa"

echo "[6/6] Build excluded table (simple, NP-safe traceability)"
# Create excluded list with a reason tag
# - excluded_len
# - remove_by_blast
# - nohit
{
  awk '{print $1"\tLEN_OUTSIDE_250_500"}' "${OUTDIR}/${FAMILY}_excluded_len.id" 2>/dev/null || true
  awk '{print $1"\tBLAST_TOP1_NOT_ACTIN_OR_ARP"}' "${OUTDIR}/${FAMILY}_remove_by_blast.id" 2>/dev/null || true
  awk '{print $1"\tNO_BLAST_HIT"}' "${OUTDIR}/${FAMILY}_nohit.id" 2>/dev/null || true
} | sort -u > "${OUTDIR}/${FAMILY}_excluded.tsv"

# QC summary
{
  echo "${FAMILY} QC"
  echo "Candidates (PF00022, iE<=1e-5): $(wc -l < "${OUTDIR}/${FAMILY}.id")"
  echo "High_len (300-450): $(wc -l < "${OUTDIR}/${FAMILY}_high_len.id")"
  echo "Rescue_len (250-500, excl high): $(wc -l < "${OUTDIR}/${FAMILY}_rescue_len.id")"
  echo "BLAST keep (actin/ARP): $(wc -l < "${OUTDIR}/${FAMILY}_keep_by_blast.id")"
  echo "BLAST remove (top1 not actin/ARP): $(wc -l < "${OUTDIR}/${FAMILY}_remove_by_blast.id")"
  echo "NoHit: $(wc -l < "${OUTDIR}/${FAMILY}_nohit.id")"
  echo "FINAL high_conf: $(wc -l < "${OUTDIR}/${FAMILY}_high_conf.id")"
  echo "FINAL rescue: $(wc -l < "${OUTDIR}/${FAMILY}_rescue.id")"
  echo "Top1 head:"
  head -n 10 "${EVDIR}/${FAMILY}_blast_top1.tsv" || true
} > "${OUTDIR}/${FAMILY}_qc.txt"

echo "[DONE] outputs in ${OUTDIR}"
