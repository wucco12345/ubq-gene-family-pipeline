#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh

F="ACT"
D="results/${F}"
PROTEOME="inputs/proteome.fa"
PFAM="results/pfam.tsv"
DB="db/uniprot/sprot"   # Swiss-Prot BLAST db prefix

mkdir -p "${D}/evidence"

# ---------- ACT-specific rules ----------
PF_RE='^PF00022([.][0-9]+)?'   # accept PF00022 or PF00022.xx
IE_COL=9
PROT_COL=3
# Length windows (you can tune, but keep consistent with your current run)
HIGH_MIN=300
HIGH_MAX=450
RESCUE_MIN=250
RESCUE_MAX=500

# ---------- 1) candidates from pfam ----------
log "[1/6] Pfam candidates: PF00022 (iE<=1e-5)"
pfam_pick_ids "$PFAM" "$PF_RE" "$IE_COL" "$PROT_COL" "${D}/${F}.id"
log "Unique proteins: $(wc -l < ${D}/${F}.id)"

# ---------- 2) extract candidate FASTA ----------
log "[2/6] Extract candidate FASTA"
extract_fa_by_id "${D}/${F}.id" "$PROTEOME" "${D}/${F}_candidates.fa"
log "Candidate FASTA seqs: $(grep -c '^>' ${D}/${F}_candidates.fa)"

# ---------- 3) length table + len filters ----------
log "[3/6] Length table + length filters"
# length table
awk '
  /^>/{if(seq){print id"\t"len}; id=substr($0,2); split(id,a,/[ \t]/); id=a[1]; len=0; seq=1; next}
  {gsub(/\r/,""); len+=length($0)}
  END{if(seq){print id"\t"len}}
' "${D}/${F}_candidates.fa" > "${D}/evidence/${F}_len.tsv"

# high/rescue IDs by length
awk -v mn="$HIGH_MIN" -v mx="$HIGH_MAX" '$2>=mn && $2<=mx{print $1}' "${D}/evidence/${F}_len.tsv" > "${D}/${F}_high_len.id"
awk -v mn="$RESCUE_MIN" -v mx="$RESCUE_MAX" '$2>=mn && $2<=mx{print $1}' "${D}/evidence/${F}_len.tsv" > "${D}/${F}_rescue_len_all.id"

# Extract length-filtered fasta (optional, but keeps outputs similar to current)
extract_fa_by_id "${D}/${F}_high_len.id" "${D}/${F}_candidates.fa" "${D}/${F}_high_len.fa" || true
extract_fa_by_id "${D}/${F}_rescue_len_all.id" "${D}/${F}_candidates.fa" "${D}/${F}_rescue_len.fa" || true

# ---------- 4) BLASTp candidates vs Swiss-Prot (top5 + top1) ----------
log "[4/6] BLASTp candidates vs Swiss-Prot (top5; derive top1)"
### BLAST replaced by blast_top5 + blast_top1_from_top5
blast_top5 "${D}/${F}_candidates.fa" "$DB" "${D}/evidence/${F}_vs_sprot.tsv" 8 blastp
blast_top1_from_top5 "${D}/evidence/${F}_vs_sprot.tsv" "${D}/evidence/${F}_blast_top1.tsv" "${D}/evidence/${F}_top1.id"



# ---------- 5) Keep/remove by BLAST top1 title (loose) ----------
# NOTE: final ACT/ARP split is done by post_ACT_standardize.v2.sh (top5 consensus)
log "[5/6] Keep/remove by BLAST top1 annotation (loose; final split in post script)"
awk -F'\t' 'BEGIN{IGNORECASE=1}
{
  q=$1; t=$7
  if(t ~ /actin/){print q > "'"${D}/${F}_keep_by_blast.id"'"}
  else{print q > "'"${D}/${F}_remove_by_blast.id"'"}
}' "${D}/evidence/${F}_blast_top1.tsv"

# Build preliminary high_conf/rescue sets (kept by blast + length windows)
LC_ALL=C sort -u "${D}/${F}_keep_by_blast.id" > "${D}/${F}_keep_by_blast.sorted.id"
LC_ALL=C sort -u "${D}/${F}_high_len.id" > "${D}/${F}_high_len.sorted.id"
LC_ALL=C sort -u "${D}/${F}_rescue_len_all.id" > "${D}/${F}_rescue_len_all.sorted.id"

comm -12 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_len.sorted.id" > "${D}/${F}_high_conf.id"
comm -12 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_rescue_len_all.sorted.id" > "${D}/${F}_rescue.id"

extract_fa_by_id "${D}/${F}_high_conf.id" "${D}/${F}_candidates.fa" "${D}/${F}_high_conf.fa" || true
extract_fa_by_id "${D}/${F}_rescue.id" "${D}/${F}_candidates.fa" "${D}/${F}_rescue.fa" || true

# ---------- 6) QC ----------
{
  echo "ACT QC (run_ACT)"
  echo "Candidates (PF00022, iE<=1e-5): $(wc -l < ${D}/${F}.id)"
  echo "High_len (${HIGH_MIN}-${HIGH_MAX}): $(wc -l < ${D}/${F}_high_len.id)"
  echo "Rescue_len (${RESCUE_MIN}-${RESCUE_MAX}): $(wc -l < ${D}/${F}_rescue_len_all.id)"
  echo "Keep_by_blast(actin*): $(wc -l < ${D}/${F}_keep_by_blast.id)"
  echo "FINAL high_conf: $(wc -l < ${D}/${F}_high_conf.id)"
  echo "FINAL rescue: $(wc -l < ${D}/${F}_rescue.id)"
} > "${D}/${F}_qc.txt"

echo "[DONE] outputs in ${D}"
