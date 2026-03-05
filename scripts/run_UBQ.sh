#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib.sh"

# ===== Inputs (fixed names, same as EF1A/ACT convention) =====
PROTEOME="inputs/proteome.fa"
PFAM_TSV="results/pfam.tsv"
DBPREFIX="db/uniprot/sprot"

# ===== Params =====
PF="PF00240"
IE_THR="1e-5"

OUTDIR="results/UBQ"
EVDIR="${OUTDIR}/evidence"
mkdir -p "$EVDIR"

# Always normalize pfam.tsv in place to avoid CRLF issues
normalize_tsv_inplace "$PFAM_TSV"

echo "[1/7] Pfam candidates: $PF (iE<=${IE_THR})"
pfam_candidates_by_pf "$PF" "$IE_THR" "$PFAM_TSV" "${OUTDIR}/UBQ_candidates.raw.id"
clean_id_list "${OUTDIR}/UBQ_candidates.raw.id" "${OUTDIR}/UBQ_candidates.id"
info "Candidates: $(wc -l < ${OUTDIR}/UBQ_candidates.id)"

echo "[2/7] Extract candidate FASTA (exact header token match)"
extract_fa_by_id "${OUTDIR}/UBQ_candidates.id" "$PROTEOME" "${OUTDIR}/UBQ_candidates.fa"
info "Candidate FASTA seqs: $(grep -c '^>' ${OUTDIR}/UBQ_candidates.fa)"

echo "[3/7] Build merged-repeat counts from PFam hit intervals (ali_from-ali_to)"
# Save pfam hits for UBQ candidates to evidence
# pfam.tsv columns: (1)model (2)pfacc (3)target(protein) ... (12)ali_from (13)ali_to ...
# We keep only UBQ candidates and PF00240 rows with iE<=thr, and store ali coords.
awk -F'\t' -v PF="$PF" -v THR="$IE_THR" '
  BEGIN{OFS="\t"}
  {
    gsub(/\r/,"",$0);
    pfacc=$2; ie=$9+0; id=$3;
    if(pfacc ~ ("^" PF "(\\.|$)") && ie<=THR){
      gsub(/\r/,"",id); sub(/[[:space:]]+$/,"",id);
      print id, $12, $13, ie, $2
    }
  }' "$PFAM_TSV" | sort -k1,1 -k2,2n > "${EVDIR}/UBQ_pfam_hits.tsv"

# Merge overlapping/adjacent intervals per protein and count repeats
# Repeat count = number of merged intervals
awk '
  BEGIN{OFS="\t"}
  {
    id=$1; s=$2+0; e=$3+0;
    if(id!=cur && cur!=""){
      # flush previous
      print cur, n;
    }
    if(id!=cur){
      cur=id; n=0; ms=s; me=e;
    }else{
      # overlap or adjacent
      if(s <= me+1){
        if(e>me) me=e;
      }else{
        n++;
        ms=s; me=e;
      }
    }
    # we count by closing segments at flush: n counts completed segments, +1 for current
    n_curr= n+1;
    n = n; # keep n for completed segments
    n_last = n_curr;
  }
  END{
    if(cur!=""){ print cur, n_last }
  }' "${EVDIR}/UBQ_pfam_hits.tsv" > "${OUTDIR}/UBQ_repeat_merged_counts.tsv"

# Split poly (>=2) and repeat1 (=1)
awk '$2>=2{print $1}' "${OUTDIR}/UBQ_repeat_merged_counts.tsv" | sort -u > "${OUTDIR}/UBQ_poly_merged.id"
awk '$2==1{print $1}' "${OUTDIR}/UBQ_repeat_merged_counts.tsv" | sort -u > "${OUTDIR}/UBQ_repeat1.id"
info "poly_merged (repeat>=2): $(wc -l < ${OUTDIR}/UBQ_poly_merged.id)"
info "repeat1 (repeat==1): $(wc -l < ${OUTDIR}/UBQ_repeat1.id)"

echo "[4/7] Extract repeat1 FASTA for BLAST"
extract_fa_by_id "${OUTDIR}/UBQ_repeat1.id" "$PROTEOME" "${OUTDIR}/UBQ_repeat1.fa"
info "repeat1 FASTA seqs: $(grep -c '^>' ${OUTDIR}/UBQ_repeat1.fa)"

echo "[5/7] BLASTp repeat1 vs Swiss-Prot (top5; derive top1)"
# UBQ is short-ish sometimes; but repeat1 here can include fusions; use normal blastp
blastp_top5 "${OUTDIR}/UBQ_repeat1.fa" "$DBPREFIX" "${EVDIR}/UBQ_vs_sprot.tsv" 8 blastp
top1_from_top5 "${EVDIR}/UBQ_vs_sprot.tsv" "${EVDIR}/UBQ_blast_top1.tsv"

echo "[6/7] Classify repeat1 by BLAST top1 (strict keep = ubiquitin-ribosomal fusion); nohit handled"
# nohit: if some repeat1 not present in top1 table (should be rare)
cut -f1 "${EVDIR}/UBQ_blast_top1.tsv" | sort -u > "${OUTDIR}/UBQ_top1.id"
comm -23 <(sort -u "${OUTDIR}/UBQ_repeat1.id") "${OUTDIR}/UBQ_top1.id" > "${OUTDIR}/UBQ_nohit.id"

# Keep strict: top1 title contains "ubiquitin" and "ribosomal" and "fusion"
awk -F'\t' 'BEGIN{IGNORECASE=1}
  $7 ~ /ubiquitin/ && $7 ~ /ribosomal/ && $7 ~ /fusion/ {print $1}
' "${EVDIR}/UBQ_blast_top1.tsv" | sort -u > "${OUTDIR}/UBQ_repeat1_keep_strict.id"

# Remove by blast: obvious UBL / non-fusion ubiquitin-related proteins (NP-safe traceability)
awk -F'\t' 'BEGIN{IGNORECASE=1}
  # hard negatives
  $7 ~ /(sumo|nedd8|rub|rad23|dsk2|ubiquitin receptor|deubiquitin|ubp\b|ubiquitin-like)/ {print $1}
' "${EVDIR}/UBQ_blast_top1.tsv" | sort -u > "${OUTDIR}/UBQ_remove_by_blast.id"

# Final strict UBQ = polyubiquitin (repeat>=2) + ribosomal fusion keep
cat "${OUTDIR}/UBQ_poly_merged.id" "${OUTDIR}/UBQ_repeat1_keep_strict.id" | sort -u > "${OUTDIR}/UBQ_high_conf.id"

# Build final_all (same as high_conf; UBQ rescue left empty but file exists for consistency)
: > "${OUTDIR}/UBQ_rescue.id"
cp "${OUTDIR}/UBQ_high_conf.id" "${OUTDIR}/UBQ_final_all.id"

echo "[7/7] Extract final FASTA + excluded_final set"
extract_fa_by_id "${OUTDIR}/UBQ_final_all.id" "$PROTEOME" "${OUTDIR}/UBQ_final_all.fa"
extract_fa_by_id "${OUTDIR}/UBQ_high_conf.id" "$PROTEOME" "${OUTDIR}/UBQ_high_conf.fa"

# excluded_final = candidates - final_all
comm -23 <(sort -u "${OUTDIR}/UBQ_candidates.id") <(sort -u "${OUTDIR}/UBQ_final_all.id") > "${OUTDIR}/UBQ_excluded_final.id"

# excluded.tsv: simple reason table
# poly and keep_strict are retained; others in excluded_final. We tag coarse reason.
awk 'BEGIN{FS=OFS="\t"}
  FNR==NR{keep[$1]=1; next}
  {if(!keep[$1]) print $1}
' "${OUTDIR}/UBQ_final_all.id" <(cat "${OUTDIR}/UBQ_candidates.id") \
  | awk '{print $1"\tEXCLUDED_BY_RULES_OR_NOHIT"}' > "${OUTDIR}/UBQ_excluded.tsv"

# QC summary
{
  echo "UBQ QC"
  echo "Candidates (PF00240, iE<=${IE_THR}): $(wc -l < ${OUTDIR}/UBQ_candidates.id)"
  echo "poly_merged (repeat>=2): $(wc -l < ${OUTDIR}/UBQ_poly_merged.id)"
  echo "repeat1: $(wc -l < ${OUTDIR}/UBQ_repeat1.id)"
  echo "repeat1 keep_strict (ribo fusion): $(wc -l < ${OUTDIR}/UBQ_repeat1_keep_strict.id)"
  echo "NoHit: $(wc -l < ${OUTDIR}/UBQ_nohit.id)"
  echo "FINAL strict: $(wc -l < ${OUTDIR}/UBQ_final_all.id)"
} > "${OUTDIR}/UBQ_qc.txt"

ok "outputs in ${OUTDIR}"
