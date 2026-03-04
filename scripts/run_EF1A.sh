#!/usr/bin/env bash
set -euo pipefail

PROTEOME="inputs/proteome.fa"
PFAM_TSV="results/pfam.tsv"
SPROT_DB="$HOME/gene_family/db/uniprot/sprot"   # reuse existing Swiss-Prot DB
OUTDIR="results/EF1A"
mkdir -p "$OUTDIR"/evidence

# 1) Pfam候选：PF00009，按 iE<=1e-5
awk -F'\t' 'BEGIN{OFS="\t"} $2 ~ /^PF00009(\.|$)/ && $9+0 <= 1e-5 {print}' "$PFAM_TSV" > "$OUTDIR/evidence/EF1A_pfam_hits.tsv"
awk -F'\t' '{print $3}' "$OUTDIR/evidence/EF1A_pfam_hits.tsv" | sort -u > "$OUTDIR/EF1A_candidates.id"

# 2) 提取候选FASTA
bash scripts/pfam_pipeline.sh fasta "$PROTEOME" "$OUTDIR/EF1A_candidates.id" "$OUTDIR/EF1A_candidates.fa"

# 3) 长度分层：high 420-480；rescue 400-520；其余excluded
seqkit fx2tab -n -l "$OUTDIR/EF1A_candidates.fa" > "$OUTDIR/evidence/EF1A_len.tsv"

awk '$2>=420 && $2<=480{print $1}' "$OUTDIR/evidence/EF1A_len.tsv" | sort -u > "$OUTDIR/EF1A_high_len.id"
awk '$2>=400 && $2<=520{print $1}' "$OUTDIR/evidence/EF1A_len.tsv" | sort -u > "$OUTDIR/EF1A_rescue_len.id"
comm -23 <(sort "$OUTDIR/EF1A_candidates.id") <(sort "$OUTDIR/EF1A_rescue_len.id") > "$OUTDIR/EF1A_excluded_len.id"

# 4) 对 rescue_len 集合做BLAST（覆盖 high+rescue 的总集合）
bash scripts/pfam_pipeline.sh fasta "$PROTEOME" "$OUTDIR/EF1A_rescue_len.id" "$OUTDIR/EF1A_rescue_len.fa"

blastp \
  -query "$OUTDIR/EF1A_rescue_len.fa" \
  -db "$SPROT_DB" \
  -evalue 1e-10 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
  -max_target_seqs 5 \
  -num_threads 10 \
  -out "$OUTDIR/evidence/EF1A_vs_sprot.tsv"

# 5) 取top1
awk '!seen[$1]++{print}' "$OUTDIR/evidence/EF1A_vs_sprot.tsv" > "$OUTDIR/evidence/EF1A_blast_top1.tsv"

# 6) 分类：keep 必须是 eEF1A；remove 必须排 EF-Tu/TUFM 等
awk -F'\t' 'BEGIN{IGNORECASE=1}
{
  t=$7
  if (t ~ /elongation factor 1-alpha/ || t ~ /eEF1A/) print $1
}' "$OUTDIR/evidence/EF1A_blast_top1.tsv" | sort -u > "$OUTDIR/EF1A_keep_by_blast.id"

awk -F'\t' 'BEGIN{IGNORECASE=1}
{
  t=$7
  if (t ~ /EF-Tu/ || t ~ /elongation factor Tu/ || t ~ /TUFM/ || t ~ /mitochondrial elongation factor Tu/ || t ~ /chloroplast elongation factor Tu/ || t ~ /bacterial EF-Tu/) print $1
}' "$OUTDIR/evidence/EF1A_blast_top1.tsv" | sort -u > "$OUTDIR/EF1A_remove_by_blast.id"

# 7) nohit：在rescue_len里但top1里没出现的
cut -f1 "$OUTDIR/evidence/EF1A_blast_top1.tsv" | sort -u > "$OUTDIR/evidence/EF1A_top1.id"
comm -23 <(sort "$OUTDIR/EF1A_rescue_len.id") <(sort "$OUTDIR/evidence/EF1A_top1.id") > "$OUTDIR/EF1A_nohit.id"

# 8) 最终高置信/补充：
#    - high_conf = high_len ∩ keep_by_blast
#    - rescue    = (rescue_len \ high_len) ∩ keep_by_blast
comm -12 <(sort "$OUTDIR/EF1A_high_len.id") <(sort "$OUTDIR/EF1A_keep_by_blast.id") > "$OUTDIR/EF1A_high_conf.id"
comm -23 <(sort "$OUTDIR/EF1A_rescue_len.id") <(sort "$OUTDIR/EF1A_high_len.id") | comm -12 /dev/fd/3 <(sort "$OUTDIR/EF1A_keep_by_blast.id") 3</dev/stdin > "$OUTDIR/EF1A_rescue.id" || true

# 9) 排除表（含原因）
{
  awk '{print $1"\tLEN_OUT_OF_RANGE"}' "$OUTDIR/EF1A_excluded_len.id"
  awk '{print $1"\tBLAST_REMOVE_EF-Tu_like"}' "$OUTDIR/EF1A_remove_by_blast.id"
  awk '{print $1"\tBLAST_NOHIT"}' "$OUTDIR/EF1A_nohit.id"
} | sort -u > "$OUTDIR/EF1A_excluded.tsv"

# 10) 导出FASTA
bash scripts/pfam_pipeline.sh fasta "$PROTEOME" "$OUTDIR/EF1A_high_conf.id" "$OUTDIR/EF1A_high_conf.fa" || true
bash scripts/pfam_pipeline.sh fasta "$PROTEOME" "$OUTDIR/EF1A_rescue.id" "$OUTDIR/EF1A_rescue.fa" || true

# 11) QC
{
  echo "EF1A QC"
  echo "Candidates (PF00009, iE<=1e-5): $(wc -l < "$OUTDIR/EF1A_candidates.id")"
  echo "High_len (420-480): $(wc -l < "$OUTDIR/EF1A_high_len.id")"
  echo "Rescue_len (400-520): $(wc -l < "$OUTDIR/EF1A_rescue_len.id")"
  echo "BLAST keep (eEF1A): $(wc -l < "$OUTDIR/EF1A_keep_by_blast.id")"
  echo "BLAST remove (EF-Tu/TUFM): $(wc -l < "$OUTDIR/EF1A_remove_by_blast.id")"
  echo "NoHit: $(wc -l < "$OUTDIR/EF1A_nohit.id")"
  echo "FINAL high_conf: $(wc -l < "$OUTDIR/EF1A_high_conf.id" 2>/dev/null || echo 0)"
  echo "FINAL rescue: $(wc -l < "$OUTDIR/EF1A_rescue.id" 2>/dev/null || echo 0)"
  echo "Top1 head:"
  head -n 10 "$OUTDIR/evidence/EF1A_blast_top1.tsv"
} > "$OUTDIR/EF1A_qc.txt"

echo "[DONE] outputs in $OUTDIR"
