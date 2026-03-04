#!/usr/bin/env bash
set -euo pipefail

# ===== User inputs =====
PROTEOME=${1:-""}          # e.g. data/DB.protein.fa
PFAM_HMM=${2:-""}          # e.g. hmm/Pfam-A.hmm
SPROT_DB=${3:-""}          # e.g. db/uniprot/sprot (makeblastdb output prefix)

if [[ -z "$PROTEOME" || -z "$PFAM_HMM" || -z "$SPROT_DB" ]]; then
  echo "Usage: bash run_ubq.sh <PROTEOME.fa> <Pfam-A.hmm> <sprot_db_prefix>"
  echo "Example: bash run_ubq.sh data/DB.protein.fa hmm/Pfam-A.hmm db/uniprot/sprot"
  exit 1
fi

# ===== Outputs =====
OUTDIR="out_ubq"
mkdir -p "$OUTDIR"/{logs,results,tmp}

echo "[1/6] hmmscan PFAM..."
hmmscan --cpu 10 --domtblout "$OUTDIR/results/pfam.domtblout" "$PFAM_HMM" "$PROTEOME" > "$OUTDIR/logs/pfam.hmmscan.out"

echo "[2/6] parse domtblout -> tsv..."
bash scripts/pfam_pipeline.sh parse "$OUTDIR/results/pfam.domtblout" "$OUTDIR/results/pfam.tsv"

echo "[3/6] list PF00240 candidates..."
bash scripts/pfam_pipeline.sh list "$OUTDIR/results/pfam.tsv" PF00240 1e-5 "$OUTDIR/UBQ"

echo "[4/6] repeat counting (merged intervals)..."
bash scripts/pfam_pipeline.sh repeat_merged "$OUTDIR/UBQ.all_hits.tsv" PF00240 "$OUTDIR/UBQ"

echo "[5/6] extract single fasta..."
bash scripts/pfam_pipeline.sh fasta "$PROTEOME" "$OUTDIR/UBQ.single_merged.id" "$OUTDIR/UBQ.single_merged.fa"

echo "[6/6] BLAST single vs Swiss-Prot..."
blastp \
  -query "$OUTDIR/UBQ.single_merged.fa" \
  -db "$SPROT_DB" \
  -evalue 1e-10 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
  -max_target_seqs 5 \
  -num_threads 10 \
  -out "$OUTDIR/UBQ_single_vs_sprot.tsv"

# top1
awk '!seen[$1]++{print}' "$OUTDIR/UBQ_single_vs_sprot.tsv" > "$OUTDIR/UBQ_single_top1.tsv"
cut -f1 "$OUTDIR/UBQ_single_top1.tsv" | sort -u > "$OUTDIR/UBQ_single_top1.id"

# nohit = single - top1
comm -23 <(sort "$OUTDIR/UBQ.single_merged.id") <(sort "$OUTDIR/UBQ_single_top1.id") > "$OUTDIR/UBQ_nohit.id"

# keep_strict: ubiquitin-ribosomal fusion only
awk -F'\t' 'BEGIN{IGNORECASE=1}{t=$7;if(t ~ /ubiquitin[- ]ribosomal/ || (t ~ /ribosomal/ && t ~ /ubiquitin/ && t ~ /fusion/)) print $1}' \
  "$OUTDIR/UBQ_single_top1.tsv" | sort -u > "$OUTDIR/UBQ_keep_strict.id"

# final strict = poly + keep_strict
cat "$OUTDIR/UBQ.poly_merged.id" "$OUTDIR/UBQ_keep_strict.id" | sort -u > "$OUTDIR/UBQ_final_strict.id"
bash scripts/pfam_pipeline.sh fasta "$PROTEOME" "$OUTDIR/UBQ_final_strict.id" "$OUTDIR/UBQ_final_strict.fa"

echo "[DONE] UBQ final:"
wc -l "$OUTDIR/UBQ.poly_merged.id" "$OUTDIR/UBQ_keep_strict.id" "$OUTDIR/UBQ_final_strict.id" "$OUTDIR/UBQ_nohit.id"
echo "Final FASTA: $OUTDIR/UBQ_final_strict.fa"
