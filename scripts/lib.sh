#!/usr/bin/env bash
set -euo pipefail

die(){ echo "[ERROR] $*" >&2; exit 1; }
info(){ echo "[INFO] $*" >&2; }
ok(){ echo "[OK] $*" >&2; }

# Normalize text file to Unix LF (strip CR). In-place.
normalize_tsv_inplace() {
  local f="$1"
  test -f "$f" || die "missing file: $f"
  tr -d '\r' < "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
}

# Clean an ID list: strip CR, trim spaces, drop empty, unique+sort.
clean_id_list() {
  local in="$1" out="$2"
  test -s "$in" || die "missing/empty id list: $in"
  awk '{
    gsub(/\r/,"",$0);
    sub(/[[:space:]]+$/,"",$0);
    sub(/^[[:space:]]+/,"",$0);
    if($0!="") print $0
  }' "$in" | sort -u > "$out"
}

# Extract FASTA records by exact match of header first token (after '>').
extract_fa_by_id() {
  local idfile="$1" fasta="$2" outfa="$3"
  test -s "$idfile" || die "missing/empty idfile: $idfile"
  test -s "$fasta"  || die "missing/empty fasta: $fasta"
  awk -v IDFILE="$idfile" '
    BEGIN{
      while((getline line < IDFILE) > 0){
        gsub(/\r/,"",line); sub(/[[:space:]]+$/,"",line); sub(/^[[:space:]]+/,"",line);
        if(line!="") want[line]=1;
      }
      close(IDFILE);
      keep=0;
    }
    /^>/{
      hdr=substr($0,2);
      split(hdr,a,/[ \t]/);
      keep = (a[1] in want);
    }
    keep{print}
  ' "$fasta" > "$outfa"
}

# Parse Pfam TSV and emit protein IDs (col3) for a Pfam accession like PF00022 or PF00240
# Works if column 2 is like PF00240.29 (with version).
pfam_candidates_by_pf() {
  local pf="$1" ie_thr="$2" pfam_tsv="$3" outid="$4"
  test -s "$pfam_tsv" || die "missing/empty pfam_tsv: $pfam_tsv"
  awk -F'\t' -v PF="$pf" -v THR="$ie_thr" '
    BEGIN{OFS="\t"}
    {
      gsub(/\r/,"",$0);
      pfacc=$2;
      ie=$9+0;
      if(pfacc ~ ("^" PF "(\\.|$)") && ie <= THR){
        id=$3; gsub(/\r/,"",id); sub(/[[:space:]]+$/,"",id);
        print id
      }
    }
  ' "$pfam_tsv" | sort -u > "$outid"
}

# Compute protein length table from FASTA: "<id>\t<len>"
fasta_len_tsv() {
  local fasta="$1" outtsv="$2"
  test -s "$fasta" || die "missing/empty fasta: $fasta"
  awk '
    /^>/{if(id!=""){print id"\t"len}; hdr=substr($0,2); split(hdr,a,/[ \t]/); id=a[1]; len=0; next}
    {gsub(/[ \t\r]/,""); len+=length($0)}
    END{if(id!=""){print id"\t"len}}
  ' "$fasta" > "$outtsv"
}

# BLASTp helper (expects db prefix exists). Output outfmt6 with stitle at col7.
blastp_top5() {
  local query_fa="$1" dbprefix="$2" outtsv="$3" threads="${4:-8}" task="${5:-blastp}"
  test -s "$query_fa" || die "missing/empty query: $query_fa"
  test -f "${dbprefix}.psq" || die "Swiss-Prot BLAST DB not found at prefix: $dbprefix (need .psq/.pin/.phr)"
  blastp -task "$task" \
    -query "$query_fa" -db "$dbprefix" \
    -evalue 1e-10 -max_target_seqs 5 -num_threads "$threads" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
    > "$outtsv"
}

# Derive top1 table from top5 (takes first hit per query).
top1_from_top5() {
  local top5="$1" outtop1="$2"
  test -s "$top5" || die "missing/empty top5: $top5"
  awk '!seen[$1]++{print}' "$top5" > "$outtop1"
}
