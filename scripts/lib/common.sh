#!/usr/bin/env bash
set -euo pipefail

die(){ echo "[ERROR] $*" >&2; exit 1; }
log(){ echo "[INFO] $*" >&2; }

# ---- CRLF hygiene ----
clean_file_inplace() {
  local f="$1"
  test -f "$f" || die "missing file: $f"
  # remove CR and trailing spaces
  perl -pi -e 's/\r//g; s/[ \t]+$//g' "$f"
}

# ---- stable sort for comm ----
sortu() { LC_ALL=C sort -u "$1"; }

# ---- set ops: A - B (both should be sorted unique already) ----
setdiff() { comm -23 "$1" "$2"; }
setinter() { comm -12 "$1" "$2"; }

# ---- check non-empty ----
must_exist_nonempty() { test -s "$1" || die "missing/empty: $1"; }

# ---- FASTA header tokens -> id list ----
fasta_ids() {
  # prints 1st token of each header without ">"
  grep '^>' "$1" | sed 's/^>//' | awk '{print $1}'
}

# ---- Extract FASTA by ID list (IDs are plain, no ">") ----
extract_fa_by_id() {
  local idfile="$1"
  local fafile="$2"
  local outfa="$3"
  must_exist_nonempty "$idfile"
  must_exist_nonempty "$fafile"
  awk -v IDS="$idfile" -f scripts/lib/extract_by_id.awk "$fafile" > "$outfa"
}

# ---- Pfam filter: match PFxxxxxx(.version optional), iEvalue threshold ----
pfam_pick_ids() {
  # args: pfam_tsv pf_regex ie_col protein_col out_id
  local pfam="$1" pf_re="$2" ie_col="$3" prot_col="$4" out="$5"
  must_exist_nonempty "$pfam"
  awk -F'\t' -v re="$pf_re" -v iec="$ie_col" -v pc="$prot_col" '
    { gsub(/\r/,"",$0) }
    ($2 ~ re) {
      ie = $(iec)+0
      if(ie <= 1e-5){
        id=$(pc)
        gsub(/\r/,"",id)
        gsub(/[ \t]+$/,"",id)
        print id
      }
    }
  ' "$pfam" | LC_ALL=C sort -u > "$out"
}

# ---- Length table from FASTA (id \t length) ----
fa_len_tsv() {
  local fa="$1" out="$2"
  must_exist_nonempty "$fa"
  awk '
    /^>/{if(id){print id "\t" len} id=substr($0,2); split(id,a,/[ \t]/); id=a[1]; len=0; next}
    {gsub(/\r/,""); len+=length($0)}
    END{if(id){print id "\t" len}}
  ' "$fa" > "$out"
}

# ---- Filter length from len.tsv ----
len_filter_ids() {
  # args: len.tsv min max out.id
  local lent="$1" min="$2" max="$3" out="$4"
  must_exist_nonempty "$lent"
  awk -F'\t' -v mn="$min" -v mx="$max" '{if($2+0>=mn && $2+0<=mx) print $1}' "$lent" \
    | LC_ALL=C sort -u > "$out"
}

# ---- BLAST helpers ----
blast_top5() {
  # args: query.fa dbprefix out.tsv threads task(optional)
  local q="$1" db="$2" out="$3" th="${4:-8}" task="${5:-blastp}"
  must_exist_nonempty "$q"
  # Expect dbprefix.{psq,pin,phr} or BLAST errors
  blastp -task "$task" -query "$q" -db "$db" -evalue 1e-10 -max_target_seqs 5 -num_threads "$th" \
    -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' > "$out"
}

blast_top1_from_top5() {
  # args: top5.tsv out_top1.tsv out_top1.id
  local top5="$1" out1="$2" outid="$3"
  must_exist_nonempty "$top5"
  awk -F'\t' '!seen[$1]++{print}' "$top5" > "$out1"
  cut -f1 "$out1" | LC_ALL=C sort -u > "$outid"
}

# ----------------------------
# Common helpers (NP-safe)
# ----------------------------

die(){ echo "[ERROR] $*" >&2; exit 1; }
info(){ echo "[INFO] $*" >&2; }

# stable sort unique
sortu(){ LC_ALL=C sort -u "$1"; }

# Extract first token of FASTA header (">ID other stuff" -> "ID")
fasta_ids(){
  grep '^>' "$1" | sed 's/^>//' | awk '{print $1}'
}

# Robust FASTA extraction by ID list (handles CR and trailing spaces)
extract_fasta_by_id(){
  local idfile="$1" fafile="$2" outfa="$3"
  test -s "$idfile" || die "ID list empty: $idfile"
  test -s "$fafile" || die "FASTA missing/empty: $fafile"
  awk -v IDS="$idfile" '
  BEGIN{
    while((getline line < IDS) > 0){
      gsub(/\r/, "", line)
      sub(/[ \t]+$/, "", line)
      if(line!="") want[line]=1
    }
    close(IDS)
    keep=0
  }
  /^>/{
    hdr=$0
    gsub(/\r/, "", hdr)
    sub(/^>/, "", hdr)
    split(hdr,a,/[ \t]/)
    id=a[1]
    keep = (id in want)
  }
  keep{
    gsub(/\r/, "", $0)
    print
  }' "$fafile" > "$outfa"
}

# Pick candidates from pfam.tsv by PF (supports versioned PF like PF00240.29)
# args: pfam_tsv pf_regex ie_col protein_col out_id
pfam_pick_candidates(){
  local pfam="$1" pfregex="$2" iecol="$3" protcol="$4" outid="$5"
  test -s "$pfam" || die "pfam missing/empty: $pfam"
  awk -F'\t' -v PFRE="$pfregex" -v IEC="$iecol" -v PC="$protcol" '
    { gsub(/\r/,"",$0) }
    $2 ~ PFRE && ($(IEC)+0) <= 1e-5 {
      id=$(PC)
      gsub(/\r/,"",id)
      gsub(/[ \t]+$/,"",id)
      print id
    }' "$pfam" | LC_ALL=C sort -u > "$outid"
}

# Run BLASTp to Swiss-Prot (top5) and derive top1 table
# args: query.fa dbprefix out_vs.tsv out_top1.tsv threads task evalue
blast_top5_top1(){
  local qfa="$1" db="$2" outvs="$3" outtop1="$4" threads="${5:-8}" task="${6:-blastp}" evalue="${7:-1e-10}"
  test -s "$qfa" || die "BLAST query missing/empty: $qfa"
  # db expects .pin/.psq/.phr
  test -s "${db}.pin" || die "BLAST DB not found: ${db}.pin"
  blastp -task "$task" -query "$qfa" -db "$db" -evalue "$evalue" -max_target_seqs 5 \
    -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
    -num_threads "$threads" > "$outvs"
  awk '!seen[$1]++{print}' "$outvs" > "$outtop1"
}
