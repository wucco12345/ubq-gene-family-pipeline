#!/usr/bin/env bash
set -euo pipefail

# =========================
# Pfam domtblout utilities
# =========================
# Designed for HMMER hmmscan --domtblout output parsed into pfam.tsv.
#
# pfam.tsv columns (tab-separated):
# 1  pfam_name
# 2  pfam_acc          (e.g. PF00240.29)   <-- filter by prefix PF00240
# 3  query_id          (your protein ID)
# 4  query_len
# 5  full_evalue
# 6  full_score
# 7  dom_num
# 8  dom_of
# 9  i_evalue          (domain i-Evalue)  <-- recommended for domain presence
# 10 i_score
# 11 ali_from          (alignment start on query protein)
# 12 ali_to            (alignment end on query protein)
# 13 env_from
# 14 env_to
# 15 hmm_from
# 16 hmm_to

cmd="${1:-}"
if [[ -z "${cmd}" ]]; then
  echo "Usage:"
  echo "  $0 parse         results/pfam.domtblout  results/pfam.tsv"
  echo "  $0 list          results/pfam.tsv PF00240 1e-5 out/UBQ"
  echo "  $0 repeat        out/UBQ.all_hits.tsv PF00240 out/UBQ"
  echo "  $0 repeat_merged out/UBQ.all_hits.tsv PF00240 out/UBQ"
  echo "  $0 fasta         data/DB.protein.fa out/UBQ.single.id out/UBQ.single.fa"
  exit 1
fi

parse_domtblout() {
  local domtblout="$1"
  local out_tsv="$2"
  mkdir -p "$(dirname "$out_tsv")"

  awk 'BEGIN{OFS="\t"}
    /^[#]/ {next}
    NF>0 {
      # domtblout columns used:
      # 1 target name
      # 2 target accession
      # 4 query name
      # 6 qlen
      # 7 full E-value
      # 8 full score
      # 10 #dom
      # 11 of
      # 13 i-Evalue
      # 14 dom score
      # 16 hmm from
      # 17 hmm to
      # 18 ali from
      # 19 ali to
      # 20 env from
      # 21 env to
      print $1,$2,$4,$6,$7,$8,$10,$11,$13,$14,$18,$19,$20,$21,$16,$17
    }' "$domtblout" > "$out_tsv"

  echo "[OK] Parsed domtblout -> $out_tsv"
  echo "[INFO] Rows: $(wc -l < "$out_tsv")"
}

list_by_pfam() {
  local pfam_tsv="$1"
  local pfam_id="$2"
  local emax="$3"
  local prefix="$4"

  mkdir -p "$(dirname "$prefix")"

  # Filter by Pfam accession prefix (col2) and i-Evalue cutoff (col9)
  awk -v PF="$pfam_id" -v EMAX="$emax" 'BEGIN{FS=OFS="\t"}
    $2 ~ ("^"PF) && ($9+0) <= (EMAX+0) {print}
  ' "$pfam_tsv" > "${prefix}.all_hits.tsv"

  awk 'BEGIN{FS="\t"} {print $3}' "${prefix}.all_hits.tsv" | sort -u > "${prefix}.id"

  echo "[OK] PF=$pfam_id (col2) iE<=${emax} (col9) -> ${prefix}.id"
  echo "[INFO] Unique proteins: $(wc -l < "${prefix}.id")"
  echo "[INFO] Hits rows: $(wc -l < "${prefix}.all_hits.tsv")"
}

repeat_count() {
  local pfam_tsv="$1"
  local pfam_id="$2"
  local prefix="$3"

  mkdir -p "$(dirname "$prefix")"

  # Count raw hit rows per query (PF in col2, query in col3)
  awk -v PF="$pfam_id" 'BEGIN{FS=OFS="\t"}
    $2 ~ ("^"PF) {cnt[$3]++}
    END{for (q in cnt) print cnt[q], q}
  ' "$pfam_tsv" | sort -nr > "${prefix}.repeat.tsv"

  awk '$1>=2{print $2}' "${prefix}.repeat.tsv" | sort -u > "${prefix}.poly.id"
  awk '$1==1{print $2}' "${prefix}.repeat.tsv" | sort -u > "${prefix}.single.id"

  echo "[OK] Repeat (raw hits) for PF=$pfam_id (col2)"
  echo "[INFO] poly (>=2): $(wc -l < "${prefix}.poly.id")"
  echo "[INFO] single (=1): $(wc -l < "${prefix}.single.id")"
  echo "[INFO] full table: ${prefix}.repeat.tsv"
}

repeat_count_merged() {
  local pfam_tsv="$1"
  local pfam_id="$2"
  local prefix="$3"

  mkdir -p "$(dirname "$prefix")"

  if ! command -v python3 >/dev/null 2>&1; then
    echo "[ERROR] python3 not found. Install:"
    echo "  sudo apt-get update && sudo apt-get install -y python3"
    exit 2
  fi

  # We merge intervals on query protein using ali_from (col11) and ali_to (col12).
  # This avoids overcounting when a single repeat is split into multiple overlapping hits.
  python3 - "$pfam_tsv" "$pfam_id" "$prefix" <<'PY'
import sys, re

tsv, pf, prefix = sys.argv[1], sys.argv[2], sys.argv[3]
pf_re = re.compile(r"^" + re.escape(pf) + r"(\.|$)")

def merge_intervals(iv):
  if not iv:
    return []
  iv.sort()
  merged = [list(iv[0])]
  for s,e in iv[1:]:
    ls,le = merged[-1]
    if s <= le + 1:   # overlap or touch
      merged[-1][1] = max(le, e)
    else:
      merged.append([s,e])
  return merged

by_q = {}
with open(tsv, "r", encoding="utf-8") as f:
  for line in f:
    line = line.rstrip("\n")
    if not line:
      continue
    cols = line.split("\t")
    if len(cols) < 12:
      continue
    pf_acc = cols[1]
    if not pf_re.match(pf_acc):
      continue
    qid = cols[2]
    try:
      a_from = int(cols[10])
      a_to   = int(cols[11])
    except ValueError:
      continue
    s,e = (a_from, a_to) if a_from <= a_to else (a_to, a_from)
    by_q.setdefault(qid, []).append((s,e))

repeat_path = prefix + ".repeat_merged.tsv"
poly_path   = prefix + ".poly_merged.id"
single_path = prefix + ".single_merged.id"

with open(repeat_path, "w", encoding="utf-8") as out:
  for qid, ivs in by_q.items():
    m = merge_intervals(ivs)
    out.write(f"{len(m)}\t{qid}\t" + ",".join(f"{s}-{e}" for s,e in m) + "\n")

with open(poly_path, "w", encoding="utf-8") as fpoly, open(single_path, "w", encoding="utf-8") as fsingle:
  for line in open(repeat_path, "r", encoding="utf-8"):
    c_str, qid, *_ = line.rstrip("\n").split("\t")
    c = int(c_str)
    if c >= 2:
      fpoly.write(qid + "\n")
    elif c == 1:
      fsingle.write(qid + "\n")

print("[OK] Repeat (merged intervals) written:")
print("  ", repeat_path)
print("  ", poly_path)
print("  ", single_path)
PY

  echo "[INFO] poly_merged (>=2): $(wc -l < "${prefix}.poly_merged.id")"
  echo "[INFO] single_merged (=1): $(wc -l < "${prefix}.single_merged.id")"
}

extract_fasta() {
  local fasta="$1"
  local idlist="$2"
  local outfa="$3"

  if ! command -v seqkit >/dev/null 2>&1; then
    echo "[ERROR] seqkit not found. Install:"
    echo "  sudo apt-get update && sudo apt-get install -y seqkit"
    exit 2
  fi

  seqkit grep -f "$idlist" "$fasta" > "$outfa"
  echo "[OK] Extracted FASTA -> $outfa"
  echo "[INFO] Seqs: $(grep -c "^>" "$outfa" || true)"
}

case "$cmd" in
  parse)
    parse_domtblout "${2:?domtblout required}" "${3:?out_tsv required}"
    ;;
  list)
    list_by_pfam "${2:?pfam_tsv required}" "${3:?pfam_id required}" "${4:?evalue_cutoff required}" "${5:?prefix required}"
    ;;
  repeat)
    repeat_count "${2:?pfam_tsv required}" "${3:?pfam_id required}" "${4:?prefix required}"
    ;;
  repeat_merged)
    repeat_count_merged "${2:?pfam_tsv required}" "${3:?pfam_id required}" "${4:?prefix required}"
    ;;
  fasta)
    extract_fasta "${2:?fasta required}" "${3:?idlist required}" "${4:?outfa required}"
    ;;
  *)
    echo "Unknown cmd: $cmd"
    exit 1
    ;;
esac
