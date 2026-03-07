#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
source scripts/lib/yaml.sh

F="NAC"
CFG="config/NAC.yml"
D="results/${F}"
EVD="${D}/evidence"
PROTEOME="inputs/proteome.fa"
PFAM="results/pfam.tsv"
DB="db/uniprot/sprot"

mkdir -p "$D" "$EVD"

DOMAIN_RE="$(yget "$CFG" pfam.domain_re)"
IE_COL="$(yget "$CFG" pfam.ie_col 9)"
PROT_COL="$(yget "$CFG" pfam.protein_col 3)"
ALI_FROM_COL="$(yget "$CFG" pfam.ali_from_col 12)"
ALI_TO_COL="$(yget "$CFG" pfam.ali_to_col 13)"

HIGH_MIN="$(yget "$CFG" length.high_min 150)"
HIGH_MAX="$(yget "$CFG" length.high_max 800)"
RESCUE_MIN="$(yget "$CFG" length.rescue_min 1)"
RESCUE_MAX="$(yget "$CFG" length.rescue_max 2000)"

BLAST_TASK="$(yget "$CFG" blast.task $BLAST_TASK_DEFAULT)"
BLAST_EVALUE="$(yget "$CFG" blast.evalue $BLAST_EVALUE_DEFAULT)"

KEEP_RE="$(yget "$CFG" rules.keep_top1_re)"
RESCUE_RE="$(yget "$CFG" rules.rescue_top1_re)"
REMOVE_RE="$(yget "$CFG" rules.remove_top1_re)"

log "[1/8] Pfam candidates: PF02365"
awk -F'\t' -v re="$DOMAIN_RE" -v IEC="$IE_COL" -v PC="$PROT_COL" '
{
  gsub(/\r/,"",$0)
  ie=$(IEC)+0
  if(ie<=1e-5 && $2 ~ re){
    id=$(PC)
    gsub(/\r/,"",id)
    sub(/[ \t]+$/,"",id)
    print id
  }
}' "$PFAM" | LC_ALL=C sort -u > "${D}/${F}_candidates.id"
log "Candidates (PF02365): $(wc -l < "${D}/${F}_candidates.id")"

log "[2/8] Extract candidate FASTA"
extract_fa_by_id "${D}/${F}_candidates.id" "$PROTEOME" "${D}/${F}_candidates.fa"
log "Candidate FASTA seqs: $(grep -c '^>' "${D}/${F}_candidates.fa")"

log "[3/8] Domain status table"
python3 - <<'PY'
from pathlib import Path

pfam = Path("results/pfam.tsv")
out = Path("results/NAC/NAC_domain_status.tsv")
re_pat = "PF02365"

hits = {}
with pfam.open() as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 13:
            continue
        hmm = parts[1]
        if re_pat not in hmm:
            continue
        try:
            ie = float(parts[8])
        except:
            continue
        if ie > 1e-5:
            continue
        pid = parts[2].strip()
        ali_from = int(float(parts[11]))
        ali_to = int(float(parts[12]))
        hits.setdefault(pid, []).append((ali_from, ali_to))

lens = {}
fa = Path("results/NAC/NAC_candidates.fa")
seq_id = None
seq = []
with fa.open() as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if seq_id is not None:
                lens[seq_id] = len("".join(seq))
            seq_id = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if seq_id is not None:
        lens[seq_id] = len("".join(seq))

with out.open("w") as w:
    for pid in sorted(hits):
        arr = sorted(hits[pid], key=lambda x: x[0])
        min_from = min(x[0] for x in arr)
        min_to = min(x[1] for x in arr)
        max_to = max(x[1] for x in arr)
        hit_n = len(arr)
        L = lens.get(pid, 0)
        domain_pass = 0
        if min_from <= 150:
            domain_pass = 1
        elif L > 0 and max_to <= (L / 3.0):
            domain_pass = 1
        w.write(f"{pid}\t{L}\t{hit_n}\t{min_from}\t{min_to}\t{max_to}\t{domain_pass}\n")
PY

awk -F'\t' '$7==1{print $1}' "${D}/${F}_domain_status.tsv" | LC_ALL=C sort -u > "${D}/${F}_domain_pass.id"

log "[4/8] Length table + QC bins"
awk -F'\t' -v min="$HIGH_MIN" -v max="$HIGH_MAX" '$2>=min && $2<=max{print $1}' "${D}/${F}_domain_status.tsv" | LC_ALL=C sort -u > "${D}/${F}_high_len.id"
awk -F'\t' -v min="$RESCUE_MIN" -v max="$RESCUE_MAX" '$2>=min && $2<=max{print $1}' "${D}/${F}_domain_status.tsv" | LC_ALL=C sort -u > "${D}/${F}_rescue_len.id"

log "[5/8] BLASTp ALL candidates vs Swiss-Prot (top5 -> top1)"
blast_top5 "${D}/${F}_candidates.fa" "$DB" "${EVD}/${F}_vs_sprot.tsv" "$THREADS_DEFAULT" "$BLAST_TASK"
blast_top1_from_top5 "${EVD}/${F}_vs_sprot.tsv" "${EVD}/${F}_blast_top1.tsv" "${D}/${F}_top1.id"

log "[6/8] nohit = candidates - top1.id"
comm -23 <(LC_ALL=C sort -u "${D}/${F}_candidates.id") <(LC_ALL=C sort -u "${D}/${F}_top1.id") > "${D}/${F}_nohit.id"

log "[7/8] Classify by top1 annotation"
: > "${D}/${F}_keep_by_blast.id"
: > "${D}/${F}_rescue_by_blast.id"
: > "${D}/${F}_remove_by_blast.id"
: > "${D}/${F}_ambiguous_by_blast.id"

awk -F'\t' -v KEEP_RE="$KEEP_RE" -v RESCUE_RE="$RESCUE_RE" -v REMOVE_RE="$REMOVE_RE" 'BEGIN{IGNORECASE=1}
{
  q=$1; t=$7
  if(t ~ REMOVE_RE){print q >> "'"${D}/${F}_remove_by_blast.id"'"; next}
  if(t ~ KEEP_RE){print q >> "'"${D}/${F}_keep_by_blast.id"'"; next}
  if(t ~ RESCUE_RE){print q >> "'"${D}/${F}_rescue_by_blast.id"'"; next}
  print q >> "'"${D}/${F}_ambiguous_by_blast.id"'"
}' "${EVD}/${F}_blast_top1.tsv"

LC_ALL=C sort -u "${D}/${F}_keep_by_blast.id" > "${D}/${F}_keep_by_blast.sorted.id"
LC_ALL=C sort -u "${D}/${F}_rescue_by_blast.id" > "${D}/${F}_rescue_by_blast.sorted.id"
LC_ALL=C sort -u "${D}/${F}_high_len.id" > "${D}/${F}_high_len.sorted.id"
LC_ALL=C sort -u "${D}/${F}_domain_pass.id" > "${D}/${F}_domain_pass.sorted.id"

# high_conf = keep_by_blast ∩ high_len ∩ domain_pass
comm -12 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_len.sorted.id" | comm -12 - "${D}/${F}_domain_pass.sorted.id" > "${D}/${F}_high_conf.id"

# rescue:
# 1) 经典NAC但不满足 high_conf（长度或N端规则不通过）
comm -23 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_conf.id" > "${D}/${F}_rescue_from_keep.id"
# 2) NAC-related atypical
cp "${D}/${F}_rescue_by_blast.sorted.id" "${D}/${F}_rescue_from_label.id"

cat "${D}/${F}_rescue_from_keep.id" "${D}/${F}_rescue_from_label.id" | LC_ALL=C sort -u > "${D}/${F}_rescue.id"

log "[8/8] QC"
{
  echo "NAC QC (run_NAC first-pass)"
  echo "Candidates (PF02365): $(wc -l < "${D}/${F}_candidates.id")"
  echo "Domain_pass: $(wc -l < "${D}/${F}_domain_pass.id")"
  echo "High_len (${HIGH_MIN}-${HIGH_MAX}): $(wc -l < "${D}/${F}_high_len.id")"
  echo "Keep_by_blast: $(wc -l < "${D}/${F}_keep_by_blast.id")"
  echo "Rescue_by_blast: $(wc -l < "${D}/${F}_rescue_by_blast.id")"
  echo "Remove_by_blast: $(wc -l < "${D}/${F}_remove_by_blast.id")"
  echo "Ambiguous_by_blast: $(wc -l < "${D}/${F}_ambiguous_by_blast.id")"
  echo "NoHit: $(wc -l < "${D}/${F}_nohit.id")"
  echo "High_conf: $(wc -l < "${D}/${F}_high_conf.id")"
  echo "Rescue: $(wc -l < "${D}/${F}_rescue.id")"
} > "${D}/${F}_qc.txt"

echo "[DONE] outputs in ${D}"
