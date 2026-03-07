#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
source scripts/lib/yaml.sh

F="MYB_like_circadian_core"
CFG="config/${F}.yml"
D="results/${F}"
EVD="${D}/evidence"
PROTEOME="inputs/proteome.fa"
PFAM="results/pfam.tsv"
SEEDS="$(yget "$CFG" seeds.fasta)"

mkdir -p "$D" "$EVD"

PF_RE="$(yget "$CFG" pfam.pf_re)"
IE_COL="$(yget "$CFG" pfam.ie_col 9)"
PROT_COL="$(yget "$CFG" pfam.protein_col 3)"
ALI_FROM_COL="$(yget "$CFG" pfam.ali_from_col 12)"
ALI_TO_COL="$(yget "$CFG" pfam.ali_to_col 13)"

HIGH_MIN="$(yget "$CFG" length.high_min 200)"
HIGH_MAX="$(yget "$CFG" length.high_max 450)"
RESCUE_MIN="$(yget "$CFG" length.rescue_min 150)"
RESCUE_MAX="$(yget "$CFG" length.rescue_max 650)"

BLAST_EVALUE="$(yget "$CFG" blast.evalue 1e-10)"
IDENTITY_MIN="$(yget "$CFG" blast.identity_min 40)"
COVERAGE_MIN="$(yget "$CFG" blast.coverage_min 0.60)"
CORE_RE="$(yget "$CFG" seeds.core_re)"

test -s "$SEEDS" || { echo "[ERROR] missing seeds fasta: $SEEDS"; exit 1; }

log "[1/9] Pfam candidates: PF00249 (iE<=1e-5)"
pfam_pick_ids "$PFAM" "$PF_RE" "$IE_COL" "$PROT_COL" "${D}/${F}_candidates.id"
cp "${D}/${F}_candidates.id" "${D}/${F}_candidates.raw.id"
log "Candidates: $(wc -l < "${D}/${F}_candidates.id")"

log "[2/9] Extract candidate FASTA"
extract_fa_by_id "${D}/${F}_candidates.id" "$PROTEOME" "${D}/${F}_candidates.fa"
log "Candidate FASTA seqs: $(grep -c '^>' "${D}/${F}_candidates.fa")"

log "[3/9] Save PFam hit intervals"
awk -F'\t' -v re="$PF_RE" -v IEC="$IE_COL" -v PC="$PROT_COL" -v AF="$ALI_FROM_COL" -v AT="$ALI_TO_COL" '
BEGIN{OFS="\t"}
{
  gsub(/\r/,"",$0)
}
($2 ~ re){
  ie=$(IEC)+0
  if(ie<=1e-5){
    id=$(PC)
    gsub(/\r/,"",id)
    sub(/[ \t]+$/,"",id)
    print id, $(AF)+0, $(AT)+0, ie, $2
  }
}' "$PFAM" | LC_ALL=C sort -k1,1 -k2,2n > "${EVD}/${F}_pfam_hits.tsv"

log "[4/9] Merge intervals and count repeats"
awk '
BEGIN{OFS="\t"}
{
  id=$1; s=$2+0; e=$3+0
  if(id!=cur && cur!=""){
    print cur, seg, min_from, max_to
  }
  if(id!=cur){
    cur=id
    seg=1
    ms=s; me=e
    min_from=s
    max_to=e
  } else {
    if(s <= me+1){
      if(e>me) me=e
    } else {
      seg++
      ms=s; me=e
    }
    if(s<min_from) min_from=s
    if(e>max_to) max_to=e
  }
}
END{
  if(cur!="") print cur, seg, min_from, max_to
}' "${EVD}/${F}_pfam_hits.tsv" > "${D}/${F}_repeat_merged_counts.tsv"

awk '$2==1{print $1}' "${D}/${F}_repeat_merged_counts.tsv" | LC_ALL=C sort -u > "${D}/${F}_repeat1.id"
awk '$2==2{print $1}' "${D}/${F}_repeat_merged_counts.tsv" | LC_ALL=C sort -u > "${D}/${F}_repeat2.id"
awk '$2>=3{print $1}' "${D}/${F}_repeat_merged_counts.tsv" | LC_ALL=C sort -u > "${D}/${F}_repeat3plus.id"

log "repeat=1: $(wc -l < "${D}/${F}_repeat1.id")"
log "repeat=2: $(wc -l < "${D}/${F}_repeat2.id")"
log "repeat>=3: $(wc -l < "${D}/${F}_repeat3plus.id")"

log "[5/9] Keep only repeat=1 for core circadian MYB-like"
extract_fa_by_id "${D}/${F}_repeat1.id" "$PROTEOME" "${D}/${F}_repeat1.fa"
fa_len_tsv "${D}/${F}_repeat1.fa" "${EVD}/${F}_len.tsv"

len_filter_ids "${EVD}/${F}_len.tsv" "$HIGH_MIN" "$HIGH_MAX" "${D}/${F}_high_len.id"
len_filter_ids "${EVD}/${F}_len.tsv" "$RESCUE_MIN" "$RESCUE_MAX" "${D}/${F}_rescue_len.id"

log "[6/9] Directed BLAST: Arabidopsis core seeds vs repeat=1 proteome subset"
blastp \
  -query "$SEEDS" \
  -subject "${D}/${F}_repeat1.fa" \
  -evalue "$BLAST_EVALUE" \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
  > "${EVD}/${F}_seed_vs_repeat1.tsv"

log "[7/9] Best-hit per subject (merged HSP coverage)"
python3 - <<'PY'
from pathlib import Path
from collections import defaultdict
import re

F = "MYB_like_circadian_core"
D = Path("results") / F
EVD = D / "evidence"

inp = EVD / f"{F}_seed_vs_repeat1.tsv"
out = D / f"{F}_best_hits.tsv"

IDENTITY_MIN = 40.0
COVERAGE_MIN = 0.60
CORE_RE = re.compile(r"(CCA1|LHY|RVE4|RVE6|RVE8)", re.I)

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return merged

def merged_len(intervals):
    return sum(e - s + 1 for s, e in merge_intervals(intervals))

groups = defaultdict(lambda: {
    "q_intervals": [],
    "s_intervals": [],
    "best_evalue": None,
    "best_bitscore": None,
    "best_pident": None,
    "qlen": None,
    "slen": None,
})

with inp.open() as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen = line.split("\t")
        pident = float(pident)
        qstart, qend = sorted((int(qstart), int(qend)))
        sstart, send = sorted((int(sstart), int(send)))
        evalue = float(evalue)
        bitscore = float(bitscore)
        qlen = int(qlen)
        slen = int(slen)

        g = groups[(sseqid, qseqid)]
        g["q_intervals"].append((qstart, qend))
        g["s_intervals"].append((sstart, send))
        g["qlen"] = qlen
        g["slen"] = slen

        if g["best_evalue"] is None or evalue < g["best_evalue"] or (evalue == g["best_evalue"] and bitscore > g["best_bitscore"]):
            g["best_evalue"] = evalue
            g["best_bitscore"] = bitscore
            g["best_pident"] = pident

rows = []
for (subject_id, seed_id), g in groups.items():
    qcov = merged_len(g["q_intervals"]) / g["qlen"] if g["qlen"] else 0.0
    scov = merged_len(g["s_intervals"]) / g["slen"] if g["slen"] else 0.0
    rows.append({
        "subject_id": subject_id,
        "best_seed": seed_id,
        "pident": g["best_pident"],
        "align_len_merged_query": merged_len(g["q_intervals"]),
        "align_len_merged_subject": merged_len(g["s_intervals"]),
        "query_coverage": qcov,
        "subject_coverage": scov,
        "evalue": g["best_evalue"],
        "bitscore": g["best_bitscore"],
        "qlen": g["qlen"],
        "slen": g["slen"],
    })

best_by_subject = {}
for r in rows:
    sid = r["subject_id"]
    cur = best_by_subject.get(sid)
    if cur is None:
        best_by_subject[sid] = r
    else:
        better = (
            (r["bitscore"] > cur["bitscore"]) or
            (r["bitscore"] == cur["bitscore"] and r["evalue"] < cur["evalue"]) or
            (r["bitscore"] == cur["bitscore"] and r["evalue"] == cur["evalue"] and r["query_coverage"] > cur["query_coverage"])
        )
        if better:
            best_by_subject[sid] = r

with out.open("w") as fo:
    fo.write("\t".join([
        "subject_id","best_seed","pident",
        "align_len_merged_query","align_len_merged_subject",
        "query_coverage","subject_coverage",
        "evalue","bitscore","qlen","slen","reason"
    ]) + "\n")

    for sid in sorted(best_by_subject):
        r = best_by_subject[sid]
        ident = r["pident"]
        qcov = r["query_coverage"]
        scov = r["subject_coverage"]
        slen = r["slen"]
        seed = r["best_seed"]
        core_seed = bool(CORE_RE.search(seed))

        if core_seed and ident >= IDENTITY_MIN and qcov >= COVERAGE_MIN and scov >= COVERAGE_MIN and 200 <= slen <= 450:
            reason = "high_conf"
        elif core_seed and ident >= IDENTITY_MIN and qcov >= COVERAGE_MIN and scov >= COVERAGE_MIN and 150 <= slen <= 650:
            reason = "rescue"
        else:
            reason = "excluded"

        fo.write(
            f"{r['subject_id']}\t{r['best_seed']}\t{r['pident']:.3f}\t"
            f"{r['align_len_merged_query']}\t{r['align_len_merged_subject']}\t"
            f"{r['query_coverage']:.3f}\t{r['subject_coverage']:.3f}\t"
            f"{r['evalue']:.3e}\t{r['bitscore']:.1f}\t{r['qlen']}\t{r['slen']}\t{reason}\n"
        )
PY

log "[8/9] Build high_conf / rescue / excluded"
awk -F'\t' 'NR>1 && $12=="high_conf"{print $1}' "${D}/${F}_best_hits.tsv" | LC_ALL=C sort -u > "${D}/${F}_high_conf.id"
awk -F'\t' 'NR>1 && $12=="rescue"{print $1}' "${D}/${F}_best_hits.tsv" | LC_ALL=C sort -u > "${D}/${F}_rescue.id"
awk -F'\t' 'NR>1 && $12=="excluded"{print $1}' "${D}/${F}_best_hits.tsv" | LC_ALL=C sort -u > "${D}/${F}_excluded.id"
: > "${D}/${F}_nohit.id"

extract_fa_by_id "${D}/${F}_high_conf.id" "$PROTEOME" "${D}/${F}_high_conf.fa" || true
extract_fa_by_id "${D}/${F}_rescue.id" "$PROTEOME" "${D}/${F}_rescue.fa" || true

python3 - <<'PY'
from pathlib import Path

F = "MYB_like_circadian_core"
D = Path("results") / F
best = D / f"{F}_best_hits.tsv"
excluded = D / f"{F}_excluded.tsv"

with best.open() as fh, excluded.open("w") as out:
    header = fh.readline().rstrip("\n").split("\t")
    idx = {k:i for i,k in enumerate(header)}
    out.write("subject_id\treason\tbest_seed\tpident\tquery_coverage\tsubject_coverage\tevalue\tbitscore\tprotein_len\n")
    for line in fh:
        arr = line.rstrip("\n").split("\t")
        if arr[idx["reason"]] == "excluded":
            out.write(
                f"{arr[idx['subject_id']]}\texcluded\t{arr[idx['best_seed']]}\t"
                f"{arr[idx['pident']]}\t{arr[idx['query_coverage']]}\t{arr[idx['subject_coverage']]}\t"
                f"{arr[idx['evalue']]}\t{arr[idx['bitscore']]}\t{arr[idx['slen']]}\n"
            )
PY

log "[9/9] QC"
{
  echo "MYB_like_circadian_core QC (standardized v2)"
  echo "repeat=1 pool: $(wc -l < "${D}/${F}_repeat1.id")"
  echo "High_conf: $(wc -l < "${D}/${F}_high_conf.id")"
  echo "Rescue: $(wc -l < "${D}/${F}_rescue.id")"
  echo "Final total: $(cat "${D}/${F}_high_conf.id" "${D}/${F}_rescue.id" | LC_ALL=C sort -u | wc -l)"
  echo "Excluded: $(wc -l < "${D}/${F}_excluded.id")"
} > "${D}/${F}_qc.txt"

echo "[DONE] outputs in ${D}"
