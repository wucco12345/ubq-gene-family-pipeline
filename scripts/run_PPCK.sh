#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
source scripts/lib/yaml.sh

F="PPCK"
CFG="config/${F}.yml"
D="results/${F}"
EVD="${D}/evidence"
PROTEOME="inputs/proteome.fa"
PFAM="results/pfam.tsv"
DB="db/uniprot/sprot"
SEEDS="$(yget "$CFG" seeds.fasta)"

mkdir -p "$D" "$EVD" "db/proteome_blast"

test -s "$SEEDS" || { echo "[ERROR] missing seeds fasta: $SEEDS"; exit 1; }
test -s "$PROTEOME" || { echo "[ERROR] missing proteome: $PROTEOME"; exit 1; }
test -s "$PFAM" || { echo "[ERROR] missing pfam table: $PFAM"; exit 1; }

DIR_EVALUE="$(yget "$CFG" directed_blast.evalue 1e-10)"
IDENTITY_MIN="$(yget "$CFG" directed_blast.identity_min 35)"
QCOV_MIN="$(yget "$CFG" directed_blast.query_coverage_min 0.60)"
SCOV_MIN="$(yget "$CFG" directed_blast.subject_coverage_min 0.60)"

CLASSIC_MIN="$(yget "$CFG" length.classic_min 280)"
CLASSIC_MAX="$(yget "$CFG" length.classic_max 320)"

BLAST_TASK="$(yget "$CFG" blast.task blastp)"
BLAST_EVALUE="$(yget "$CFG" blast.evalue 1e-10)"

KEEP_RE="$(yget "$CFG" rules.keep_top1_re)"
RESCUE_RE="$(yget "$CFG" rules.rescue_top1_re)"
REMOVE_RE="$(yget "$CFG" rules.remove_top1_re)"

log "[1/9] Build local BLAST database for proteome"
if [ ! -s "db/proteome_blast/proteome.pin" ]; then
  makeblastdb -in "$PROTEOME" -dbtype prot -out db/proteome_blast/proteome >/dev/null
fi

log "[2/9] Directed BLAST: PPCK seeds vs proteome"
blastp \
  -query "$SEEDS" \
  -db db/proteome_blast/proteome \
  -evalue "$DIR_EVALUE" \
  -max_target_seqs 20000 \
  -num_threads "${THREADS_DEFAULT}" \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
  > "${EVD}/${F}_seed_vs_proteome.tsv"

log "[3/9] Best-hit per subject using merged HSP coverage"
python3 - <<'PY'
from pathlib import Path
from collections import defaultdict

F = "PPCK"
D = Path("results") / F
EVD = D / "evidence"
inp = EVD / f"{F}_seed_vs_proteome.tsv"
out = D / f"{F}_directed_besthit.merged.tsv"

def merge_intervals(intervals):
    if not intervals:
        return []
    arr = sorted((min(a,b), max(a,b)) for a,b in intervals)
    merged = [list(arr[0])]
    for s, e in arr[1:]:
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
    "hsp_count": 0,
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
        g["hsp_count"] += 1

        if (
            g["best_evalue"] is None or
            evalue < g["best_evalue"] or
            (evalue == g["best_evalue"] and bitscore > g["best_bitscore"])
        ):
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
        "seed_len": g["qlen"],
        "protein_len": g["slen"],
        "hsp_count": g["hsp_count"],
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
            (r["bitscore"] == cur["bitscore"] and r["evalue"] == cur["evalue"] and r["query_coverage"] > cur["query_coverage"]) or
            (r["bitscore"] == cur["bitscore"] and r["evalue"] == cur["evalue"] and r["query_coverage"] == cur["query_coverage"] and r["subject_coverage"] > cur["subject_coverage"])
        )
        if better:
            best_by_subject[sid] = r

with out.open("w") as fo:
    fo.write("\t".join([
        "subject_id","best_seed","pident",
        "align_len_merged_query","align_len_merged_subject",
        "query_coverage","subject_coverage",
        "evalue","bitscore","seed_len","protein_len","hsp_count"
    ]) + "\n")
    for sid in sorted(best_by_subject):
        r = best_by_subject[sid]
        fo.write(
            f"{r['subject_id']}\t{r['best_seed']}\t{r['pident']:.3f}\t"
            f"{r['align_len_merged_query']}\t{r['align_len_merged_subject']}\t"
            f"{r['query_coverage']:.3f}\t{r['subject_coverage']:.3f}\t"
            f"{r['evalue']:.3e}\t{r['bitscore']:.1f}\t"
            f"{r['seed_len']}\t{r['protein_len']}\t{r['hsp_count']}\n"
        )
PY

log "[4/9] Build directed candidate pool"
awk -F'\t' -v PID="$IDENTITY_MIN" -v QC="$QCOV_MIN" -v SC="$SCOV_MIN" '
NR==1{next}
($3+0)>=PID && ($6+0)>=QC && ($7+0)>=SC {print $1}
' "${D}/${F}_directed_besthit.merged.tsv" | LC_ALL=C sort -u > "${D}/${F}_candidates.id"

cp "${D}/${F}_candidates.id" "${D}/${F}_candidates.raw.id"
extract_fa_by_id "${D}/${F}_candidates.id" "$PROTEOME" "${D}/${F}_candidates.fa"
log "Directed merged candidates: $(wc -l < "${D}/${F}_candidates.id")"

log "[5/9] Domain status: kinase required, EF-hand forbidden"
python3 - <<'PY'
from pathlib import Path
from collections import defaultdict

F = "PPCK"
D = Path("results") / F
cand = {x.strip() for x in (D / f"{F}_candidates.id").read_text().splitlines() if x.strip()}
pfam = Path("results/pfam.tsv")
out = D / f"{F}_domain_status.tsv"

KINASE = {"PF00069", "PF07714"}
EFHAND = {"PF00036", "PF13499", "PF13202"}

def norm_pf(x):
    return x.split(".")[0]

hits = defaultdict(list)

with pfam.open() as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 13:
            continue
        pid = parts[2].strip()
        if pid not in cand:
            continue
        pfacc = norm_pf(parts[1].strip())
        try:
            ie = float(parts[8])
        except:
            continue
        if ie > 1e-5:
            continue
        hits[pid].append(pfacc)

with out.open("w") as fo:
    fo.write("protein_id\tkinase_domain\tefhand_domain\tpfam_list\n")
    for pid in sorted(cand):
        pfset = sorted(set(hits.get(pid, [])))
        has_kinase = 1 if any(x in KINASE for x in pfset) else 0
        has_efhand = 1 if any(x in EFHAND for x in pfset) else 0
        fo.write(f"{pid}\t{has_kinase}\t{has_efhand}\t{','.join(pfset) if pfset else 'NA'}\n")
PY

awk -F'\t' 'NR>1 && $2==1 && $3==0{print $1}' "${D}/${F}_domain_status.tsv" | LC_ALL=C sort -u > "${D}/${F}_structure_pass.id"
extract_fa_by_id "${D}/${F}_structure_pass.id" "$PROTEOME" "${D}/${F}_structure_pass.fa"
log "Structure pass: $(wc -l < "${D}/${F}_structure_pass.id")"

log "[6/9] Length QC table"
fa_len_tsv "${D}/${F}_structure_pass.fa" "${EVD}/${F}_len.tsv"
len_filter_ids "${EVD}/${F}_len.tsv" "$CLASSIC_MIN" "$CLASSIC_MAX" "${D}/${F}_classic_len.id"

log "[7/9] Swiss-Prot reviewed BLAST for structure-pass candidates"
blast_top5 "${D}/${F}_structure_pass.fa" "$DB" "${EVD}/${F}_vs_sprot.tsv" "$THREADS_DEFAULT" "$BLAST_TASK" "$BLAST_EVALUE"
blast_top1_from_top5 "${EVD}/${F}_vs_sprot.tsv" "${EVD}/${F}_blast_top1.tsv" "${D}/${F}_top1.id"

log "[8/9] nohit = structure_pass - top1.id"
comm -23 <(LC_ALL=C sort -u "${D}/${F}_structure_pass.id") <(LC_ALL=C sort -u "${D}/${F}_top1.id") > "${D}/${F}_nohit.id"

log "[9/9] Classify by top1 annotation"
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
LC_ALL=C sort -u "${D}/${F}_classic_len.id" > "${D}/${F}_classic_len.sorted.id"

comm -12 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_classic_len.sorted.id" > "${D}/${F}_high_conf.id"
comm -23 "${D}/${F}_keep_by_blast.sorted.id" "${D}/${F}_high_conf.id" > "${D}/${F}_keep_not_classic.id"
cat "${D}/${F}_rescue_by_blast.sorted.id" "${D}/${F}_keep_not_classic.id" 2>/dev/null | LC_ALL=C sort -u > "${D}/${F}_rescue.id"

extract_fa_by_id "${D}/${F}_high_conf.id" "$PROTEOME" "${D}/${F}_high_conf.fa" || true
extract_fa_by_id "${D}/${F}_rescue.id" "$PROTEOME" "${D}/${F}_rescue.fa" || true

{
  echo "PPCK QC (run_PPCK directed+structure+membership)"
  echo "Directed merged candidates: $(wc -l < "${D}/${F}_candidates.id")"
  echo "Structure pass: $(wc -l < "${D}/${F}_structure_pass.id")"
  echo "Classic_len (${CLASSIC_MIN}-${CLASSIC_MAX}): $(wc -l < "${D}/${F}_classic_len.id")"
  echo "Keep_by_blast: $(wc -l < "${D}/${F}_keep_by_blast.id")"
  echo "Rescue_by_blast: $(wc -l < "${D}/${F}_rescue_by_blast.id")"
  echo "Remove_by_blast: $(wc -l < "${D}/${F}_remove_by_blast.id")"
  echo "Ambiguous_by_blast: $(wc -l < "${D}/${F}_ambiguous_by_blast.id")"
  echo "NoHit: $(wc -l < "${D}/${F}_nohit.id")"
  echo "High_conf: $(wc -l < "${D}/${F}_high_conf.id")"
  echo "Rescue: $(wc -l < "${D}/${F}_rescue.id")"
} > "${D}/${F}_qc.txt"

echo "[DONE] outputs in ${D}"
