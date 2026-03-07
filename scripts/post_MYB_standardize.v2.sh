#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh

F="MYB"
D="results/${F}"
EVD="${D}/evidence"
PROTEOME="inputs/proteome.fa"

REPEAT2_ID="${D}/${F}_repeat2.id"
HIGH_ID="${D}/${F}_high_conf.id"
RESC_ID="${D}/${F}_rescue.id"
NOHIT_ID="${D}/${F}_nohit.id"
REM_ID="${D}/${F}_remove_by_blast.id"
AMB_ID="${D}/${F}_ambiguous_by_blast.id"
TOP1="${EVD}/${F}_blast_top1.tsv"
LEN_TSV="${EVD}/${F}_len.tsv"
REPEAT_TSV="${D}/${F}_repeat_merged_counts.tsv"

for f in "$PROTEOME" "$REPEAT2_ID" "$HIGH_ID" "$RESC_ID" "$LEN_TSV" "$REPEAT_TSV"; do
  test -f "$f" || { echo "[ERROR] missing file: $f"; exit 1; }
done

test -f "$NOHIT_ID" || : > "$NOHIT_ID"
test -f "$REM_ID" || : > "$REM_ID"
test -f "$AMB_ID" || : > "$AMB_ID"
test -f "$TOP1" || : > "$TOP1"

LC_ALL=C sort -u "$REPEAT2_ID" > "${D}/${F}_repeat2.sorted.id"
LC_ALL=C sort -u "$HIGH_ID" > "${D}/${F}_high_conf.sorted.id"
LC_ALL=C sort -u "$RESC_ID" > "${D}/${F}_rescue.sorted.id"
LC_ALL=C sort -u "$NOHIT_ID" > "${D}/${F}_nohit.sorted.id"

cp "${D}/${F}_high_conf.sorted.id" "${D}/${F}_high_conf.id"
cp "${D}/${F}_rescue.sorted.id" "${D}/${F}_rescue.id"
cp "${D}/${F}_nohit.sorted.id" "${D}/${F}_nohit.id"

cat "${D}/${F}_high_conf.id" "${D}/${F}_rescue.id" | LC_ALL=C sort -u > "${D}/${F}_final_union.id"
comm -23 "${D}/${F}_repeat2.sorted.id" "${D}/${F}_final_union.id" > "${D}/${F}_excluded.id"

: > "${D}/${F}_high_conf.fa"
: > "${D}/${F}_rescue.fa"

if test -s "${D}/${F}_high_conf.id"; then
  extract_fa_by_id "${D}/${F}_high_conf.id" "$PROTEOME" "${D}/${F}_high_conf.fa"
fi

if test -s "${D}/${F}_rescue.id"; then
  extract_fa_by_id "${D}/${F}_rescue.id" "$PROTEOME" "${D}/${F}_rescue.fa"
fi

grep -Ff "${D}/${F}_high_conf.id" "$TOP1" > "${D}/${F}_high_conf.top1.tsv" || true
grep -Ff "${D}/${F}_rescue.id" "$TOP1" > "${D}/${F}_rescue.top1.tsv" || true
grep -Ff "${D}/${F}_excluded.id" "$TOP1" > "${D}/${F}_excluded.top1.tsv" || true

python3 - <<'PY'
from pathlib import Path
import csv

F = "MYB"
D = Path("results") / F
EVD = D / "evidence"

excluded_ids = [x.strip() for x in (D / f"{F}_excluded.id").read_text().splitlines() if x.strip()]

def read_idset(path: Path):
    if not path.exists():
        return set()
    return {x.strip() for x in path.read_text().splitlines() if x.strip()}

nohit = read_idset(D / f"{F}_nohit.id")
rem = read_idset(D / f"{F}_remove_by_blast.id")
amb = read_idset(D / f"{F}_ambiguous_by_blast.id")

top1_map = {}
top1_tsv = EVD / f"{F}_blast_top1.tsv"
if top1_tsv.exists():
    with top1_tsv.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 7:
                top1_map[parts[0]] = parts[6]

len_map = {}
len_tsv = EVD / f"{F}_len.tsv"
if len_tsv.exists():
    with len_tsv.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                len_map[parts[0]] = parts[1]

rep_map = {}
rep_tsv = D / f"{F}_repeat_merged_counts.tsv"
if rep_tsv.exists():
    with rep_tsv.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                rep_map[parts[0]] = {
                    "merged_repeat": parts[1],
                    "min_ali_from": parts[2],
                    "max_ali_to": parts[3],
                }

out = D / f"{F}_excluded.tsv"
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow([
        "ID", "reason", "top1_stitle", "length_aa",
        "merged_repeat", "min_ali_from", "max_ali_to"
    ])
    for pid in excluded_ids:
        top1 = top1_map.get(pid, "NA")
        length = len_map.get(pid, "NA")
        info = rep_map.get(pid, {})
        rep = info.get("merged_repeat", "NA")
        min_from = info.get("min_ali_from", "NA")
        max_to = info.get("max_ali_to", "NA")

        reason = "not_in_final_set"
        if pid in nohit:
            reason = "blast_nohit"
        elif pid in rem:
            reason = "hard_remove_non_myb"
        elif pid in amb:
            reason = "ambiguous_top1"

        w.writerow([pid, reason, top1, length, rep, min_from, max_to])
PY

{
  echo "MYB QC (standardized v2)"
  echo "Candidates (PF00249): $(wc -l < "${D}/${F}_repeat2.sorted.id")"
  echo "High_conf: $(wc -l < "${D}/${F}_high_conf.id")"
  echo "Rescue: $(wc -l < "${D}/${F}_rescue.id")"
  echo "Final total: $(cat "${D}/${F}_high_conf.id" "${D}/${F}_rescue.id" | LC_ALL=C sort -u | wc -l)"
  echo "NoHit: $(wc -l < "${D}/${F}_nohit.id")"
  echo "Excluded: $(wc -l < "${D}/${F}_excluded.id")"
} > "${D}/${F}_qc.txt"

echo "[DONE] standardized outputs in ${D}"
