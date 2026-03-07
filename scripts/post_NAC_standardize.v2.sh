#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh

F="NAC"
D="results/${F}"
EVD="${D}/evidence"
PROTEOME="inputs/proteome.fa"

CAND_ID="${D}/${F}_candidates.id"
HIGH_ID="${D}/${F}_high_conf.id"
RESC_ID="${D}/${F}_rescue.id"
NOHIT_ID="${D}/${F}_nohit.id"
REM_ID="${D}/${F}_remove_by_blast.id"
AMB_ID="${D}/${F}_ambiguous_by_blast.id"
TOP1="${EVD}/${F}_blast_top1.tsv"
LEN_TSV="${D}/${F}_domain_status.tsv"

for f in "$PROTEOME" "$CAND_ID" "$HIGH_ID" "$RESC_ID" "$LEN_TSV"; do
  test -f "$f" || { echo "[ERROR] missing file: $f"; exit 1; }
done

test -f "$NOHIT_ID" || : > "$NOHIT_ID"
test -f "$REM_ID" || : > "$REM_ID"
test -f "$AMB_ID" || : > "$AMB_ID"
test -f "$TOP1" || : > "$TOP1"

LC_ALL=C sort -u "$CAND_ID" > "${D}/${F}_candidates.sorted.id"
LC_ALL=C sort -u "$HIGH_ID" > "${D}/${F}_high_conf.sorted.id"
LC_ALL=C sort -u "$RESC_ID" > "${D}/${F}_rescue.sorted.id"
LC_ALL=C sort -u "$NOHIT_ID" > "${D}/${F}_nohit.sorted.id"

cp "${D}/${F}_high_conf.sorted.id" "${D}/${F}_high_conf.id"
cp "${D}/${F}_rescue.sorted.id" "${D}/${F}_rescue.id"
cp "${D}/${F}_nohit.sorted.id" "${D}/${F}_nohit.id"

cat "${D}/${F}_high_conf.id" "${D}/${F}_rescue.id" | LC_ALL=C sort -u > "${D}/${F}_final_union.id"
comm -23 "${D}/${F}_candidates.sorted.id" "${D}/${F}_final_union.id" > "${D}/${F}_excluded.id"

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

F = "NAC"
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
resc = read_idset(D / f"{F}_rescue.id")

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

domain_map = {}
domain_tsv = D / f"{F}_domain_status.tsv"
if domain_tsv.exists():
    with domain_tsv.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 6:
                domain_map[parts[0]] = {
                    "length_aa": parts[0] if False else None
                }

len_map = {}
with (D / f"{F}_domain_status.tsv").open() as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 6:
            len_map[parts[0]] = {
                "length_aa": parts[1] if False else "NA"
            }

# 直接从 domain_status.tsv 读取完整列:
full_map = {}
with (D / f"{F}_domain_status.tsv").open() as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 6:
            # 兼容当前格式: ID length_aa pfam_hit_count min_ali_from min_ali_to max_ali_to domain_pass
            if len(parts) >= 7:
                full_map[parts[0]] = {
                    "length_aa": parts[1],
                    "pfam_hit_count": parts[2],
                    "min_ali_from": parts[3],
                    "min_ali_to": parts[4],
                    "max_ali_to": parts[5],
                    "domain_pass": parts[6],
                }

out = D / f"{F}_excluded.tsv"
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow([
        "ID", "reason", "top1_stitle", "length_aa",
        "pfam_hit_count", "min_ali_from", "min_ali_to", "max_ali_to", "domain_pass"
    ])
    for pid in excluded_ids:
        info = full_map.get(pid, {})
        reason = "not_in_final_set"
        if pid in nohit:
            reason = "blast_nohit"
        elif pid in rem:
            reason = "hard_remove_non_nac"
        elif pid in amb:
            reason = "ambiguous_top1"

        w.writerow([
            pid,
            reason,
            top1_map.get(pid, "NA"),
            info.get("length_aa", "NA"),
            info.get("pfam_hit_count", "NA"),
            info.get("min_ali_from", "NA"),
            info.get("min_ali_to", "NA"),
            info.get("max_ali_to", "NA"),
            info.get("domain_pass", "NA"),
        ])
PY

{
  echo "NAC QC (standardized v2)"
  echo "Candidates (PF02365): $(wc -l < "${D}/${F}_candidates.sorted.id")"
  echo "High_conf: $(wc -l < "${D}/${F}_high_conf.id")"
  echo "Rescue: $(wc -l < "${D}/${F}_rescue.id")"
  echo "Final total: $(cat "${D}/${F}_high_conf.id" "${D}/${F}_rescue.id" | LC_ALL=C sort -u | wc -l)"
  echo "NoHit: $(wc -l < "${D}/${F}_nohit.id")"
  echo "Excluded: $(wc -l < "${D}/${F}_excluded.id")"
} > "${D}/${F}_qc.txt"

echo "[DONE] standardized outputs in ${D}"
