#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh

F="PPCK"
D="results/${F}"
EVD="${D}/evidence"
PROTEOME="inputs/proteome.fa"

CAND_ID="${D}/${F}_candidates.id"
HIGH_ID="${D}/${F}_high_conf.id"
RESC_ID="${D}/${F}_rescue.id"
NOHIT_ID="${D}/${F}_nohit.id"
REM_ID="${D}/${F}_remove_by_blast.id"
AMB_ID="${D}/${F}_ambiguous_by_blast.id"
DOMAIN_TSV="${D}/${F}_domain_status.tsv"
TOP1="${EVD}/${F}_blast_top1.tsv"
LEN_TSV="${EVD}/${F}_len.tsv"

for f in "$PROTEOME" "$CAND_ID" "$HIGH_ID" "$RESC_ID" "$DOMAIN_TSV"; do
  test -f "$f" || { echo "[ERROR] missing file: $f"; exit 1; }
done

test -f "$NOHIT_ID" || : > "$NOHIT_ID"
test -f "$REM_ID" || : > "$REM_ID"
test -f "$AMB_ID" || : > "$AMB_ID"
test -f "$TOP1" || : > "$TOP1"
test -f "$LEN_TSV" || : > "$LEN_TSV"

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

extract_fa_by_id "${D}/${F}_excluded.id" "$PROTEOME" "${D}/${F}_excluded.fa" || true

grep -Ff "${D}/${F}_high_conf.id" "$TOP1" > "${D}/${F}_high_conf.top1.tsv" || true
grep -Ff "${D}/${F}_rescue.id" "$TOP1" > "${D}/${F}_rescue.top1.tsv" || true
grep -Ff "${D}/${F}_excluded.id" "$TOP1" > "${D}/${F}_excluded.top1.tsv" || true

python3 - <<'PY'
from pathlib import Path
import csv

F = "PPCK"
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
            if len(parts) >= 13:
                top1_map[parts[0]] = parts[12]

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

domain_map = {}
domain_tsv = D / f"{F}_domain_status.tsv"
if domain_tsv.exists():
    with domain_tsv.open() as fh:
        next(fh, None)
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                domain_map[parts[0]] = {
                    "kinase_domain": parts[1],
                    "efhand_domain": parts[2],
                    "pfam_list": parts[3],
                }

out = D / f"{F}_excluded.tsv"
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow([
        "ID", "reason", "top1_stitle", "length_aa",
        "kinase_domain", "efhand_domain", "pfam_list"
    ])
    for pid in excluded_ids:
        title = top1_map.get(pid, "NA")
        title_u = title.upper()
        reason = "not_in_final_set"
        if pid in nohit:
            reason = "blast_nohit"
        elif "AURORA" in title_u:
            reason = "hard_remove_aurora_kinase"
        elif "CBL-INTERACTING PROTEIN KINASE" in title_u or "CIPK" in title_u:
            reason = "hard_remove_cipk"
        elif "CALCIUM-DEPENDENT PROTEIN KINASE" in title_u or "CDPK" in title_u:
            reason = "hard_remove_cdpk"
        elif "MAPK" in title_u or "MITOGEN-ACTIVATED PROTEIN KINASE" in title_u:
            reason = "hard_remove_mapk"
        elif pid in rem:
            reason = "hard_remove_non_ppck"
        elif "PEPKR" in title_u:
            reason = "ppck_related_noncore"
        elif pid in amb:
            reason = "ambiguous_top1"

        d = domain_map.get(pid, {})
        w.writerow([
            pid,
            reason,
            title,
            len_map.get(pid, "NA"),
            d.get("kinase_domain", "NA"),
            d.get("efhand_domain", "NA"),
            d.get("pfam_list", "NA"),
        ])
PY

{
  echo "PPCK QC (standardized v2)"
  echo "Candidates: $(wc -l < "${D}/${F}_candidates.sorted.id")"
  echo "High_conf: $(wc -l < "${D}/${F}_high_conf.id")"
  echo "Rescue: $(wc -l < "${D}/${F}_rescue.id")"
  echo "Final total: $(cat "${D}/${F}_high_conf.id" "${D}/${F}_rescue.id" | LC_ALL=C sort -u | wc -l)"
  echo "NoHit: $(wc -l < "${D}/${F}_nohit.id")"
  echo "Excluded: $(wc -l < "${D}/${F}_excluded.id")"
  echo
  echo "[QC] high_conf ∩ rescue (must be 0):"
  comm -12 "${D}/${F}_high_conf.sorted.id" "${D}/${F}_rescue.sorted.id" | wc -l
  echo
  echo "[QC] final_union ∩ excluded (must be 0):"
  comm -12 "${D}/${F}_final_union.id" "${D}/${F}_excluded.id" | wc -l
} > "${D}/${F}_qc.txt"

echo "[DONE] standardized outputs in ${D}"
