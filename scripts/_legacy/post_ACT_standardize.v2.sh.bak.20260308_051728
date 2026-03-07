#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh

F="ACT"
D="results/${F}"
PROTEOME="inputs/proteome.fa"
TOP1="${D}/evidence/${F}_blast_top1.tsv"
TOP5="${D}/evidence/${F}_vs_sprot.tsv"
CAND_RAW="${D}/${F}.id"
PFAM="results/pfam.tsv"
LEN_TSV="${D}/evidence/${F}_len.tsv"

for f in "$PROTEOME" "$TOP1" "$TOP5" "$CAND_RAW" "$PFAM" "$LEN_TSV"; do
  test -s "$f" || { echo "[ERROR] missing/empty: $f"; exit 1; }
done

# ---------- 0) normalize candidates ----------
LC_ALL=C sort -u "$CAND_RAW" > "${D}/${F}_candidates.sorted.id"

# ---------- 1) nohit = candidates - top1.id ----------
cut -f1 "$TOP1" | LC_ALL=C sort -u > "${D}/${F}_top1.id"
comm -23 "${D}/${F}_candidates.sorted.id" "${D}/${F}_top1.id" > "${D}/${F}_nohit.id"

# ---------- 2) top5 consensus per qid ----------
awk -F'\t' 'BEGIN{IGNORECASE=1; OFS="\t"}
{
  q=$1; t=$7
  if (t ~ /actin/ && t !~ /actin-related|arp/) act[q]++
  if (t ~ /actin-related|arp/) arp[q]++
  n[q]++
}
END{
  for (q in n){
    if (!(q in act)) act[q]=0
    if (!(q in arp)) arp[q]=0
    print q, act[q], arp[q], n[q]
  }
}' "$TOP5" | LC_ALL=C sort -k1,1 > "${D}/${F}_top5_consensus.tsv"

# ---------- 3) join with top1 title ----------
cut -f1,7 "$TOP1" | LC_ALL=C sort -k1,1 > "${D}/${F}_top1_title.tsv"

join -t $'\t' -a 1 -e 'NA' -o '0,1.2,2.2,2.3,2.4' \
  "${D}/${F}_top1_title.tsv" "${D}/${F}_top5_consensus.tsv" \
  > "${D}/${F}_classify.tsv"

# ---------- 4) classify: keep ONLY true ACT ----------
: > "${D}/${F}_final_ACT.id"
: > "${D}/${F}_excluded_final.id"

awk -F'\t' 'BEGIN{IGNORECASE=1}
{
  q=$1; title=$2; act=$3+0; arp=$4+0

  is_act_top1 = (title ~ /actin/ && title !~ /actin-related|arp/)
  is_arp_top1 = (title ~ /actin-related|arp/)

  # keep only true ACT
  if (is_act_top1 && act>=3) {print q >> "'"${D}/${F}_final_ACT.id"'"; next}
  if (!is_act_top1 && !is_arp_top1 && act>=4) {print q >> "'"${D}/${F}_final_ACT.id"'"; next}

  # everything else excluded, including ARP
  print q >> "'"${D}/${F}_excluded_final.id"'"
}' "${D}/${F}_classify.tsv"

LC_ALL=C sort -u "${D}/${F}_final_ACT.id" -o "${D}/${F}_final_ACT.id"
LC_ALL=C sort -u "${D}/${F}_excluded_final.id" -o "${D}/${F}_excluded_final.id"

# ---------- 5) second-round unified outputs ----------
cp "${D}/${F}_final_ACT.id" "${D}/${F}_high_conf.id"
: > "${D}/${F}_rescue.id"

extract_fa_by_id "${D}/${F}_high_conf.id" "$PROTEOME" "${D}/${F}_high_conf.fa"
: > "${D}/${F}_rescue.fa"

# keep ACT-specific output
cp "${D}/${F}_final_ACT.id" "${D}/${F}_final_all.id"
extract_fa_by_id "${D}/${F}_final_all.id" "$PROTEOME" "${D}/${F}_final_all.fa"
extract_fa_by_id "${D}/${F}_excluded_final.id" "$PROTEOME" "${D}/${F}_excluded_final.fa"
extract_fa_by_id "${D}/${F}_final_ACT.id" "$PROTEOME" "${D}/${F}_final_ACT.fa"

# remove old ARP outputs if they exist, to avoid confusion
rm -f "${D}/${F}_final_ARP.id" "${D}/${F}_final_ARP.fa" "${D}/${F}_final_ARP.top1.tsv"

# ---------- 6) evidence snapshots ----------
grep -Ff "${D}/${F}_final_ACT.id" "$TOP1" > "${D}/${F}_final_ACT.top1.tsv" || true
grep -Ff "${D}/${F}_excluded_final.id" "$TOP1" > "${D}/${F}_excluded_final.top1.tsv" || true
grep -Ff "${D}/${F}_high_conf.id" "$TOP1" > "${D}/${F}_high_conf.top1.tsv" || true

# ---------- 7) build excluded.tsv ----------
python3 - <<'PY'
from pathlib import Path
import csv

F = "ACT"
D = Path("results") / F
TOP1 = D / "evidence" / f"{F}_blast_top1.tsv"
LEN = D / "evidence" / f"{F}_len.tsv"
PFAM = Path("results/pfam.tsv")

excluded_ids = [x.strip() for x in (D / f"{F}_excluded_final.id").read_text().splitlines() if x.strip()]
nohit_ids = {x.strip() for x in (D / f"{F}_nohit.id").read_text().splitlines() if x.strip()}

# length
len_map = {}
with LEN.open() as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        pid, L = line.split("\t")[:2]
        len_map[pid] = L

# top1
top1_map = {}
with TOP1.open() as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 7:
            top1_map[parts[0]] = parts[6]

# pfam info
pfam_cnt = {}
pfam_hits = {}
with PFAM.open() as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        pf = parts[1]
        pid = parts[2]
        if pf.startswith("PF00022"):
            pfam_cnt[pid] = pfam_cnt.get(pid, 0) + 1
            pfam_hits.setdefault(pid, []).append(pf)

out = D / f"{F}_excluded.tsv"
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["ID", "reason", "top1_stitle", "length_aa", "pfam_hit_count", "pfam_hits"])
    for pid in excluded_ids:
        top1 = top1_map.get(pid, "NA")
        length = len_map.get(pid, "NA")
        pfc = str(pfam_cnt.get(pid, "NA"))
        pfh = ";".join(pfam_hits.get(pid, [])) if pid in pfam_hits else "NA"

        reason = "excluded_non_act"
        if pid in nohit_ids:
            reason = "blast_nohit"
        elif "actin-related" in top1.lower() or "arp" in top1.lower():
            reason = "excluded_arp"

        w.writerow([pid, reason, top1, length, pfc, pfh])
PY

# ---------- 8) unified QC ----------
{
  echo "ACT QC (standardized v2, ACT-only)"
  echo "Candidates (PF00022): $(wc -l < "${D}/${F}_candidates.sorted.id")"
  echo "High_conf (true ACT): $(wc -l < "${D}/${F}_high_conf.id")"
  echo "Rescue: $(wc -l < "${D}/${F}_rescue.id")"
  echo "NoHit: $(wc -l < "${D}/${F}_nohit.id")"
  echo "Excluded final: $(wc -l < "${D}/${F}_excluded_final.id")"
  echo ""
  echo "[QC] high_conf ∩ excluded_final (must be 0):"
  comm -12 <(LC_ALL=C sort -u "${D}/${F}_high_conf.id") <(LC_ALL=C sort -u "${D}/${F}_excluded_final.id") | wc -l
} > "${D}/${F}_qc.txt"

echo "[DONE] ACT standardized v2 outputs updated (ACT-only):"
ls -lh \
  "${D}/${F}_high_conf.id" "${D}/${F}_high_conf.fa" \
  "${D}/${F}_rescue.id" "${D}/${F}_rescue.fa" \
  "${D}/${F}_nohit.id" "${D}/${F}_excluded.tsv" "${D}/${F}_qc.txt" \
  "${D}/${F}_final_ACT.id" "${D}/${F}_final_ACT.fa"
