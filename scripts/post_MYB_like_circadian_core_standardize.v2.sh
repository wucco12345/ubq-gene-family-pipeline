#!/usr/bin/env bash
set -euo pipefail

F="MYB_like_circadian_core"
D="results/${F}"
PROTEOME="inputs/proteome.fa"
TOP="${D}/${F}_best_hits.tsv"

HIGH_ID="${D}/${F}_high_conf.id"
RESCUE_ID="${D}/${F}_rescue.id"
NOHIT_ID="${D}/${F}_nohit.id"
EXCLUDED_ID="${D}/${F}_excluded.id"

for f in "$PROTEOME" "$HIGH_ID" "$EXCLUDED_ID" "$TOP"; do
  test -f "$f" || { echo "[ERROR] missing file: $f"; exit 1; }
done

# allow empty optional files
test -f "$RESCUE_ID" || : > "$RESCUE_ID"
test -f "$NOHIT_ID" || : > "$NOHIT_ID"

LC_ALL=C sort -u "$HIGH_ID" > "${D}/${F}_high_conf.sorted.id"
LC_ALL=C sort -u "$RESCUE_ID" > "${D}/${F}_rescue.sorted.id"
LC_ALL=C sort -u "$NOHIT_ID" > "${D}/${F}_nohit.sorted.id"
LC_ALL=C sort -u "$EXCLUDED_ID" > "${D}/${F}_excluded.sorted.id"

cat "${D}/${F}_high_conf.sorted.id" "${D}/${F}_rescue.sorted.id" | LC_ALL=C sort -u > "${D}/${F}_final_union.id"

cat > "${D}/_extract_by_id.awk" <<'AWK'
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
}
AWK

awk -v IDS="${D}/${F}_high_conf.sorted.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_high_conf.fa"
awk -v IDS="${D}/${F}_rescue.sorted.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_rescue.fa"
awk -v IDS="${D}/${F}_excluded.sorted.id" -f "${D}/_extract_by_id.awk" "$PROTEOME" > "${D}/${F}_excluded.fa"

python3 - <<'PY'
from pathlib import Path

F = "MYB_like_circadian_core"
D = Path("results") / F
top = D / f"{F}_best_hits.tsv"
excluded_ids = {x.strip() for x in (D / f"{F}_excluded.sorted.id").read_text().splitlines() if x.strip()}

out = D / f"{F}_excluded.tsv"
with top.open() as fh, out.open("w") as fo:
    header = fh.readline().rstrip("\n").split("\t")
    idx = {k:i for i,k in enumerate(header)}
    fo.write("subject_id\treason\tbest_seed\tpident\tquery_coverage\tsubject_coverage\tevalue\tbitscore\tprotein_len\n")
    for line in fh:
        arr = line.rstrip("\n").split("\t")
        sid = arr[idx["subject_id"]]
        if sid in excluded_ids:
            fo.write(
                f"{sid}\t{arr[idx['reason']]}\t{arr[idx['best_seed']]}\t"
                f"{arr[idx['pident']]}\t{arr[idx['query_coverage']]}\t{arr[idx['subject_coverage']]}\t"
                f"{arr[idx['evalue']]}\t{arr[idx['bitscore']]}\t{arr[idx['slen']]}\n"
            )
PY

grep -Ff "${D}/${F}_high_conf.sorted.id" "$TOP" > "${D}/${F}_high_conf.top1.tsv" || true
grep -Ff "${D}/${F}_rescue.sorted.id" "$TOP" > "${D}/${F}_rescue.top1.tsv" || true
grep -Ff "${D}/${F}_excluded.sorted.id" "$TOP" > "${D}/${F}_excluded.top1.tsv" || true

{
  echo "MYB_like_circadian_core QC (standardized v2)"
  echo "repeat=1 pool: $(wc -l < "${D}/${F}_repeat1.id")"
  echo "High_conf: $(wc -l < "${D}/${F}_high_conf.sorted.id")"
  echo "Rescue: $(wc -l < "${D}/${F}_rescue.sorted.id")"
  echo "Final total: $(wc -l < "${D}/${F}_final_union.id")"
  echo "NoHit: $(wc -l < "${D}/${F}_nohit.sorted.id")"
  echo "Excluded: $(wc -l < "${D}/${F}_excluded.sorted.id")"
} > "${D}/${F}_qc.txt"

echo "[DONE] standardized outputs in ${D}"
