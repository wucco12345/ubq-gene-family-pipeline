#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import math

INFILE = Path("results/MYB_like_circadian/At_seed_vs_proteome.tsv")
OUTDIR = Path("results/MYB_like_circadian_core")
OUTDIR.mkdir(parents=True, exist_ok=True)

PAIR_OUT = OUTDIR / "MYB_like_circadian_core_pairwise_fixed.tsv"
BEST_OUT = OUTDIR / "MYB_like_circadian_core_best_hits.tsv"
EXC_OUT = OUTDIR / "MYB_like_circadian_core_excluded.tsv"
HC_OUT = OUTDIR / "MYB_like_circadian_core_high_conf.id"
RESC_OUT = OUTDIR / "MYB_like_circadian_core_rescue.id"
QC_OUT = OUTDIR / "MYB_like_circadian_core_qc.txt"

# 只保留核心 circadian clade seeds
CORE_SEEDS = {
    "sp|P92973|CCA1_ARATH",
    "sp|Q6R0H1|LHY_ARATH",
    "sp|Q6R0G4|RVE4_ARATH",
    "sp|Q8H0W3|RVE6_ARATH",
    "sp|Q8RWU3|RVE8_ARATH",
}

# 阈值（严格版）
MIN_IDENT = 40.0
MIN_QCOV = 0.60

def merge_intervals(intervals):
    if not intervals:
        return []
    xs = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [list(xs[0])]
    for s, e in xs[1:]:
        if s <= merged[-1][1] + 1:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return merged

def merged_span(intervals):
    return sum(e - s + 1 for s, e in merge_intervals(intervals))

# pair_data[(qseqid, sseqid)] = {...}
pair_data = {}

with INFILE.open() as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 14:
            continue

        qseqid = parts[0]
        sseqid = parts[1]

        if qseqid not in CORE_SEEDS:
            continue

        pident = float(parts[2])
        length = int(parts[3])
        qstart = int(parts[6])
        qend = int(parts[7])
        sstart = int(parts[8])
        send = int(parts[9])
        evalue = float(parts[10])
        bitscore = float(parts[11])
        qlen = int(parts[12])
        slen = int(parts[13])

        key = (qseqid, sseqid)
        if key not in pair_data:
            pair_data[key] = {
                "qseqid": qseqid,
                "sseqid": sseqid,
                "qlen": qlen,
                "slen": slen,
                "q_intervals": [],
                "s_intervals": [],
                "best_evalue": evalue,
                "best_bitscore": bitscore,
                "sum_len": 0,
                "sum_ident_len": 0.0,
            }

        d = pair_data[key]
        d["q_intervals"].append((qstart, qend))
        d["s_intervals"].append((sstart, send))
        d["sum_len"] += length
        d["sum_ident_len"] += pident * length
        if evalue < d["best_evalue"]:
            d["best_evalue"] = evalue
        if bitscore > d["best_bitscore"]:
            d["best_bitscore"] = bitscore

# pairwise summary
pair_rows = []
for key, d in pair_data.items():
    qspan = merged_span(d["q_intervals"])
    sspan = merged_span(d["s_intervals"])
    qcov = qspan / d["qlen"] if d["qlen"] > 0 else 0.0
    scov = sspan / d["slen"] if d["slen"] > 0 else 0.0
    wpident = d["sum_ident_len"] / d["sum_len"] if d["sum_len"] > 0 else 0.0

    pair_rows.append({
        "qseqid": d["qseqid"],
        "sseqid": d["sseqid"],
        "pident": wpident,
        "align_len_merged_query": qspan,
        "align_len_merged_subject": sspan,
        "query_coverage": qcov,
        "subject_coverage": scov,
        "evalue": d["best_evalue"],
        "bitscore": d["best_bitscore"],
        "qlen": d["qlen"],
        "slen": d["slen"],
    })

pair_rows.sort(key=lambda x: (x["qseqid"], x["sseqid"]))

with PAIR_OUT.open("w") as out:
    out.write(
        "seed_id\tsubject_id\tpident\talign_len_merged_query\talign_len_merged_subject\t"
        "query_coverage\tsubject_coverage\tevalue\tbitscore\tqlen\tslen\n"
    )
    for r in pair_rows:
        out.write(
            f'{r["qseqid"]}\t{r["sseqid"]}\t'
            f'{r["pident"]:.3f}\t{r["align_len_merged_query"]}\t{r["align_len_merged_subject"]}\t'
            f'{r["query_coverage"]:.3f}\t{r["subject_coverage"]:.3f}\t'
            f'{r["evalue"]:.3e}\t{r["bitscore"]:.1f}\t{r["qlen"]}\t{r["slen"]}\n'
        )

# choose best seed per subject
by_subject = defaultdict(list)
for r in pair_rows:
    by_subject[r["sseqid"]].append(r)

best_rows = []
for sseqid, rows in by_subject.items():
    rows.sort(key=lambda x: (-x["bitscore"], x["evalue"], -x["pident"], -x["query_coverage"]))
    best = rows[0]

    # strict core rule:
    # high_conf: pident >= 40 and qcov >= 0.60
    # rescue: 30<=pident<40 and qcov>=0.60, or pident>=40 and 0.40<=qcov<0.60
    # excluded: everything else
    reason = "excluded"
    if best["pident"] >= MIN_IDENT and best["query_coverage"] >= MIN_QCOV:
        reason = "high_conf"
    elif ((30.0 <= best["pident"] < MIN_IDENT and best["query_coverage"] >= MIN_QCOV) or
          (best["pident"] >= MIN_IDENT and 0.40 <= best["query_coverage"] < MIN_QCOV)):
        reason = "rescue"

    best["reason"] = reason
    best_rows.append(best)

best_rows.sort(key=lambda x: (x["reason"], x["subject_id"] if "subject_id" in x else x["sseqid"]))

with BEST_OUT.open("w") as out:
    out.write(
        "subject_id\tbest_seed\tpident\talign_len_merged_query\talign_len_merged_subject\t"
        "query_coverage\tsubject_coverage\tevalue\tbitscore\tqlen\tslen\treason\n"
    )
    for r in sorted(best_rows, key=lambda x: x["sseqid"]):
        out.write(
            f'{r["sseqid"]}\t{r["qseqid"]}\t{r["pident"]:.3f}\t'
            f'{r["align_len_merged_query"]}\t{r["align_len_merged_subject"]}\t'
            f'{r["query_coverage"]:.3f}\t{r["subject_coverage"]:.3f}\t'
            f'{r["evalue"]:.3e}\t{r["bitscore"]:.1f}\t{r["qlen"]}\t{r["slen"]}\t{r["reason"]}\n'
        )

high_conf_ids = sorted(r["sseqid"] for r in best_rows if r["reason"] == "high_conf")
rescue_ids = sorted(r["sseqid"] for r in best_rows if r["reason"] == "rescue")
excluded_rows = sorted((r for r in best_rows if r["reason"] == "excluded"), key=lambda x: x["sseqid"])

HC_OUT.write_text("".join(x + "\n" for x in high_conf_ids))
RESC_OUT.write_text("".join(x + "\n" for x in rescue_ids))

with EXC_OUT.open("w") as out:
    out.write(
        "subject_id\treason\tbest_seed\tpident\tquery_coverage\tsubject_coverage\tevalue\tbitscore\tprotein_len\n"
    )
    for r in excluded_rows:
        out.write(
            f'{r["sseqid"]}\texcluded\t{r["qseqid"]}\t{r["pident"]:.3f}\t'
            f'{r["query_coverage"]:.3f}\t{r["subject_coverage"]:.3f}\t'
            f'{r["evalue"]:.3e}\t{r["bitscore"]:.1f}\t{r["slen"]}\n'
        )

with QC_OUT.open("w") as out:
    out.write("MYB_like_circadian_core QC (fixed coverage)\n")
    out.write(f"repeat=1 pool searched: 520\n")
    out.write(f"High_conf: {len(high_conf_ids)}\n")
    out.write(f"Rescue: {len(rescue_ids)}\n")
    out.write(f"Final total: {len(high_conf_ids) + len(rescue_ids)}\n")
    out.write(f"Excluded: {len(excluded_rows)}\n")

print("[OK] wrote:", PAIR_OUT)
print("[OK] wrote:", BEST_OUT)
print("[OK] wrote:", EXC_OUT)
print("[OK] wrote:", HC_OUT)
print("[OK] wrote:", RESC_OUT)
print("[OK] wrote:", QC_OUT)
