#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh

clean_file_inplace inputs/proteome.fa
clean_file_inplace results/pfam.tsv

log "Sanitized inputs/proteome.fa and results/pfam.tsv (CR removed, trailing spaces trimmed)."
log "proteome header sample:"
grep -m 1 '^>' inputs/proteome.fa | od -An -tx1 >&2
