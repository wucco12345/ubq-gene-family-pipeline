#!/usr/bin/env bash
set -euo pipefail

# yget <yaml_file> <dot.key> [default]
# Example: yget config/global.yml blast.evalue 1e-10
yget() {
  local f="$1" key="$2" def="${3-}"
  python3 - "$f" "$key" "$def" <<'PY'
import sys, yaml
from pathlib import Path

f, key, default = sys.argv[1], sys.argv[2], sys.argv[3]
data = yaml.safe_load(Path(f).read_text(encoding="utf-8")) or {}

cur = data
ok = True
for part in key.split("."):
    if isinstance(cur, dict) and part in cur:
        cur = cur[part]
    else:
        ok = False
        break

if not ok or cur is None:
    print(default)
else:
    print(cur)
PY
}
