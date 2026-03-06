#!/usr/bin/env bash
set -euo pipefail
yget() {
  local f="$1" key="$2" def="${3-}"
  python3 - "$f" "$key" "$def" <<'PY'
import sys
from pathlib import Path
try:
    import yaml
except Exception:
    raise SystemExit("Missing python3-yaml. Install: sudo apt-get install -y python3-yaml")

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
