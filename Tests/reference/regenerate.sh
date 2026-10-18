#!/usr/bin/env bash
# Regenerates every reference file from Where It Isn't's own modules (Node 22, Python 3 + numpy).
#   Tests/reference/regenerate.sh /path/to/where-it-isnt
# CI runs this against a fresh checkout and fails if any .tsv changes.
set -euo pipefail
WII="${1:?usage: regenerate.sh /path/to/where-it-isnt}"
HERE="$(cd "$(dirname "$0")" && pwd)"
for gen in "$HERE"/gen_*.mjs; do
  name="$(basename "$gen" .mjs)"; name="${name#gen_}"
  node "$gen" "$WII" > "$HERE/$name.tsv"
  echo "$name.tsv  $(grep -vc '^#' "$HERE/$name.tsv") cases"
done
for gen in "$HERE"/gen_*.py; do
  name="$(basename "$gen" .py)"; name="${name#gen_}"
  # A .py next to a .mjs of the same name is that generator's helper (gen_ops.mjs runs gen_ops.py
  # and appends its cases, and likewise for control); the .mjs already wrote the file.
  [ -e "$HERE/gen_$name.mjs" ] && continue
  python3 "$gen" "$WII" > "$HERE/$name.tsv"
  echo "$name.tsv  $(grep -vc '^#' "$HERE/$name.tsv") cases"
done
