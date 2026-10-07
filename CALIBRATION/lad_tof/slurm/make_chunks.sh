#!/bin/bash
# Split a full runlist into chunk runlists of roughly <chunk_gb> GB of input each
# (files kept in runlist order, never split).
#   usage: make_chunks.sh <sigma> [chunk_gb=20]
# Writes $WORK_BASE/<sigma>/chunks/chunk_NNN.dat and prints the number of chunks.
source "$(dirname "$0")/config.sh" || exit 1
sig=$1
chunk_gb=${2:-20}
[ -n "$sig" ] || { echo "usage: $0 <sigma> [chunk_gb]" >&2; exit 1; }
full=$(runlist_for "$sig")
[ -f "$full" ] || { echo "no runlist $full" >&2; exit 1; }

out=$WORK_BASE/$sig/chunks
mkdir -p "$out"
rm -f "$out"/chunk_*.dat

files=$(grep -vE '^\s*(#|$)' "$full" | sed -E 's/^\s+|\s+$//g')
missing=$(for f in $files; do [ -e "$f" ] || echo "$f"; done)
[ -z "$missing" ] || { echo "missing inputs:" $missing >&2; exit 1; }

stat -c '%s %n' $files | awk -v lim=$((chunk_gb * 1000000000)) -v out="$out" '
  { if (n == 0 || (tot > 0 && tot + $1 > lim)) { n++; tot = 0 }
    tot += $1; printf "%s\n", $2 > sprintf("%s/chunk_%03d.dat", out, n - 1) }
  END { print n }'
