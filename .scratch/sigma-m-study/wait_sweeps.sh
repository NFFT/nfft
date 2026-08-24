#!/bin/sh
# Emit one line when the long-double sweep (the slowest) has all its rows.
OUT="$1"
WANT="${2:-4680}"
while true; do
  n=$(wc -l < "$OUT/sweep-l.csv" 2>/dev/null || echo 0)
  if [ "$n" -ge "$WANT" ]; then
    echo "sweeps complete: long-double has $n rows"
    break
  fi
  sleep 30
done
