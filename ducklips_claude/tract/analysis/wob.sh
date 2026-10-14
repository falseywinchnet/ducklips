#!/bin/bash
# Render variants of the 29 ARCTIC scores and print their wobble.
# usage: wob.sh NAME "options" [NAME "options" ...]
while [ $# -gt 0 ]; do
  name=$1; opts=$2; shift 2
  out=/tmp/es/rv_$name; mkdir -p $out
  export OPTS="$opts" OUT=$out
  ls /tmp/es/rs/*.score | xargs -P 8 -I{} bash -c 'b=$(basename {} .score); ~/ducklips/tract/build/sheet_speak {} $OUT/$b.wav $OPTS >/dev/null 2>&1'
  cd /tmp && ~/miniforge3/bin/python3 f0mod.py "$name" "$out/*.wav" 2>/dev/null
done
