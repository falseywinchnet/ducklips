#!/bin/bash
# Render the 29 ARCTIC scores per option set (binary SHEET_BIN) and pool pitch error.
# usage: r29.sh NAME "options" [NAME "options" ...]
bin=${SHEET_BIN:-$HOME/ducklips/tract/build/sheet_speak}
dirs=""
while [ $# -gt 0 ]; do
  name=$1; opts=$2; shift 2
  out=/tmp/es/rv_$name; mkdir -p $out
  export OPTS="$opts" OUT=$out BIN=$bin
  ls /tmp/es/rs/*.score | xargs -P 8 -I{} bash -c 'b=$(basename {} .score); $BIN {} $OUT/$b.wav $OPTS >/dev/null 2>&1'
  dirs="$dirs $out"
done
~/miniforge3/bin/python3 /tmp/pe29.py /tmp/es/rs $dirs 2>/dev/null
