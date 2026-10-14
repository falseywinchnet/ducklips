#!/bin/bash
# Render the 29 scores per option set; pooled pitch/voicing, LTAS vs bdl, H1-H2, chosen phone levels.
# usage: full29.sh NAME "options" [NAME "options" ...]
R=$HOME/ducklips/voice/ref/arctic_bdl; P=~/miniforge3/bin/python3
args=("$@")
bash /tmp/r29.sh "${args[@]}"
while [ $# -gt 0 ]; do
  n=$1; shift 2
  echo "== $n"
  $P /tmp/ltas_diff.py /tmp/bdl29.txt /tmp/es/rv_$n/*.wav 2>/dev/null
  (cd /tmp && $P h1h2.py ours /tmp/es/rs /tmp/es/rv_$n 2>/dev/null)
  (cd /tmp && $P phone_level.py $R/lab $R/wav /tmp/es/rs /tmp/es/rv_$n 2>/dev/null | grep -E "^(iy|aa|w |b |d |g |v |dh|z |m |n )" | awk '{printf "%s %s  ", $1, $7}'); echo
done
