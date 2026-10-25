#!/bin/bash
# Closure share of voiced cycles per phone over arctic a0001-a0015 for speak env variants.
src=/tmp/es/clo/scores
mkdir -p $src; [ -f $src/arctic_a0015.score ] || cp ~/ducklips/voice/out/arctic/arctic_a00{01,02,03,04,05,06,07,08,09,10,11,12,13,14,15}.score $src/
while [ $# -gt 1 ]; do
  name=$1; envs=$2; shift 2
  d=/tmp/es/clo/$name; mkdir -p $d; cp $src/*.score $d/
  for i in $(seq -w 1 15); do printf "%s\n" "arctic_a00$i"; done | xargs -P 8 -n 1 bash /tmp/clojob.sh $d "$envs" 2>/dev/null
  echo "== $name [$envs]"; ~/miniforge3/bin/python3 /tmp/closure.py $d/*.score | head -14 | tr "\n" ";"; echo
done
