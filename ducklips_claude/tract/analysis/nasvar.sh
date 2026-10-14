#!/bin/bash
# Nasal band levels over arctic a0001-a0030 for speak environment variants.
# usage: nasvar.sh NAME "ENV=..." ...  (pairs)
src=~/ducklips/voice/out/arctic
while [ $# -gt 1 ]; do
  name=$1; envs=$2; shift 2
  d=/tmp/es/nas/$name; mkdir -p $d
  for i in $(seq -w 1 30); do
    k=arctic_a00$i; cp $src/$k.score $d/$k.score
    echo "env ROOM=0 $envs /tmp/es/tbuild/speak $d/$k.score $d/$k.wav 1 >/dev/null 2>&1"
  done | xargs -P 8 -I{} bash -c "{}"
  ~/miniforge3/bin/python3 /tmp/nasalband.py ours $d | sed "s/^/$name [$envs] /"
done
