#!/bin/bash
# Three-seed identification for each sheet option set given; mean last.
cd ~/ducklips/voice
for opts in "$@"; do
  total=0
  line=""
  for s in 1 2 3; do
    ENGINE=sheet SHEET="$opts" SAY_ARGS="--seed $s" ~/miniforge3/bin/python3 tools/render_arctic.py --analyse >/dev/null 2>&1
    v=$(~/miniforge3/bin/python3 tools/phone_test.py --views stft,cq,causal,wigner | tail -1 | awk '{print $NF}')
    line="$line $v"
    total=$(python3 -c "print($total + $v)")
  done
  printf "%-44s %s  mean %.1f\n" "[$opts]" "$line" $(python3 -c "print($total / 3)")
done
