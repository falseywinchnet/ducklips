#!/bin/bash
# Three-seed identification for each phone table given (sheet defaults).
cd ~/ducklips/voice
for table in "$@"; do
  total=0; line=""
  for s in 1 2 3; do
    ENGINE=sheet PHONES="$table" SAY_ARGS="--seed $s" ~/miniforge3/bin/python3 tools/render_arctic.py --analyse >/dev/null 2>&1
    v=$(~/miniforge3/bin/python3 tools/phone_test.py --views stft,cq,causal,wigner | tail -1 | awk "{print \$NF}")
    line="$line $v"; total=$(python3 -c "print($total + $v)")
  done
  printf "%-44s %s  mean %.1f\n" "[$table]" "$line" $(python3 -c "print($total / 3)")
done
