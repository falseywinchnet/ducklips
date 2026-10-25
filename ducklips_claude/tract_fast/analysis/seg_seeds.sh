#!/bin/bash
# Segment-engine identification for each seed given, with say options from OPTS.
cd ~/ducklips/voice
for s in "$@"; do
  ENGINE=tract SAY_ARGS="--seed $s $OPTS" ~/miniforge3/bin/python3 tools/render_arctic.py --analyse >/dev/null 2>&1
  out=$(~/miniforge3/bin/python3 tools/phone_test.py --views stft,cq,causal,wigner 2>/dev/null)
  echo "[seed $s $OPTS] $(echo "$out" | tail -1 | awk "{print \$NF}")"
done
