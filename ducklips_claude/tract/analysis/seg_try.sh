#!/bin/bash
# One-seed segment-engine identification per say-option set; overall and chosen phones (WANT regex).
cd ~/ducklips/voice
for o in "$@"; do
  ENGINE=tract SAY_ARGS="--seed 1 $o" ~/miniforge3/bin/python3 tools/render_arctic.py --analyse >/dev/null 2>&1
  out=$(~/miniforge3/bin/python3 tools/phone_test.py --views stft,cq,causal,wigner 2>/dev/null)
  echo "[$o] $(echo "$out" | tail -1 | awk "{print \$NF}")  $(echo "$out" | grep "per-phone" -A1 | tail -1 | tr "," "\n" | grep -E "^ *(${WANT:-hh}) " | tr "\n" ",")"
done
