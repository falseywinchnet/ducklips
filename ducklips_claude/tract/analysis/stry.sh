#!/bin/bash
# One-seed identification with sheet options; prints overall and rows of chosen phones.
cd ~/ducklips/voice
for o in "$@"; do
  ENGINE=sheet SHEET="$o" SAY_ARGS="--seed 1" ~/miniforge3/bin/python3 tools/render_arctic.py --analyse >/dev/null 2>&1
  out=$(ROWS=${ROWS:-} ~/miniforge3/bin/python3 tools/phone_test.py --views stft,cq,causal,wigner 2>/dev/null)
  echo "[$o] $(echo "$out" | tail -1 | awk "{print \$NF}")  $(echo "$out" | grep "per-phone" -A1 | tail -1 | tr "," "\n" | grep -E "^ *(${WANT:-z|v|dh}) " | tr "\n" ",")"
done
