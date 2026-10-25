#!/bin/bash
# One-seed identification with a modified phone table; prints overall and chosen phones.
# usage: ptry.sh NAME "sed expression on the table" "phones regex"
name=$1; expr=$2; want=$3
sed -E "$expr" ~/ducklips/tract/data/phones_sheet.txt > /tmp/es/pt_$name.txt
cd ~/ducklips/voice
ENGINE=sheet PHONES=/tmp/es/pt_$name.txt SAY_ARGS="--seed 1" ~/miniforge3/bin/python3 tools/render_arctic.py --analyse >/dev/null 2>&1
out=$(~/miniforge3/bin/python3 tools/phone_test.py --views stft,cq,causal,wigner 2>/dev/null)
echo "$name: $(echo "$out" | tail -1 | awk '{print $NF}')  $(echo "$out" | grep 'per-phone' -A1 | tail -1 | tr ',' '\n' | grep -E "^ *($want) " | tr '\n' ',')"
