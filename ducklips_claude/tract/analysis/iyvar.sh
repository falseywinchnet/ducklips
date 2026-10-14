#!/bin/bash
# Re-plan a0001-a0015 with IY varea variants, render, report closure share and IY formants.
R=~/ducklips
for va in "$@"; do
  t=/tmp/es/clo/phones_iy$va.txt
  sed "s/^IY inherent 200 vplace 11.95 varea 0.138/IY inherent 200 vplace 11.95 varea $va/" $R/tract/data/phones_tract.txt > $t
  d=/tmp/es/clo/iy$va; mkdir -p $d
  for i in $(seq -w 1 15); do
    k=arctic_a00$i; text=$(grep "( $k " $R/voice/ref/arctic_bdl/arctic.data | sed 's/.*"\(.*\)".*/\1/')
    $R/voice/build/say "$text" /dev/null --score $d/$k.score --plan-only --pressure 750 --phones $t --seed 1 >/dev/null 2>&1
  done
  for i in $(seq -w 1 15); do printf "%s\n" "arctic_a00$i"; done | xargs -P 8 -n 1 bash /tmp/clojob.sh $d "" 2>/dev/null
  echo "== IY varea $va"; ~/miniforge3/bin/python3 /tmp/closure.py $d/*.score | grep -E "^(IY|IH|EY|W) " | tr "\n" ";"; echo
done
