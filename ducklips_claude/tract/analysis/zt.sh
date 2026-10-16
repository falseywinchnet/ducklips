#!/bin/bash
# /z/ spectra (4-8 kHz mean re its own peak) for voiced-fricative gaps.
cd /tmp/es/zt
for g in 0.012 0.025 0.04 0.06; do
  (cd ~/ducklips/voice && ./build/say "His zoo was easy. Is she busy? The size was amazing." /dev/null --score /tmp/es/zt/z_$g.score --plan-only --pressure 750 --phones ~/ducklips/tract/data/phones_tract.txt --organs ~/ducklips/tract/data/organs.txt --voiced-fricative-gap $g >/dev/null)
  ROOM=0 KINEMATIC=1 VOICEPOST=0.015,0.01 ~/ducklips/tract/build/speak z_$g.score z_$g.wav 2>/dev/null
  grep "^# Z" z_$g.score | while read h n s d r; do ~/miniforge3/bin/python3 /tmp/fricspec.py z_$g.wav $s $(echo "$s + $d" | bc) | tail -1; done | awk -v g=$g '{h=0; for(i=10;i<=16;i++) h+=$i; s+=h/7; n++} END {printf "gap %s: z 4-7.5 kHz mean %.1f dB re peak (n %d)\n", g, s/n, n}'
done
