#!/bin/bash
# Render the ARCTIC set under speak env variants; print vowel band
# differences (ours - bdl) averaged over the main vowels, and selected phones.
cd ~/ducklips/voice
while [ $# -gt 1 ]; do
  name=$1; envs=$2; shift 2
  ENGINE=tract SAY_ARGS="--seed 1 --organs $HOME/ducklips/tract/data/organs.txt" SPEAK_ENV="KINEMATIC=1 VOICEPOST=0.015,0.01 $envs" ~/miniforge3/bin/python3 tools/render_arctic.py --analyse >/dev/null 2>&1
  ~/miniforge3/bin/python3 /tmp/phonebands.py > /tmp/es/bands_$name.txt
  echo "== $name [$envs]"
  awk '$1 ~ /^(iy|ih|eh|ae|aa|ah|ay|aw|ey|ao)$/ {for(i=3;i<=9;i++) s[i]+=$i; n++} END {printf "vowels mean ours-bdl:"; for(i=3;i<=9;i++) printf " %+5.1f", s[i]/n; print ""}' /tmp/es/bands_$name.txt
  grep -E "^(l|sh|dh|jh) " /tmp/es/bands_$name.txt
done
