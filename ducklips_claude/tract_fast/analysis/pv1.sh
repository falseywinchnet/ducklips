#!/bin/bash
# Per-vowel voicing on sentence 1 for option sets.
# usage: pv1.sh "options" ...
cd /tmp/es; files=""
for o in "$@"; do
  k=$((k+1)); n=pv_$k
  ${SHEET_BIN:-~/ducklips/tract/build/sheet_speak} v1t.score $n.wav $o >/dev/null 2>&1; files="$files $n.wav"
  echo "col: [$o]"
done
~/miniforge3/bin/python3 /tmp/pervowel.py v1t.score $files 2>/dev/null
