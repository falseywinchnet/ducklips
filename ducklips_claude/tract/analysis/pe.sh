#!/bin/bash
# Pitch error, swing and vowel voicing on the two test sentences per option set.
# usage: pe.sh NAME "options" [NAME "options" ...]
cd /tmp/es; P=~/miniforge3/bin/python3; S=~/ducklips/tract/build/sheet_speak
while [ $# -gt 0 ]; do
  n=$1; o=$2; shift 2
  $S v1t.score pe1_$n.wav $o >/dev/null 2>&1; $S v2t.score pe2_$n.wav $o >/dev/null 2>&1
  echo "== $n [$o]"
  $P /tmp/pitch_error.py v1t.score pe1_$n.wav 2>/dev/null
  $P /tmp/pitch_error.py v2t.score pe2_$n.wav 2>/dev/null
  $P /tmp/vowel_voicing.py v1t.score pe1_$n.wav 2>/dev/null | tail -1
  $P /tmp/vowel_voicing.py v2t.score pe2_$n.wav 2>/dev/null | tail -1
done
