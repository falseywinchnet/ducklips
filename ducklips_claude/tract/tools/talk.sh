#!/bin/bash
# Speak English text with the current main line: the planner (say, breath
# 750 Pa, the segment engine's phone table) and the segment-module engine
# (speak, with its room). Run on the machine that holds the builds (the
# Mini: ~/ducklips).
# usage: talk.sh "text" OUT.wav [seed]
#   VOICE_ROOT overrides the ducklips root (default: this script's ../..).
set -euo pipefail
if [ $# -lt 2 ]; then
  echo "usage: talk.sh \"text\" OUT.wav [seed]" >&2
  exit 2
fi
root="${VOICE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
seed="${3:-73021}"
score="$(mktemp -t talk).score"
"$root/voice/build/say" "$1" /dev/null --score "$score" --plan-only --pressure 750 \
  --phones "$root/tract/data/phones_tract.txt" >/dev/null
"$root/tract/build/speak" "$score" "$2" "$seed"
rm -f "$score"
