#!/bin/bash
# From the laptop: speak English text on the Mini (where the builds live,
# ~/ducklips) with tools/talk.sh, copy the WAV back here and play it.
# usage: talk_remote.sh "text" [OUT.wav] [seed]
#   MINI overrides the ssh host (default m4mini); PLAY=0 skips playback.
set -euo pipefail
if [ $# -lt 1 ]; then
  echo "usage: talk_remote.sh \"text\" [OUT.wav] [seed]" >&2
  exit 2
fi
host="${MINI:-m4mini}"
out="${2:-talk.wav}"
seed="${3:-73021}"
remote="/tmp/talk_remote_$$.wav"
text_quoted=$(printf '%q' "$1")
ssh "$host" "~/ducklips/tract/tools/talk.sh $text_quoted $remote $seed"
scp -q "$host:$remote" "$out"
ssh "$host" "rm -f $remote"
echo "$out"
if [ "${PLAY:-1}" != "0" ] && command -v afplay >/dev/null 2>&1; then
  afplay "$out"
fi
