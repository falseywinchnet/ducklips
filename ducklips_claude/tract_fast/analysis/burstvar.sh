#!/bin/bash
# Re-speak the 40-sentence corpus with speak env settings; print P/T/K burst rows (ours) vs bdl.
for pair in "$@"; do
  name=${pair%%:*}; envs=${pair#*:}
  mkdir -p /tmp/es/cv_$name; cd /tmp/es/corpus
  export ENVS="$envs" OUT=/tmp/es/cv_$name
  ls *.score | xargs -P 8 -I{} bash -c 'b=$(basename {} .score); cp {} $OUT/; env $ENVS ~/ducklips/tract/build/speak {} $OUT/$b.wav >/dev/null 2>&1'
  sed "s|/tmp/es/corpus/|$OUT/|g" /tmp/bursts2.py > /tmp/bursts_v.py
  echo "== $name [$envs]"; (cd ~/ducklips/voice && ~/miniforge3/bin/python3 /tmp/bursts_v.py 2>/dev/null | grep -E "^(P|T|K) ")
  rm -rf $OUT
done
