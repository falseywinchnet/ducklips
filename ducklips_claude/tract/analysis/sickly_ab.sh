#!/bin/bash
# A/B set for "sickly": the same two sentences under one change each.
# usage: sickly_ab.sh OUTDIR
R=~/ducklips; out=$1; mkdir -p $out; cd $out
T1="I hope I sound a little more human today."
T2="Hello, my name is Papa, and I speak with a body."
where=$R/voice   # say finds its learned pitch and duration models here
render() {  # name "say options" "speak env"
  name=$1; sayopts=$2; envs=$3; parts=""
  for k in 1 2; do
    t=T$k
    (cd $where && $R/voice/build/say "${!t}" /dev/null --score $out/$name.$k.score --plan-only --pressure 750 --phones $R/tract/data/phones_tract.txt $sayopts >/dev/null 2>&1)
    env $envs $R/tract/build/speak $name.$k.score $name.$k.wav 73021 >/dev/null 2>&1
    parts="$parts $name.$k.wav"
  done
  ~/miniforge3/bin/python3 - $name $parts <<'PY'
import sys, numpy
from scipy.io import wavfile
parts = []
for f in sys.argv[2:]:
    r, x = wavfile.read(f); parts += [x, numpy.zeros(int(0.5 * r), x.dtype)]
wavfile.write(sys.argv[1] + ".wav", r, numpy.concatenate(parts))
PY
  rm -f $name.?.wav $name.?.score
}
where=/tmp
render A0_as_heard_hand_rules "" ""
where=$R/voice
render A_learned_models "" ""
render B_firmer_folds "" "CORDDAMP=0.15,0.8"
render C_stronger_breath "--pressure 950" ""
render D_perfectly_steady "" "TREMOR=0 BREATHVAR=0 TENSIONVAR=0 POSTUREVAR=0 WALLVAR=0 WANDER=0"
render E_higher_pitch "--pitch-median 140" ""
ls -la
