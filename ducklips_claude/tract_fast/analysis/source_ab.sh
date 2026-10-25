#!/bin/bash
# Source A/B: the same three sentences through our self-oscillating folds and
# through the kinematic folds. usage: source_ab.sh OUTDIR
R=~/ducklips; out=$1; mkdir -p $out; cd $out
T1="Hello, my name is Papa, and I speak with a body."
T2="The quick brown fox jumps over the lazy dog."
T3="I hope I sound a little more human today."
for k in 1 2 3; do
  t=T$k
  (cd $R/voice && ./build/say "${!t}" /dev/null --score $out/s$k.score --plan-only --pressure 750 --phones $R/tract/data/phones_tract.txt >/dev/null 2>&1)
done
render() {  # name "speak env"
  name=$1; envs=$2; parts=""
  for k in 1 2 3; do env $envs $R/tract/build/speak s$k.score $name.$k.wav 73021 >/dev/null 2>&1; parts="$parts $name.$k.wav"; done
  ~/miniforge3/bin/python3 - $name $parts <<'PY'
import sys, numpy
from scipy.io import wavfile
parts = []
for f in sys.argv[2:]:
    r, x = wavfile.read(f); parts += [x, numpy.zeros(int(0.5 * r), x.dtype)]
wavfile.write(sys.argv[1] + ".wav", r, numpy.concatenate(parts))
PY
  rm -f $name.?.wav
}
render A_self_oscillating_folds ""
render B_moved_folds "KINEMATIC=1 VOICEPOST=0.015,0.01"
render C_moved_folds_breathier "KINEMATIC=1 VOICEPOST=0.015,0.01 CHINK=0.03"
ls -la *.wav
