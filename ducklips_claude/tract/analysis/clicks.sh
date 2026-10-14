#!/bin/bash
# Peak over speech rms (dB) and single-sample spikes for NATO word scores under option sets.
cd /tmp/es; S=~/ducklips/tract/build/sheet_speak
for o in "$@"; do
  for w in b q d k p; do GAIN=0.0005 $S /tmp/nato_ours_$w.wav.score ck_$w.wav $o --room 0 >/dev/null 2>&1; done
  ~/miniforge3/bin/python3 - "$o" <<'PY'
import numpy, sys
from scipy.io import wavfile
out=[]
for w in "bqdkp":
    r,x=wavfile.read("/tmp/es/ck_%s.wav"%w); x=x.astype(float)
    rms=numpy.sqrt(numpy.mean(x[numpy.abs(x)>numpy.abs(x).max()*0.05]**2))
    d2=numpy.abs(numpy.diff(x,2)); spikes=(d2>8*numpy.median(d2[d2>0])*20).sum() if (d2>0).any() else 0
    out.append("%s %.1f/%d" % (w, 20*numpy.log10(numpy.abs(x).max()/rms), spikes))
print("[%s] peak re rms dB / spikes: %s" % (sys.argv[1], "  ".join(out)))
PY
done
