#!/bin/bash
# Skin share of the output in voiced frames for option sets (a0001, fixed gain).
cd /tmp/es; S=~/ducklips/tract/build/sheet_speak
for o in "$@"; do
  GAIN=0.0005 $S rs/arctic_a0001.score sk1.wav $o --room 0 >/dev/null 2>&1
  GAIN=0.0005 $S rs/arctic_a0001.score sk0.wav $o --skin 0 --room 0 >/dev/null 2>&1
  ~/miniforge3/bin/python3 - "$o" <<'PY'
import numpy, sys
from scipy.io import wavfile
r,a=wavfile.read("/tmp/es/sk1.wav"); r,b=wavfile.read("/tmp/es/sk0.wav")
d=a.astype(float)-b.astype(float); b=b.astype(float); n=480
fr=numpy.array([(20*numpy.log10(d[i:i+n].std()+1e-9), 20*numpy.log10(b[i:i+n].std()+1e-9)) for i in range(0,len(d)-n,n)]); live=fr[:,1]>fr[:,1].max()-30
print("[%s] skin re mouth, voiced frames: median %.1f dB, 90th %.1f" % (sys.argv[1], numpy.median(fr[live,0]-fr[live,1]), numpy.percentile(fr[live,0]-fr[live,1],90)))
PY
done
rm -f sk0.wav sk1.wav
