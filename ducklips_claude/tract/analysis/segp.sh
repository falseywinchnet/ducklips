#!/bin/bash
# Segment engine on m2i's four scores with lungs scaled; assembles and measures modulation.
# usage: segp.sh NAME SCALE [speak seed]
cd /tmp/es; name=$1; scale=$2
for i in 1 2 3 4; do
  awk -v s=$scale '$2=="lungs"{$3=sprintf("%.3f",$3*s)}1' listen/s_$i.score > sp_$i.score
  ~/ducklips/tract/build/speak sp_$i.score sp_$i.wav >/dev/null 2>&1
done
~/miniforge3/bin/python3 - $name <<'PY'
import numpy, sys
from scipy.io import wavfile
from scipy.signal import resample_poly
def load(p):
    r, x = wavfile.read(p); x = x.astype(float)
    if r != 48000: x = resample_poly(x, 48000, r)
    f = numpy.array([numpy.sqrt(numpy.mean(x[i:i + 960]**2)) for i in range(0, len(x) - 960, 960)])
    ref = numpy.mean(numpy.sort(f)[-max(1, len(f) * 3 // 10):])
    return numpy.clip(x / ref * 0.12, -0.98, 0.98)
gap = numpy.zeros(int(0.6 * 48000)); parts = []
for i in range(1, 5): parts += [load("/tmp/es/sp_%d.wav" % i), gap]
wavfile.write("/tmp/es/%s.wav" % sys.argv[1], 48000, (numpy.concatenate(parts) * 32767).astype(numpy.int16))
PY
rm -f sp_*.wav sp_*.score
~/miniforge3/bin/python3 /tmp/vibrator.py $name /tmp/es/$name.wav 2>/dev/null
