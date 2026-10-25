#!/bin/bash
# Re-speak the four planned scores of four/p1_* with environment settings for speak; assemble.
# usage: four2.sh NAME [VAR=value ...]
name=$1; shift
cd /tmp/es/four
for i in 0 1 2 3; do env "$@" ~/ducklips/tract/build/speak p1_$i.score ${name}_$i.wav >/dev/null 2>&1; done
~/miniforge3/bin/python3 - $name <<'PY'
import numpy, sys
from scipy.io import wavfile
from scipy.signal import resample_poly
name = sys.argv[1]
def load(p):
    r, x = wavfile.read(p); x = x.astype(float)
    if r != 48000: x = resample_poly(x, 48000, r)
    f = numpy.array([numpy.sqrt(numpy.mean(x[i:i + 960]**2)) for i in range(0, len(x) - 960, 960)])
    ref = numpy.mean(numpy.sort(f)[-max(1, len(f) * 3 // 10):])
    return numpy.clip(x / ref * 0.12, -0.98, 0.98)
gap = numpy.zeros(int(0.6 * 48000)); parts = []
for i in range(4): parts += [load("/tmp/es/four/%s_%d.wav" % (name, i)), gap]
wavfile.write("/tmp/es/four/%s.wav" % name, 48000, (numpy.concatenate(parts) * 32767).astype(numpy.int16))
r, x = wavfile.read("/tmp/es/four/%s_0.wav" % name); x = resample_poly(x.astype(float), 16000, r)[:int(2.4 * 16000)]
x = x / numpy.abs(x).max() * 0.8
wavfile.write("/tmp/es/tfv/%s.wav" % name, 16000, (x * 32767).astype(numpy.int16))
PY
~/ducklips/voice/build/tf_analyze /tmp/es/tfv/$name.wav /tmp/es/tfv/$name >/dev/null 2>&1
