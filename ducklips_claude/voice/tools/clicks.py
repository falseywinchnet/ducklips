"""Lists click-like transients: 1 ms peaks far above the surrounding 30 ms RMS.
usage: clicks.py FILE.wav [ratio]"""
import sys
import wave

import numpy

h = wave.open(sys.argv[1])
r = h.getframerate()
d = numpy.frombuffer(h.readframes(h.getnframes()), dtype="<i2").astype(float) / 32768
ratio = float(sys.argv[2]) if len(sys.argv) > 2 else 6.0
hp = numpy.diff(d, prepend=0.0)  # emphasize the step content
step = r // 1000
for i in range(15, len(d) // step - 15):
    peak = numpy.abs(hp[i * step:(i + 1) * step]).max()
    around = numpy.concatenate([hp[(i - 15) * step:(i - 1) * step], hp[(i + 2) * step:(i + 15) * step]])
    rms = numpy.sqrt(numpy.mean(around ** 2)) + 1e-5
    if peak / rms > ratio and peak > 0.01:
        print("%.3f s  peak/rms %.1f" % (i / 1000, peak / rms))
