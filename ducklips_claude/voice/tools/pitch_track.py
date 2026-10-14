"""F0 contour (autocorrelation, 10 ms hop). usage: pitch_track.py FILE.wav[:start:dur]"""
import sys
import wave

import numpy

parts = sys.argv[1].split(":")
h = wave.open(parts[0])
r = h.getframerate()
d = numpy.frombuffer(h.readframes(h.getnframes()), dtype="<i2").astype(float).reshape(-1, h.getnchannels()).mean(axis=1)
if len(parts) == 3:
    d = d[int(float(parts[1]) * r):int((float(parts[1]) + float(parts[2])) * r)]
size = int(0.04 * r)
out = []
peak = numpy.abs(d).max()
for start in range(0, len(d) - size, int(0.01 * r)):
    seg = d[start:start + size]
    if numpy.sqrt(numpy.mean(seg ** 2)) < peak * 0.03:
        out.append("  .")
        continue
    seg = seg - seg.mean()
    c = numpy.correlate(seg, seg, "full")[size - 1:]
    lo, hi = int(r / 350), int(r / 60)
    lag = lo + int(numpy.argmax(c[lo:hi]))
    out.append("%3.0f" % (r / lag) if c[lag] / c[0] > 0.45 else "  ?")
print(" ".join(out))
