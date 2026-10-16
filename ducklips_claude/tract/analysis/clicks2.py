"""Local clicks: the signal high-passed at 5 kHz; samples whose magnitude
exceeds K times the rms of the high-passed signal over the surrounding 20 ms
(excluding 1 ms around the sample). Prints count per second and the times.
usage: clicks2.py WAV [K]"""
import sys, numpy
from scipy.io import wavfile
from scipy.signal import butter, sosfiltfilt
r, x = wavfile.read(sys.argv[1]); x = x.astype(float)
if x.ndim > 1: x = x[:, 0]
K = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
h = sosfiltfilt(butter(4, 5000, btype="high", fs=r, output="sos"), x)
w = int(0.02 * r); e = numpy.sqrt(numpy.convolve(h ** 2, numpy.ones(w) / w, "same")) + 1e-9
level = numpy.sqrt(numpy.convolve(x ** 2, numpy.ones(w) / w, "same"))
loud = level > level.max() * 10 ** (-40 / 20)
idx = numpy.where((numpy.abs(h) > K * e) & loud)[0]
times = []
for i in idx:
    t = i / r
    if not times or t - times[-1] > 0.003: times.append(t)
print("%s: %d clicks (%.1f/s)" % (sys.argv[1].split("/")[-1], len(times), len(times) / (len(x) / r)))
print(" ".join("%.4f" % t for t in times))
