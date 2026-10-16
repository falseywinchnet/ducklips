"""Narrow high lines (ringing), time-resolved: 21 ms frames (hop 5 ms);
in each frame, bins above 3 kHz standing 12 dB over the median of +-400 Hz
around them, local maxima; counts frames per 50 Hz band. A ring repeats at
one frequency; noise and moving formants do not.
usage: ring2.py WAV..."""
import sys, numpy
from scipy.io import wavfile
from scipy.ndimage import median_filter
for p in sys.argv[1:]:
    r, x = wavfile.read(p); x = x.astype(float)
    if x.ndim > 1: x = x[:, 0]
    n = 1024 if r < 30000 else 2048 if r < 60000 else 4096
    hop = int(0.005 * r); win = numpy.hanning(n)
    f = numpy.fft.rfftfreq(n, 1 / r); df = f[1] - f[0]
    k = int(800 / df) | 1
    counts = {}
    frames = 0
    for s in range(0, len(x) - n, hop):
        seg = x[s:s + n]
        if numpy.sqrt(numpy.mean(seg ** 2)) < 1e-4 * (numpy.abs(x).max() + 1e-9): continue
        frames += 1
        L = 20 * numpy.log10(numpy.abs(numpy.fft.rfft(seg * win)) + 1e-12)
        prom = L - median_filter(L, size=k)
        for i in numpy.where((f > 3000) & (f < 0.45 * r) & (prom > 12))[0]:
            if prom[i] >= prom[i - 1] and prom[i] >= prom[i + 1]:
                b = int(f[i] // 50) * 50; counts[b] = counts.get(b, 0) + 1
    top = sorted(counts.items(), key=lambda kv: -kv[1])[:6]
    print("%-24s frames %4d  " % (p.split("/")[-1][:24], frames) + "  ".join("%5d Hz %3d" % (b, c) for b, c in top))
