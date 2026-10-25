"""Tonal lines at exact kHz multiples in a stretch: level of the 1, 2, 3 kHz
bins (+-15 Hz) re the median level of the surrounding 200-Hz neighbourhood
(dB), and the stretch's level re the file's loudest 50 ms (dB).
usage: tones.py WAV START END"""
import sys, numpy
from scipy.io import wavfile
r, x = wavfile.read(sys.argv[1]); x = x.astype(float)
a, b = float(sys.argv[2]), float(sys.argv[3])
seg = x[int(a * r):int(b * r)]
s = numpy.abs(numpy.fft.rfft(seg * numpy.hanning(len(seg)))) ** 2
f = numpy.fft.rfftfreq(len(seg), 1 / r)
n = int(0.05 * r)
top = max(numpy.sqrt(numpy.mean(x[i:i + n] ** 2)) for i in range(0, len(x) - n, n // 2))
out = []
for k in (1000, 2000, 3000):
    line = s[(f > k - 15) & (f < k + 15)].max()
    around = numpy.median(s[((f > k - 200) & (f < k - 40)) | ((f > k + 40) & (f < k + 200))])
    out.append("%d Hz +%.1f dB (abs %.1f)" % (k, 10 * numpy.log10(line / around), 10 * numpy.log10(line / (top ** 2 * len(seg) ** 2 / 4) + 1e-30)))
print("%s: stretch %.1f dB re loudest; lines %s" % (sys.argv[1].split("/")[-1], 20 * numpy.log10(seg.std() / top + 1e-12), "  ".join(out)))
