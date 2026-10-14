"""Band levels (dB re 0-1.5 kHz) of a segment. usage: bands.py FILE.wav start dur"""
import sys
import wave

import numpy

h = wave.open(sys.argv[1])
rate = h.getframerate()
d = numpy.frombuffer(h.readframes(h.getnframes()), dtype="<i2").astype(float).reshape(-1, h.getnchannels()).mean(axis=1)
seg = d[int(float(sys.argv[2]) * rate):int((float(sys.argv[2]) + float(sys.argv[3])) * rate)]
S = numpy.abs(numpy.fft.rfft(seg * numpy.hanning(len(seg)))) ** 2
f = numpy.fft.rfftfreq(len(seg), 1 / rate)


def band(a, b):
    return 10 * numpy.log10(S[(f >= a) & (f < b)].sum() + 1e-20)


ref = band(0, 1500)
print("%-28s" % (sys.argv[1].split("/")[-1] + "@" + sys.argv[2]),
      " ".join("%d-%d:%.0f" % (a, b, band(a, b) - ref) for a, b in [(1500, 3000), (3000, 4500), (4500, 6000), (6000, 8000), (8000, 12000)]))
