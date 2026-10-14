"""Spectrum (dB, 500 Hz bands 1-8 kHz, re its own max) of the middle 50% of spans.
usage: fricspec.py WAV t0 t1 [WAV t0 t1 ...]"""
import sys, numpy, librosa
a = sys.argv[1:]
print("band kHz     " + " ".join("%4.1f" % (b / 1000) for b in range(1000, 8000, 500)))
for i in range(0, len(a), 3):
    y, sr = librosa.load(a[i], sr=16000); s0, s1 = float(a[i + 1]), float(a[i + 2]); d = s1 - s0
    s = y[int((s0 + 0.25 * d) * sr):int((s1 - 0.25 * d) * sr)]
    S = numpy.abs(numpy.fft.rfft(s * numpy.hanning(len(s)))) ** 2; f = numpy.fft.rfftfreq(len(s), 1 / sr)
    b = numpy.array([10 * numpy.log10(S[(f >= x) & (f < x + 500)].sum() + 1e-12) for x in range(1000, 8000, 500)])
    print("%-10s %.2f " % (a[i].split("/")[-1][:10], s0) + " ".join("%4.0f" % v for v in b - b.max()), " peak %.1f kHz" % (1 + 0.5 * numpy.argmax(b)))
