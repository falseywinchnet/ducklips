"""Long-term spectrum of voiced loud frames of the given files together,
against a saved reference row (bands as ltas.py, dB re 300-1000 Hz):
prints the per-band difference and the mean |difference| over 160 Hz-10 kHz.
usage: ltas_diff.py REF.txt WAV... ; ltas_diff.py --save REF.txt WAV..."""
import sys, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly
BANDS = [(80, 160), (160, 300), (300, 600), (600, 1000), (1000, 1600), (1600, 2500), (2500, 4000), (4000, 6300), (6300, 10000)]
def spectrum(paths):
    acc = None
    for p in paths:
        r, x = wavfile.read(p); x = x.astype(float)
        if x.ndim > 1: x = x[:, 0]
        if r != 32000: x = resample_poly(x, 32000, r); r = 32000
        n = int(0.03 * r); frames = [x[i:i + n] for i in range(0, len(x) - n, n // 2)]
        level = numpy.array([numpy.sqrt(numpy.mean(f ** 2)) for f in frames]); top = level.max()
        if acc is None: acc = numpy.zeros(n // 2 + 1)
        for f, l in zip(frames, level):
            if l < top * 10 ** (-25 / 20): continue
            s = f - f.mean(); c = numpy.correlate(s, s, "full")[n - 1:]; lo, hi = int(r / 300), int(r / 60)
            if c[0] <= 0 or c[lo + numpy.argmax(c[lo:hi])] < 0.3 * c[0]: continue
            acc += numpy.abs(numpy.fft.rfft(s * numpy.hanning(n))) ** 2
    freqs = numpy.fft.rfftfreq(len(acc) * 2 - 2, 1 / 32000)
    ref = 10 * numpy.log10(acc[(freqs >= 300) & (freqs < 1000)].sum() + 1e-30)
    return numpy.array([10 * numpy.log10(acc[(freqs >= a) & (freqs < b)].sum() + 1e-30) - ref for a, b in BANDS])
if sys.argv[1] == "--save":
    numpy.savetxt(sys.argv[2], spectrum(sys.argv[3:])); sys.exit()
reference = numpy.loadtxt(sys.argv[1]); ours = spectrum(sys.argv[2:]); d = ours - reference
print("ltas diff %s  mean|d| %.1f" % (" ".join("%+.0f" % v for v in d), numpy.mean(numpy.abs(d[1:]))))
