"""Harmonic-to-noise ratio by band in steady voiced frames: in each 40 ms
frame with a clear period (autocorrelation over 0.6), the energy within a
quarter of F0 of each harmonic against the energy in the quarter-F0 middle
between harmonics. Median dB per band over all such frames. Real voices are
noisy between their high harmonics; a too-clean source reads high there.
usage: hnr.py LABEL GLOB [LABEL GLOB ...]"""
import sys, glob, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly

RATE = 16000
BANDS = [(100, 1000), (1000, 2000), (2000, 3000), (3000, 4500), (4500, 7000)]


def frames_of(path):
    r, x = wavfile.read(path)
    x = x.astype(float)
    if x.ndim > 1:
        x = x[:, 0]
    if r != RATE:
        x = resample_poly(x, RATE, r)
    n = int(0.04 * RATE)
    frames = [x[i:i + n] for i in range(0, len(x) - n, n // 2)]
    if not frames:
        return
    level = numpy.array([numpy.sqrt(numpy.mean(f ** 2)) for f in frames])
    top = level.max()
    for f, l in zip(frames, level):
        if l < top * 10 ** (-20 / 20):
            continue
        s = f - f.mean()
        c = numpy.correlate(s, s, "full")[n - 1:]
        lo, hi = int(RATE / 300), int(RATE / 70)
        lag = lo + numpy.argmax(c[lo:hi])
        if c[0] <= 0 or c[lag] < 0.6 * c[0]:
            continue
        # parabolic refinement of the lag
        a, b, d = c[lag - 1], c[lag], c[lag + 1]
        shift = 0.5 * (a - d) / (a - 2 * b + d) if (a - 2 * b + d) != 0 else 0.0
        yield s, RATE / (lag + shift)


def bands_hnr(paths):
    rows = []
    for p in paths:
        for s, f0 in frames_of(p):
            m = 8 * len(s)
            spec = numpy.abs(numpy.fft.rfft(s * numpy.hanning(len(s)), m)) ** 2
            freqs = numpy.fft.rfftfreq(m, 1 / RATE)
            phase = (freqs / f0) % 1.0
            on = (phase < 0.125) | (phase > 0.875)
            off = (phase > 0.375) & (phase < 0.625)
            row = []
            for lo, hi in BANDS:
                sel = (freqs >= lo) & (freqs < hi)
                row.append(10 * numpy.log10((spec[sel & on].mean() + 1e-30) / (spec[sel & off].mean() + 1e-30)))
            rows.append(row)
    return numpy.array(rows)


for label, pattern in zip(sys.argv[1::2], sys.argv[2::2]):
    rows = bands_hnr(sorted(glob.glob(pattern)))
    med = numpy.median(rows, axis=0)
    print("%-8s frames %5d  hnr dB %s" % (label, len(rows), "  ".join(
        "%d-%dk:%5.1f" % (lo // 1000, hi // 1000, v) if lo >= 1000 else "<1k:%5.1f" % v for (lo, hi), v in zip(BANDS, med))))
