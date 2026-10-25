"""Pitch and loudness wobble inside voiced stretches: F0 (normalized
autocorrelation, 30 ms window, 5 ms hop) and level, over voiced runs of at
least 150 ms. Each run's track (cents re its mean; dB) minus a straight-line
trend; then the RMS of what remains in 3-8 Hz and 8-20 Hz (FFT band sums).
Medians over runs. Also median |cycle-to-cycle F0 step| in cents per 5 ms.
usage: f0mod.py LABEL GLOB [LABEL GLOB ...]"""
import sys, glob, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly

RATE = 16000
HOP = 80
WIN = 480


def track(path):
    r, x = wavfile.read(path)
    x = x.astype(float)
    if x.ndim > 1:
        x = x[:, 0]
    if r != RATE:
        x = resample_poly(x, RATE, r)
    lo, hi = int(RATE / 300), int(RATE / 60)
    f0, lev = [], []
    top = numpy.sqrt(numpy.mean(x ** 2)) * 4
    for i in range(0, len(x) - WIN - hi, HOP):
        a = x[i:i + WIN] - x[i:i + WIN].mean()
        energy = numpy.dot(a, a)
        level = 10 * numpy.log10(energy / WIN + 1e-9)
        best, lag = 0.0, 0
        for k in range(lo, hi):
            b = x[i + k:i + k + WIN] - x[i + k:i + k + WIN].mean()
            c = numpy.dot(a, b) / numpy.sqrt(energy * numpy.dot(b, b) + 1e-30)
            if c > best:
                best, lag = c, k
        if best > 0.6 and energy / WIN > (top / 30) ** 2:
            cs = []
            for k in (lag - 1, lag, lag + 1):
                b = x[i + k:i + k + WIN] - x[i + k:i + k + WIN].mean()
                cs.append(numpy.dot(a, b) / numpy.sqrt(energy * numpy.dot(b, b) + 1e-30))
            den = cs[0] - 2 * cs[1] + cs[2]
            shift = 0.5 * (cs[0] - cs[2]) / den if den < 0 else 0.0
            f0.append(RATE / (lag + shift))
        else:
            f0.append(0.0)
        lev.append(level)
    return numpy.array(f0), numpy.array(lev)


def band_rms(v, lo, hi):
    t = numpy.arange(len(v))
    v = v - numpy.polyval(numpy.polyfit(t, v, 1), t)
    s = numpy.fft.rfft(v * numpy.hanning(len(v)))
    f = numpy.fft.rfftfreq(len(v), HOP / RATE)
    sel = (f >= lo) & (f < hi)
    return numpy.sqrt(2 * numpy.sum(numpy.abs(s[sel]) ** 2)) / numpy.sqrt(numpy.sum(numpy.hanning(len(v)) ** 2))


def measure(paths):
    p38, p820, a38, a820, steps = [], [], [], [], []
    for p in paths:
        f0, lev = track(p)
        voiced = f0 > 0
        i = 0
        while i < len(f0):
            if not voiced[i]:
                i += 1
                continue
            j = i
            while j < len(f0) and voiced[j]:
                j += 1
            c = 1200 * numpy.log2(f0[i:j] / numpy.median(f0[i:j])) if j > i else None
            if j - i >= 30 and numpy.max(numpy.abs(numpy.diff(c))) < 250:
                steps.extend(numpy.abs(numpy.diff(c)))
                p38.append(band_rms(c, 3, 8))
                p820.append(band_rms(c, 8, 20))
                a38.append(band_rms(lev[i:j], 3, 8))
                a820.append(band_rms(lev[i:j], 8, 20))
            i = j
    return len(p38), [numpy.median(v) for v in (p38, p820, a38, a820, steps)]


for label, pattern in zip(sys.argv[1::2], sys.argv[2::2]):
    runs, (p38, p820, a38, a820, step) = measure(sorted(glob.glob(pattern)))
    print("%-8s runs %4d  pitch 3-8Hz %5.1f c  8-20Hz %5.1f c  level 3-8Hz %.2f dB  8-20Hz %.2f dB  step %.1f c" % (
        label, runs, p38, p820, a38, a820, step))
