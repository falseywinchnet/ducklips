"""Pitch warble: F0 every 5 ms (autocorrelation, 40 ms frames) over voiced
stretches of at least 150 ms; the contour in cents re its own 3 Hz low-pass
(a 120 ms moving average); rms of that residual in 3-8 Hz and 8-20 Hz
(by differencing moving averages), median over stretches.
usage: warble.py WAV..."""
import sys, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly
def f0track(x, r):
    hop = int(0.005 * r); n = int(0.04 * r); out = []
    for s in range(0, len(x) - n, hop):
        seg = x[s:s + n] - x[s:s + n].mean(); e = numpy.dot(seg, seg)
        if e < 1e-9: out.append(0.0); continue
        c = numpy.correlate(seg, seg, "full")[n - 1:]; lo, hi = int(r / 250), int(r / 65)
        k = lo + numpy.argmax(c[lo:hi])
        if c[k] < 0.5 * c[0]: out.append(0.0); continue
        a, b, g = c[k - 1], c[k], c[k + 1]; d = 0.5 * (a - g) / (a - 2 * b + g) if (a - 2 * b + g) != 0 else 0
        out.append(r / (k + d))
    return numpy.array(out)
def ma(v, w):
    k = numpy.ones(w) / w; return numpy.convolve(v, k, "same")
for p in sys.argv[1:]:
    rs = []
    files = [p]
    for f in files:
        r, x = wavfile.read(f); x = x.astype(float)
        if x.ndim > 1: x = x[:, 0]
        x = resample_poly(x, 16000, r); r = 16000
        f0 = f0track(x, r)
        voiced = f0 > 0; i = 0
        while i < len(f0):
            if not voiced[i]: i += 1; continue
            j = i
            while j < len(f0) and voiced[j]: j += 1
            if j - i >= 30:
                c = 1200 * numpy.log2(f0[i:j] / numpy.median(f0[i:j]))
                slow = ma(c, 24); mid = ma(c, 6) - slow; fast = c - ma(c, 6)
                e = slice(12, j - i - 12)
                if j - i > 30: rs.append((numpy.sqrt(numpy.mean(mid[e] ** 2)), numpy.sqrt(numpy.mean(fast[e] ** 2))))
            i = j
    rs = numpy.array(rs)
    print("%-26s stretches %3d  warble 3-8 Hz %.1f cents  8-20 Hz %.1f cents" % (p.split("/")[-1][:26], len(rs), numpy.median(rs[:, 0]), numpy.median(rs[:, 1])))
