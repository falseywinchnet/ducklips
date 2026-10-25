"""Median LPC F1-F3 in the middle 40% of each instance of the given vowels
(at least 70 ms), bdl from labels or ours from scores.
usage: vformants.py bdl LABDIR WAVDIR N VOWELS... | vformants.py ours SCOREDIR WAVDIR 0 VOWELS..."""
import sys, os, numpy
sys.path.insert(0, "/tmp")
from vw import load
from dips import phones_bdl, phones_ours
from scipy.signal import resample_poly, lfilter
from scipy.linalg import solve_toeplitz


def lpc(y):
    out = []
    for i in range(0, len(y) - 250, 100):
        f = y[i:i + 250] * numpy.hamming(250); c = numpy.correlate(f, f, "full")[249:262]
        if c[0] <= 0:
            continue
        a = solve_toeplitz(c[:12], -c[1:13]); roots = numpy.roots(numpy.r_[1, a]); roots = roots[numpy.imag(roots) > 0]
        fr = numpy.angle(roots) * 10000 / (2 * numpy.pi); bw = -numpy.log(numpy.abs(roots)) * 10000 / numpy.pi
        q = sorted((a_, b_) for a_, b_ in zip(fr, bw) if a_ > 150 and b_ < 800)[:3]
        if len(q) == 3:
            out.append([v for pair in q for v in pair])
    return numpy.median(out, axis=0) if out else None


kind, d1, d2, n = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
want = set(sys.argv[5:])
src = phones_bdl(d1, d2, n) if kind == "bdl" else phones_ours(d1, d2)
got = {}
for wav, phones in src:
    if not os.path.exists(wav):
        continue
    x = resample_poly(load(wav), 10000, 16000); x = lfilter([1, -0.63], [1], x)
    for name, s, e in phones:
        if name in want and e - s >= 0.07:
            a, b = s + 0.3 * (e - s), e - 0.3 * (e - s)
            r = lpc(x[int(a * 10000):int(b * 10000)])
            if r is not None:
                got.setdefault(name, []).append(r)
for name in sorted(got):
    m = numpy.median(got[name], axis=0)
    print("%-6s %-3s n %3d  F1 %4.0f (B %3.0f) F2 %4.0f (B %3.0f) F3 %4.0f (B %3.0f)" % (kind, name, len(got[name]), m[0], m[1], m[2], m[3], m[4], m[5]))
