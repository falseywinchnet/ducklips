"""LPC formants (median over 30 ms frames, 1.0-2.0 s) and RMS of a held vowel.
usage: lpcf.py WAV LABEL"""
import sys, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, lfilter
from scipy.linalg import solve_toeplitz
r, x = wavfile.read(sys.argv[1]); x = x.astype(float)
rms = 20 * numpy.log10(x[int(1.0 * r):int(2.0 * r)].std() + 1e-9)
y = resample_poly(x, 10000, r)[10000:20000]; y = lfilter([1, -0.63], [1], y)
res = []
for i in range(0, len(y) - 300, 300):
    f = y[i:i + 300] * numpy.hamming(300); c = numpy.correlate(f, f, "full")[299:312]
    a = solve_toeplitz(c[:12], -c[1:13]); roots = numpy.roots(numpy.r_[1, a]); roots = roots[numpy.imag(roots) > 0]
    fr = numpy.angle(roots) * 10000 / (2 * numpy.pi); bw = -numpy.log(numpy.abs(roots)) * 10000 / numpy.pi
    res.append(sorted((a_, b_) for a_, b_ in zip(fr, bw) if a_ > 90 and b_ < 600)[:3])
fs = [numpy.median([q[k][0] for q in res if len(q) > k]) for k in range(3)]
bs = [numpy.median([q[k][1] for q in res if len(q) > k]) for k in range(3)]
print("%-34s F1 %4.0f (B %3.0f) F2 %4.0f (B %3.0f) F3 %4.0f  rms %.1f dB" % (sys.argv[2], fs[0], bs[0], fs[1], bs[1], fs[2], rms))
