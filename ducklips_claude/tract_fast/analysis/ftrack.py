"""LPC formant track (F1-F3, every 5 ms, 25 ms window, 10 kHz) between t0 and t1, with score phone labels.
usage: ftrack.py WAV t0 t1 [SCORE]"""
import sys, re, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, lfilter
from scipy.linalg import solve_toeplitz
r, x = wavfile.read(sys.argv[1]); x = x.astype(float)
if x.ndim > 1: x = x[:, 0]
x = lfilter([1, -0.63], [1], resample_poly(x, 10000, r))
t0, t1 = float(sys.argv[2]), float(sys.argv[3])
ph = []
if len(sys.argv) > 4:
    ph = [(m.group(1), float(m.group(2))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sys.argv[4])) if m]
ref = numpy.sqrt(numpy.mean(x ** 2))
for t in numpy.arange(t0, t1, 0.005):
    i = int(t * 10000); f = x[i:i + 250] * numpy.hamming(250)
    lev = 20 * numpy.log10(numpy.sqrt(numpy.mean(f ** 2)) / ref + 1e-9)
    c = numpy.correlate(f, f, "full")[249:262]
    lab = " ".join(n for n, s in ph if abs(s - t) < 0.0025)
    if c[0] <= 0: continue
    a = solve_toeplitz(c[:12], -c[1:13]); rt = numpy.roots(numpy.r_[1, a]); rt = rt[numpy.imag(rt) > 0]
    fr = numpy.angle(rt) * 10000 / (2 * numpy.pi); bw = -numpy.log(numpy.abs(rt)) * 10000 / numpy.pi
    q = sorted((a_, b_) for a_, b_ in zip(fr, bw) if a_ > 150 and b_ < 600)[:3]
    print("%.3f %5.1f dB  %s  %s" % (t, lev, "  ".join("%4.0f(%3.0f)" % fb for fb in q), lab))
