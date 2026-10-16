"""Find clicks: samples whose high-passed (3 kHz) magnitude exceeds 8x the
local (20 ms) rms of the high-passed signal; report time, size re the file's
speech level, and the phone there.
usage: clicks.py WAV [SCORE]"""
import sys, re, numpy
from scipy.io import wavfile
from scipy.signal import butter, sosfiltfilt
r, x = wavfile.read(sys.argv[1]); x = x.astype(float)
if x.ndim > 1: x = x[:, 0]
ph = []
if len(sys.argv) > 2:
    ph = [(m.group(1), float(m.group(2)), float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sys.argv[2])) if m]
peak = numpy.abs(x).max(); speech = numpy.percentile(numpy.abs(x), 99)
print("peak / 99th percentile %.1f dB" % (20 * numpy.log10(peak / speech)))
order = numpy.argsort(-numpy.abs(x))
seen = []
for i in order[:4000]:
    t = i / r
    if any(abs(t - s) < 0.01 for s in seen): continue
    seen.append(t)
    lab = " ".join(n for n, s, d in ph if s <= t < s + d)
    print("%.4f s  %.1f dB re 99th  %s" % (t, 20 * numpy.log10(abs(x[i]) / speech), lab))
    if len(seen) >= 12: break
