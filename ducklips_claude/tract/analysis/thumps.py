"""Thumps: the signal low-passed at 120 Hz, its envelope (10 ms), bumps that
rise more than 6 dB over the local 100 ms median of the 120 Hz-band envelope
AND stand at least -20 dB re the file's 99th percentile broadband level.
Reports count per second, times, and the phone.
usage: thumps.py WAV [SCORE]"""
import sys, re, numpy
from scipy.io import wavfile
from scipy.signal import butter, sosfiltfilt, medfilt
r, x = wavfile.read(sys.argv[1]); x = x.astype(float)
if x.ndim > 1: x = x[:, 0]
low = sosfiltfilt(butter(4, 120, btype="low", fs=r, output="sos"), x)
w = int(0.01 * r)
env = numpy.sqrt(numpy.convolve(low ** 2, numpy.ones(w) / w, "same"))[::w]
ref = numpy.percentile(numpy.abs(x), 99)
base = medfilt(env, 11)
idx = numpy.where((env > 2.0 * base + 1e-9) & (env > ref * 0.1))[0]
ph = []
if len(sys.argv) > 2:
    ph = [(m.group(1), float(m.group(2)), float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sys.argv[2])) if m]
times = []
for i in idx:
    t = i * 0.01
    if not times or t - times[-1] > 0.03: times.append(t)
print("%s: %d thumps (%.1f/s); 0-120 Hz share of energy %.1f dB" % (sys.argv[1].split("/")[-1], len(times), len(times) / (len(x) / r), 10 * numpy.log10(numpy.sum(low ** 2) / numpy.sum(x ** 2))))
print(" ".join("%.2f%s" % (t, ("(" + "".join(n for n, s, d in ph if s <= t < s + d) + ")") if ph else "") for t in times[:40]))
