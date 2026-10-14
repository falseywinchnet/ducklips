"""Per-vowel voicing of renders of one score: share of 10 ms frames with a
period (normalised autocorrelation >= 0.5), and median F0, side by side.
usage: pervowel.py SCORE WAV [WAV...]"""
import sys, re, numpy
sys.path.insert(0, "/tmp")
from vw import load
vowels = [(m.group(1), float(m.group(2)), float(m.group(3))) for m in
          (re.match(r"# ([A-Z]+\d) ([\d.]+) ([\d.]+) ([\d.]+)", l) for l in open(sys.argv[1])) if m and float(m.group(4)) > 0]
xs = [load(w) for w in sys.argv[2:]]
def frame(x, t):
    a = x[int(t * 16000):int((t + 0.03) * 16000)]; a = a - a.mean()
    c = numpy.correlate(a, a, "full")[len(a) - 1:]; lo, hi = 64, 267
    if c[0] <= 0: return 0.0
    k = lo + numpy.argmax(c[lo:hi]); return 16000 / k if c[k] >= 0.5 * c[0] else 0.0
for name, s, d in vowels:
    cells = []
    for x in xs:
        f = [frame(x, t) for t in numpy.arange(s, s + d - 0.02, 0.01)]
        v = [q for q in f if q > 0]
        cells.append("%3d%% %5.0f" % (100 * len(v) / max(len(f), 1), numpy.median(v) if v else 0))
    print("%-5s %.3f  %s" % (name, s, "   ".join(cells)))
