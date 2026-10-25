"""Tongue ringing: from AREALOG, the oral nodes' log areas inside vowels (score phones),
the rms of their 5-30 Hz part (difference of 6 ms and 35 ms moving averages), median over vowels.
usage: ring.py SCORE AREALOG..."""
import sys, re, numpy
phones = [(m.group(1), float(m.group(2)), float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sys.argv[1])) if m]
vowels = [(a, a + d) for n, a, d in phones if n.split("-")[0].rstrip("012") in ("AA","AE","AH","AO","AW","AY","EH","ER","EY","IH","IY","OW","OY","UH","UW","AX")]
def ma(x, k): return numpy.convolve(x, numpy.ones(k) / k, "same")
for path in sys.argv[2:]:
    rows = [l.split() for l in open(path)]
    w = min(len(r) for r in rows)
    a = numpy.array([[float(v) for v in r[:w]] for r in rows])
    t = a[:, 0]; A = numpy.log(numpy.maximum(a[:, 1 + int(7 / 0.3646):], 1e-3))
    band = numpy.array([ma(c, 6) - ma(c, 35) for c in A.T]).T
    vals = []
    for s, e in vowels:
        m = (t >= s + 0.01) & (t <= e - 0.01)
        if m.sum() > 10: vals.append(numpy.sqrt((band[m] ** 2).mean()))
    print("%-22s ringing (5-30 Hz, log area, in vowels): median %.4f  max %.4f" % (path.split("/")[-1], numpy.median(vals), numpy.max(vals)))
