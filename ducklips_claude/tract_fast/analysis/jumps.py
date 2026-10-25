"""Pitch jumps at closures and releases: from a CYCLES log, the largest
cycle-to-cycle F0 change (%) within 20 ms of each listed boundary, and the
median cycle-to-cycle change elsewhere. usage: jumps.py CYCLES t1 t2 ..."""
import sys, numpy
d = numpy.loadtxt(sys.argv[1], ndmin=2); t = d[:, 0]; T = d[:, 1]
ok = (T > 3) & (T < 15); t = t[ok]; f = 1000.0 / T[ok]
ch = numpy.abs(numpy.diff(f)) / f[:-1] * 100; tc = t[1:]
cont = numpy.diff(t) < 0.015
res = []
for b in [float(v) for v in sys.argv[2:]]:
    near = cont & (numpy.abs(tc - b) < 0.02)
    res.append(ch[near].max() if near.any() else float("nan"))
far = cont.copy()
for b in [float(v) for v in sys.argv[2:]]:
    far &= numpy.abs(tc - b) >= 0.02
print("boundary jumps %s %%   elsewhere median %.2f %%  F0 median %.0f" % (" ".join("%.1f" % v for v in res), numpy.median(ch[far]), numpy.median(f)))
