"""Share of voiced cycles (audio waveform matching, stem_wm.txt) per phone
in which the glottis shut (the CYCLES log's smallest area per cycle under
0.001 cm^2, the chink being 0.0002). Cycles whose area never fell below
0.01 have no CYCLES line and count as open.
usage: closure.py SCORE..."""
import sys, re, numpy, collections
tot = collections.Counter(); hit = collections.Counter(); mins = collections.defaultdict(list)
for sc in sys.argv[1:]:
    stem = sc[:-6]
    wm = numpy.loadtxt(stem + "_wm.txt", ndmin=2); cy = numpy.loadtxt(stem + ".cyc", ndmin=2)
    ph = [(m.group(1).rstrip("-0123456789"), float(m.group(2)), float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sc)) if m]
    for t, per, c in wm:
        if c < 0.8: continue
        mid = t + per / 2000
        for n, s, d in ph:
            if s <= mid < s + d:
                tot[n] += 1
                near = cy[numpy.abs(cy[:, 0] - mid) < per / 1000] if len(cy) else cy
                if len(near) and near[:, 7].min() < 0.001: hit[n] += 1
                if len(near): mins[n].append(near[:, 7].min())
                else: mins[n].append(0.05)
for n in sorted(tot, key=lambda k: -tot[k]):
    if tot[n] >= 10: print("%-4s cycles %4d shut %3.0f%% min area median %.4f" % (n, tot[n], 100 * hit[n] / tot[n], numpy.median(mins[n])))
