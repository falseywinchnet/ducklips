import sys, numpy, re
d = numpy.loadtxt(sys.argv[1]); ph = []
for l in open(sys.argv[2]):
    m = re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l)
    if m: ph.append((m.group(1), float(m.group(2)), float(m.group(3))))
def at(t):
    for n, s, du in ph:
        if s <= t < s + du: return "%s+%.0fms" % (n, 1000 * (t - s))
    return "-"
p = d[:, 1]
for i in range(1, len(d) - 1):
    if all(2.5 < x < 15 for x in p[i - 1:i + 2]):
        m = (p[i - 1] + p[i] + p[i + 1]) / 3; r = abs(p[i] - m) / m * 100
        if r > 3.33: print("%.4f %5.1f Hz (nbrs %5.1f %5.1f) rap %.1f%% %s lung %.0f flow %.0f" % (d[i, 0], 1000 / p[i], 1000 / p[i - 1], 1000 / p[i + 1], r, at(d[i, 0]), d[i, 6], d[i, 2]))
