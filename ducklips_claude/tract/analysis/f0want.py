import sys, re, numpy, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
out, score = sys.argv[1], sys.argv[2]
ph = [(m.group(1), float(m.group(2)), float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(score)) if m]
fig, ax = plt.subplots(len(sys.argv) - 3, 1, figsize=(16, 3.2 * (len(sys.argv) - 3)), squeeze=False)
for k, p in enumerate(sys.argv[3:]):
    c = numpy.loadtxt(p); c = c[(c[:, 1] > 2.5) & (c[:, 1] < 15)]
    a = ax[k, 0]; a.plot(c[:, 0], 1000 / c[:, 1], ".", ms=4, label="F0"); a.plot(c[:, 0], c[:, 4], "-", lw=1, label="intended")
    a.set_ylim(80, 180); a.grid(); a.set_title(p)
    for n, s, d in ph:
        a.axvline(s, color="gray", lw=0.3); a.text(s + 0.003, 172, n.rstrip("-0123456789"), fontsize=8)
plt.tight_layout(); plt.savefig(out, dpi=70)
