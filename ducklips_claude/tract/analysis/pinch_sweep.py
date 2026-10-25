"""Sweep pinch-tube settings on a held /ah/: per variant the mean flow, the
pinch's area range, its oscillation (autocorrelation period 2.5-20 ms and
strength) from WAVEDUMP at 1.00-1.10 s.
usage: pinch_sweep.py BASE.params SCORE key=v1,v2 [key=...]  (pinch.* keys)"""
import itertools, os, subprocess, sys, numpy
base, score = sys.argv[1], sys.argv[2]
axes = [(a.split("=")[0], a.split("=")[1].split(",")) for a in sys.argv[3:]]
speak = os.path.expanduser("~/ducklips/tract/build/speak")
for combo in itertools.product(*[values for _, values in axes]):
    lines = open(base).read().rstrip("\n") + "\npinch {\n  on 1\n" + "".join("  %s %s\n" % (k.split(".")[1], v) for (k, _), v in zip(axes, combo)) + "}\n"
    path = "/tmp/es/pinch/sweep.params"
    open(path, "w").write(lines)
    env = dict(os.environ, FOLDTUBE="1", PARAMS=path, WAVEDUMP="1.00,1.10")
    out = subprocess.run([speak, score, "/dev/null", "73021"], env=env, capture_output=True, text=True).stdout
    rows = [l.split("cord ")[1].split() for l in out.splitlines() if l.startswith("wave") and "cord " in l]
    t = numpy.array([float(l.split()[1]) for l in out.splitlines() if l.startswith("wave")])
    if len(rows) < 100:
        print(combo, "failed"); continue
    u = numpy.array([float(r[0]) for r in rows]); a = numpy.array([float(r[1]) for r in rows])
    x = a - a.mean(); ac = numpy.correlate(x, x, "full")[len(x) - 1:]
    rate = 1 / numpy.median(numpy.diff(t)); lags = numpy.arange(len(ac)) / rate
    m = (lags > 0.0025) & (lags < 0.02)
    k = numpy.argmax(ac[m]); strength = ac[m][k] / max(ac[0], 1e-12)
    print("%-40s flow %5.0f (%4.0f-%4.0f)  area %.4f-%.4f  period %.2f ms (%3.0f Hz) strength %.2f" % (
        " ".join("%s=%s" % (k_.split(".")[1], v) for (k_, _), v in zip(axes, combo)), u.mean(), u.min(), u.max(), a.min(), a.max(), 1000 * lags[m][k], 1 / lags[m][k], strength))
