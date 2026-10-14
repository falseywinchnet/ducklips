"""Random search over fold materials for a firm, stable source.
Scores: OQ near 0.5, jitter small, F0 near target, on two vowels x two pitches."""
import os, random, subprocess, sys, numpy
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
names = ["F_LM", "F_UM", "F_BM", "F_LK", "F_UK", "F_SK", "F_BK", "F_LD", "F_UD", "F_BD", "F_CONTACT", "THICK1", "THICK2"]
base = [0.153, 0.122, 0.898, 179621, 1013, 20457, 284701, 0.171, 0.282, 0.227, 15.94, 1.03, 0.0275]
cases = [("7.75", "0.582", "0.71", "105", "0.0", 0.5), ("4.0", "1.35", "6.75", "145", "0.0", 0.5), ("4.0", "1.35", "6.75", "120", "-0.02", 0.35), ("7.75", "0.582", "0.71", "120", "0.03", 0.7)]

def measure(env):
    total = 0.0
    detail = []
    for place, area, lips, pitch, gap, want in cases:
        pressure = str(750 + 5 * (float(pitch) - 125))
        e = dict(os.environ, **env, TRACE="/tmp/fs_tr.txt", ASP="0")
        subprocess.run(["./build/vowel_probe", "/tmp/fs.wav", place, area, lips, pitch, pressure, gap, "0.8"], env=e, stdout=subprocess.DEVNULL)
        d = numpy.loadtxt("/tmp/fs_tr.txt")[int(0.45 * 48000):]
        a = d[:, 0]
        lo, hi = a.min(), a.max()
        if hi - lo < 0.02:
            total += 5.0
            detail.append("dead")
            continue
        level = lo + 0.3 * (hi - lo)
        rising = numpy.where((a[:-1] < level) & (a[1:] >= level))[0]
        if len(rising) < 5:
            total += 5.0
            detail.append("dead")
            continue
        periods = numpy.diff(rising) / 48000.0
        f0 = 1 / periods.mean()
        jitter = numpy.mean(numpy.abs(numpy.diff(periods))) / periods.mean()
        oq = (a > lo + 0.05 * (hi - lo)).mean()
        total += 4 * (oq - want) ** 2 + 20 * jitter + 2 * numpy.log(f0 / float(pitch)) ** 2
        detail.append("%.0f/%.2f/%.1f%%" % (f0, oq, jitter * 100))
    return total, " ".join(detail)

best = (measure({})[0], base)
print("base", best[0], measure({})[1], flush=True)
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 60):
    scale = 0.5 if trial < 30 else 0.25
    values = [v * numpy.exp(random.gauss(0, scale)) for v in best[1]]
    env = {n: "%g" % v for n, v in zip(names, values)}
    score, info = measure(env)
    if score < best[0]:
        best = (score, values)
        print("%3d %.3f %s | %s" % (trial, score, info, " ".join("%s=%g" % (n, v) for n, v in zip(names, values))), flush=True)
