"""Tune each vowel's organ posture to bdl's F1-F3 by coordinate search,
starting from the MRI fit and penalised for leaving it (0.5 st per cm of
tongue movement). Each trial renders a held vowel (0.4 s) through speak with
the articulatory body and measures LPC formants over its middle.
usage: organ_vopt.py FIT.txt OUT.txt ROUNDS [vowel...]   (SPEAK_ENV for speak)"""
import sys, os, subprocess, numpy, tempfile
from concurrent.futures import ThreadPoolExecutor
from scipy.io import wavfile
from scipy.signal import resample_poly, lfilter
from scipy.linalg import solve_toeplitz
SPEAK = "/Users/joshuahkuttenkuler/ducklips/tract/build/speak"
# bdl's running-speech medians (LPC, vowel middles, a0001-a0060; analysis/f123.py).
BDL = {"iy": (319, 2333, 2822), "ih": (433, 2046, 2746), "eh": (598, 1673, 2641), "ae": (638, 1833, 2619),
       "aa": (724, 1221, 2585), "ah": (578, 1231, 2651), "ao": (460, 918, 2549), "ow": (467, 974, 2623),
       "uh": (429, 1223, 2717), "uw": (337, 1238, 2496), "ey": (485, 2053, 2666),
       "er": (419, 1438, 1853)}
STEPS = {"tx": 0.25, "ty": 0.2, "tipx": 0.2, "tipy": 0.2, "ap": 0.15, "pro": 0.15, "lw": 0.3, "jaw": 1.5, "drop": 0.2}
fit, outpath, rounds = sys.argv[1], sys.argv[2], int(sys.argv[3])
post = {}
head = []
for l in open(fit):
    f = l.split("#")[0].split()
    if len(f) > 1 and "=" in f[1] and not f[0].startswith(("rule", "anat")):
        post[f[0]] = {k: float(x) for k, x in (y.split("=") for y in f[1:])}
    elif f:
        head.append(l.rstrip())
want = sys.argv[4:] or [v for v in BDL if v in post]
env = dict(os.environ, ROOM="0"); env.update(dict(x.split("=") for x in os.environ.get("SPEAK_ENV", "").split() if "=" in x))


def formants(p, tag):
    lines = ["0.0000 lungs 0 0 0.05", "0.0000 glottis 0.15 0.3 0.05", "0.0000 pitch 115 0 0.05", "0.0000 velum 0 0 0.05",
             "0.0500 lungs 750 0 0.06", "0.0800 glottis -0.006 0.6 0.04",
             "0.0000 body %.4f %.4f 0.05" % (p["tx"], p["ty"]), "0.0000 apex %.4f %.4f 0.05" % (p["tipx"], p["tipy"]),
             "0.0000 mouth %.4f %.4f 0.05" % (p["ap"], p["pro"]), "0.0000 spread %.4f 0 0.05" % p["lw"],
             "0.0000 larynx %.4f 0 0.05" % p["drop"], "0.0000 angle %.4f 0 0.05" % p["jaw"], "0.4500 mark stop 0"]
    d = tempfile.mkdtemp(prefix="ov_", dir="/tmp/es")
    open(d + "/s.score", "w").write("\n".join(lines) + "\n")
    subprocess.run([SPEAK, d + "/s.score", d + "/s.wav"], env=env, capture_output=True)
    r, x = wavfile.read(d + "/s.wav"); subprocess.run(["rm", "-rf", d])
    x = lfilter([1, -0.63], [1], resample_poly(x.astype(float), 10000, r))[1800:4200]
    out = []
    for i in range(0, len(x) - 250, 100):
        f = x[i:i + 250] * numpy.hamming(250); c = numpy.correlate(f, f, "full")[249:262]
        if c[0] <= 0: continue
        a = solve_toeplitz(c[:12], -c[1:13]); rt = numpy.roots(numpy.r_[1, a]); rt = rt[numpy.imag(rt) > 0]
        fr = numpy.angle(rt) * 10000 / (2 * numpy.pi); bw = -numpy.log(numpy.abs(rt)) * 10000 / numpy.pi
        q = sorted(a_ for a_, b_ in zip(fr, bw) if a_ > 150 and b_ < 600)[:3]
        if len(q) == 3: out.append(q)
    return numpy.median(out, axis=0) if out else numpy.array([1.0, 1.0, 1.0])


def cost(v, p, start):
    f = formants(p, v); b = BDL[v]
    err = numpy.sqrt((12 * numpy.log2(f[0] / b[0])) ** 2 + (12 * numpy.log2(f[1] / b[1])) ** 2 + (12 * numpy.log2(f[2] / b[2])) ** 2)
    move = numpy.hypot(p["tx"] - start["tx"], p["ty"] - start["ty"])
    return err + 0.5 * move, f


def tune(v):
    start = dict(post[v]); best = dict(start); c, f = cost(v, best, start); log = ["%s start %.2f F %s" % (v, c, f.astype(int))]
    for r in range(rounds):
        for k, step in STEPS.items():
            for sign in (1, -1):
                trial = dict(best); trial[k] = trial[k] + sign * step
                if k in ("ap", "pro", "lw", "drop", "jaw") and trial[k] < 0: continue
                tc, tf = cost(v, trial, start)
                if tc < c - 0.05: best, c, f = trial, tc, tf
        log.append("%s round %d %.2f F %s" % (v, r, c, f.astype(int)))
    return v, best, c, f, log


with ThreadPoolExecutor(8) as pool:
    results = list(pool.map(tune, want))
with open(outpath, "w") as out:
    for h in head:
        out.write(h + "\n")
    for v in post:
        p = post[v]; note = ""
        for rv, rb, rc, rf, log in results:
            if rv == v: p = rb; note = "  # cost %.2f F %d %d %d" % (rc, rf[0], rf[1], rf[2])
        out.write("%s %s%s\n" % (v, " ".join("%s=%.3f" % (k, p[k]) for k in ["jaw", "tx", "ty", "tipx", "tipy", "ap", "pro", "lw", "drop"]), note))
for rv, rb, rc, rf, log in results:
    print("\n".join(log))
