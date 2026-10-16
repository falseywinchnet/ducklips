"""Tune hard vowels' organ postures to bdl's F1-F3 with Powell's method
(bounded), from the current posture and one perturbed start; cost as
organ_vopt.py. /er/ also opens the pocket under the tongue (floor 2).
usage: organ_vopt2.py IN.txt OUT.txt vowel...   (SPEAK_ENV for speak)"""
import sys, os, numpy
sys.path.insert(0, "/tmp")
import importlib.util
spec = importlib.util.spec_from_file_location("ov", "/tmp/organ_vopt.py")
from scipy.optimize import minimize
from concurrent.futures import ThreadPoolExecutor
inp, outpath, want = sys.argv[1], sys.argv[2], sys.argv[3:]
sys.argv = ["organ_vopt.py", inp, "/dev/null", "0", "aa"]
ov = importlib.util.module_from_spec(spec)
import builtins
_print = builtins.print; builtins.print = lambda *a, **k: None
spec.loader.exec_module(ov)
builtins.print = _print
KEYS = ["jaw", "tx", "ty", "tipx", "tipy", "ap", "pro", "lw", "drop"]
LOW = [0.5, 2.2, 5.5, 5.5, 6.0, 0.1, 0.0, 1.5, 0.0]
HIGH = [14.0, 5.6, 8.6, 7.8, 9.8, 2.5, 1.2, 5.0, 1.2]
floor_line = {"er": "0.0000 floor 2.0 0 0.05"}
orig_formants = ov.formants


def formants_with_floor(p, tag):
    if tag in floor_line:
        # write the score ourselves with the floor opened
        import tempfile, subprocess
        from scipy.io import wavfile
        d = tempfile.mkdtemp(prefix="ov2_", dir="/tmp/es")
        lines = ["0.0000 lungs 0 0 0.05", "0.0000 glottis 0.15 0.3 0.05", "0.0000 pitch 115 0 0.05", "0.0000 velum 0 0 0.05",
                 "0.0500 lungs 750 0 0.06", "0.0800 glottis -0.006 0.6 0.04", floor_line[tag],
                 "0.0000 body %.4f %.4f 0.05" % (p["tx"], p["ty"]), "0.0000 apex %.4f %.4f 0.05" % (p["tipx"], p["tipy"]),
                 "0.0000 mouth %.4f %.4f 0.05" % (p["ap"], p["pro"]), "0.0000 spread %.4f 0 0.05" % p["lw"],
                 "0.0000 larynx %.4f 0 0.05" % p["drop"], "0.0000 angle %.4f 0 0.05" % p["jaw"], "0.4500 mark stop 0"]
        open(d + "/s.score", "w").write("\n".join(lines) + "\n")
        subprocess.run([ov.SPEAK, d + "/s.score", d + "/s.wav"], env=ov.env, capture_output=True)
        r, x = wavfile.read(d + "/s.wav"); subprocess.run(["rm", "-rf", d])
        from scipy.signal import resample_poly, lfilter
        from scipy.linalg import solve_toeplitz
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
    return orig_formants(p, tag)


ov.formants = formants_with_floor


def tune(v):
    start = dict(ov.post[v]); best = None
    for k, jitter in enumerate([0.0, 1.0]):
        x0 = numpy.array([start[key] for key in KEYS])
        if jitter:
            x0 = x0 + numpy.array([1.5, 0.3, 0.3, 0.3, 0.3, 0.15, 0.15, 0.3, 0.2]) * numpy.random.default_rng(len(v)).uniform(-1, 1, 9)
        x0 = numpy.clip(x0, LOW, HIGH)
        def cost(x):
            p = dict(zip(KEYS, numpy.clip(x, LOW, HIGH)))
            c, f = ov.cost(v, p, start)
            return c
        res = minimize(cost, x0, method="Powell", bounds=list(zip(LOW, HIGH)), options={"maxfev": 220, "xtol": 0.02, "ftol": 0.02})
        if best is None or res.fun < best.fun: best = res
    p = dict(zip(KEYS, numpy.clip(best.x, LOW, HIGH)))
    c, f = ov.cost(v, p, start)
    return v, p, c, f


with ThreadPoolExecutor(6) as pool:
    results = list(pool.map(tune, want))
with open(outpath, "w") as out:
    for v in ov.post:
        p = ov.post[v]; note = ""
        for rv, rp, rc, rf in results:
            if rv == v: p = rp; note = "  # cost %.2f F %d %d %d" % (rc, rf[0], rf[1], rf[2])
        out.write("%s %s%s\n" % (v, " ".join("%s=%.3f" % (k, p[k]) for k in KEYS), note))
for rv, rp, rc, rf in results:
    print(rv, "cost %.2f F %s" % (rc, rf.astype(int)))
