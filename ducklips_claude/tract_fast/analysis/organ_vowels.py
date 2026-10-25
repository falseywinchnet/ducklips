"""Write a score holding vowels through the articulatory body (organ parts),
from a geofit result, and report each vowel's formants (LPC, middle 50%) of
the rendered sound against Story's tube and bdl's averages.
usage: organ_vowels.py FIT.txt OUT_PREFIX [vowel...]   (SPEAK_ENV for speak)"""
import sys, os, subprocess, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, lfilter
from scipy.linalg import solve_toeplitz
fit, prefix = sys.argv[1], sys.argv[2]
want = sys.argv[3:] or ["aa", "iy", "uw", "eh", "ow", "ae", "ih", "ah"]
post = {}
for l in open(fit):
    f = l.split("#")[0].split()
    if len(f) > 1 and "=" in f[1] and not f[0].startswith(("rule", "anat")):
        post[f[0]] = dict(x.split("=") for x in f[1:])
lines = ["0.0000 lungs 0 0 0.05", "0.0000 glottis 0.15 0.3 0.05", "0.0000 pitch 115 0 0.05", "0.0000 velum 0 0 0.05",
         "0.1000 lungs 750 0 0.08", "0.1500 glottis -0.006 0.6 0.05"]
t = 0.0; spans = []
for v in want:
    p = {k: float(x) for k, x in post[v].items()}
    move = 0.12 if t > 0 else 0.05
    lines += ["%.4f body %.3f %.3f %.3f" % (t, p["tx"], p["ty"], move), "%.4f apex %.3f %.3f %.3f" % (t, p["tipx"], p["tipy"], move),
              "%.4f mouth %.3f %.3f %.3f" % (t, p["ap"], p["pro"], move), "%.4f spread %.3f 0 %.3f" % (t, p["lw"], move),
              "%.4f larynx %.3f 0 %.3f" % (t, p["drop"], move), "%.4f angle %.3f 0 %.3f" % (t, p["jaw"], move)]
    start = max(t, 0.25); spans.append((v, start + 0.12, t + 0.5)); t += 0.5
lines += ["%.4f lungs 0 0 0.08" % t, "%.4f glottis 0.3 0.3 0.05" % t, "%.4f mark stop 0" % (t + 0.2)]
open(prefix + ".score", "w").write("\n".join(lines) + "\n")
env = dict(os.environ, ROOM="0"); env.update(dict(x.split("=") for x in os.environ.get("SPEAK_ENV", "").split() if "=" in x))
print(subprocess.run(["/Users/joshuahkuttenkuler/ducklips/tract/build/speak", prefix + ".score", prefix + ".wav"], env=env, capture_output=True, text=True).stderr.strip())
r, x = wavfile.read(prefix + ".wav"); x = lfilter([1, -0.63], [1], resample_poly(x.astype(float), 10000, r))
BDL = {"iy": (294, 2213), "ih": (432, 1780), "eh": (557, 1609), "ae": (661, 1640), "aa": (708, 1139), "ah": (590, 1240),
       "ao": (620, 990), "ow": (510, 1060), "uh": (450, 1300), "uw": (340, 1239), "ey": (430, 1950)}
for v, s, e in spans:
    seg = x[int(s * 10000):int(e * 10000)]; out = []
    for i in range(0, len(seg) - 250, 100):
        f = seg[i:i + 250] * numpy.hamming(250); c = numpy.correlate(f, f, "full")[249:262]
        a = solve_toeplitz(c[:12], -c[1:13]); rt = numpy.roots(numpy.r_[1, a]); rt = rt[numpy.imag(rt) > 0]
        fr = numpy.angle(rt) * 10000 / (2 * numpy.pi); bw = -numpy.log(numpy.abs(rt)) * 10000 / numpy.pi
        q = sorted(a_ for a_, b_ in zip(fr, bw) if a_ > 150 and b_ < 600)[:3]
        if len(q) == 3: out.append(q)
    m = numpy.median(out, axis=0) if out else [0, 0, 0]
    b = BDL.get(v, (0, 0))
    print("%-3s F1 %4.0f F2 %4.0f F3 %4.0f   bdl F1 %4.0f F2 %4.0f  dist %.1f st" % (v, m[0], m[1], m[2], b[0], b[1],
          numpy.hypot(12 * numpy.log2(m[0] / b[0]), 12 * numpy.log2(m[1] / b[1])) if b[0] else 0))
