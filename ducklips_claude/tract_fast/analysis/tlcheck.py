"""Formants of one vowel posture three ways: (1) an ideal transmission line
(chain matrices, lossless sections with a small wall/viscous loss, closed
glottis, piston radiation) through Story's (2008) MRI area function;
(2) the same line through the model's area profile for the posture
(geometry_probe); (3) our engine speaking the posture (held, LPC median).
If (2) and (3) disagree, the engine's acoustics move formants and posture
tuning bends the shapes to make up for it.
usage: tlcheck.py POSTURES.txt [vowel...]"""
import sys, subprocess, tempfile, os, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, lfilter, find_peaks
from scipy.linalg import solve_toeplitz
PROBE = os.path.expanduser("~/ducklips/tract/build/geometry_probe")
SPEAK = os.path.expanduser("~/ducklips/tract/build/speak")
RHO, C = 1.14e-3, 35000.0
story = {}
for l in open(os.path.expanduser("~/ducklips/tract/data/story_2008_areas.txt")):
    f = l.split()
    if f and not f[0].startswith("#"): story[f[0]] = (float(f[1]), numpy.array([float(v) for v in f[2:]]))
def tl_formants(lengths, areas):
    fr = numpy.arange(100, 5000, 5.0); H = []
    for f in fr:
        w = 2 * numpy.pi * f; k = w / C * (1 - 0.0j) - 1j * 0.0004 * numpy.sqrt(f) / 1.0
        a = areas[-1]; rad = numpy.sqrt(a / numpy.pi)
        Zr = RHO * C / a * ((k.real * rad) ** 2 / 2 + 1j * 8 * k.real * rad / (3 * numpy.pi))
        M = numpy.eye(2, dtype=complex)
        for L, A in zip(lengths, areas):
            Z = RHO * C / A; kl = k * L
            M = M @ numpy.array([[numpy.cos(kl), 1j * Z * numpy.sin(kl)], [1j * numpy.sin(kl) / Z, numpy.cos(kl)]])
        # U_glottis = C*P_lips + D*U_lips ; P_lips = Zr U_lips
        H.append(1 / abs(M[1, 0] * Zr + M[1, 1]))
    H = 20 * numpy.log10(numpy.array(H)); pk, _ = find_peaks(H)
    return [fr[i] for i in pk[:4]]
def profile(post):
    out = subprocess.run([PROBE], input=post + "\n", capture_output=True, text=True).stdout.split("\n")
    arcs, areas = [], []
    for l in out[1:]:
        if l.startswith("cut"): f = l.split(); arcs.append(float(f[1])); areas.append(max(float(f[2]), 0.05))
    L = float(out[0].split()[1]); arcs = numpy.array(arcs)
    edges = numpy.r_[0, (arcs[1:] + arcs[:-1]) / 2, L]
    return numpy.diff(edges), numpy.array(areas)
def engine(p):
    lines = ["0.0000 lungs 0 0 0.05", "0.0000 glottis 0.15 0.3 0.05", "0.0000 pitch 115 0 0.05", "0.0000 velum 0 0 0.05",
             "0.0500 lungs 750 0 0.06", "0.0800 glottis -0.006 0.6 0.04",
             "0.0000 body %s %s 0.05" % (p["tx"], p["ty"]), "0.0000 apex %s %s 0.05" % (p["tipx"], p["tipy"]),
             "0.0000 mouth %s %s 0.05" % (p["ap"], p["pro"]), "0.0000 spread %s 0 0.05" % p["lw"],
             "0.0000 larynx %s 0 0.05" % p["drop"], "0.0000 angle %s 0 0.05" % p["jaw"], "0.4500 mark stop 0"]
    d = tempfile.mkdtemp(prefix="tl_", dir="/tmp/es"); open(d + "/s.score", "w").write("\n".join(lines) + "\n")
    subprocess.run([SPEAK, d + "/s.score", d + "/s.wav"], env=dict(os.environ, ROOM="0"), capture_output=True)
    r, x = wavfile.read(d + "/s.wav"); subprocess.run(["rm", "-rf", d])
    x = lfilter([1, -0.63], [1], resample_poly(x.astype(float), 10000, r))[1800:4200]; out = []
    for i in range(0, len(x) - 250, 100):
        f = x[i:i + 250] * numpy.hamming(250); c = numpy.correlate(f, f, "full")[249:262]
        if c[0] <= 0: continue
        a = solve_toeplitz(c[:12], -c[1:13]); rt = numpy.roots(numpy.r_[1, a]); rt = rt[numpy.imag(rt) > 0]
        fr = numpy.angle(rt) * 10000 / (2 * numpy.pi); bw = -numpy.log(numpy.abs(rt)) * 10000 / numpy.pi
        q = sorted(a_ for a_, b_ in zip(fr, bw) if a_ > 150 and b_ < 600)[:3]
        if len(q) == 3: out.append(q)
    return numpy.median(out, axis=0)
post = {}
for l in open(sys.argv[1]):
    f = l.split("#")[0].split()
    if f and not l.startswith("#") and "=" in l and not f[0].startswith(("rule", "anat")):
        post[f[0].lower()] = " ".join(f[1:])
want = sys.argv[2:] or sorted(post)
print("vowel  Story-TL F1-F3     model-TL F1-F3     engine F1-F3")
for v in want:
    if v not in post: continue
    p = dict(kv.split("=") for kv in post[v].split())
    dl, A = story.get(v, (None, None))
    s = tl_formants(numpy.full(len(A), dl), A)[:3] if A is not None else [0, 0, 0]
    m = tl_formants(*profile(post[v]))[:3]
    e = engine(p)
    print("%-4s  %5.0f %5.0f %5.0f   %5.0f %5.0f %5.0f   %5.0f %5.0f %5.0f" % (v, *s, *m, *e))
