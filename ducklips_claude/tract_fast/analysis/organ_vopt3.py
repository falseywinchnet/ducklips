"""Tune vowels' organ postures to bdl's F1-F3 while keeping their shapes
near Story's (2008) MRI area functions: cost = formant distance (semitones,
F1-F3) + W x rms log-area error against Story at 16 relative positions along
the tract (glottis to lips). Powell's method (bounded), from START's posture.
(organ_vopt.py's cost, 0.5 per cm the tongue moved, let the tuner flatten
every vowel toward a uniform tube: /aa/'s pharynx 0.6-0.8 cm^2 against
Story's 0.25, its mouth 5.6 against 9.3.)
usage: organ_vopt3.py START.txt CURRENT.txt OUT.txt W vowel...
   START: postures to start from (geofit); CURRENT: the file whose other
   lines (vowels not tuned, comments) are kept.  SPEAK_ENV for speak;
   SPEAK, the speak binary (default ~/ducklips/tract/build/speak)."""
import sys, os, subprocess, tempfile, shutil, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, lfilter
from scipy.linalg import solve_toeplitz
from scipy.optimize import minimize
from concurrent.futures import ThreadPoolExecutor
start_path, current_path, outpath, W, want = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5:]
SPEAK = os.path.expanduser(os.environ.get("SPEAK", "~/ducklips/tract/build/speak"))
ENV = dict(os.environ, ROOM="0")
ENV.update(dict(x.split("=", 1) for x in os.environ.get("SPEAK_ENV", "").split() if "=" in x))
# bdl's running-speech medians (LPC, vowel middles, a0001-a0060; analysis/f123.py).
BDL = {"iy": (319, 2333, 2822), "ih": (433, 2046, 2746), "eh": (598, 1673, 2641), "ae": (638, 1833, 2619),
       "aa": (724, 1221, 2585), "ah": (578, 1231, 2651), "ao": (460, 918, 2549), "ow": (467, 974, 2623),
       "uh": (429, 1223, 2717), "uw": (337, 1238, 2496), "ey": (485, 2053, 2666),
       "ax": (549, 1449, 2660)}  # bdl's schwa, his most frequent phone (158 in a0001-a0060)
# TARGETS=file.json: {vowel: [F1, F2, F3]} in place of bdl's for those vowels
# (2026-10-18: the NATO reader's stressed vowels, scaled by 1 / 1.12 for his
# shorter tract, /tmp/es/readervowels.py; the owner chose the NATO voice as
# the reference).
if os.environ.get("TARGETS"):
    import json
    for vowel, values in json.load(open(os.environ["TARGETS"])).items():
        BDL[vowel] = tuple(values)


def formants(p):
    """The engine speaking the posture held: median LPC F1-F3."""
    lines = ["0.0000 lungs 0 0 0.05", "0.0000 glottis 0.15 0.3 0.05", "0.0000 pitch %s 0 0.05" % os.environ.get("VOWEL_PITCH", "115"), "0.0000 velum 0 0 0.05",
             "0.0500 lungs 750 0 0.06", "0.0800 glottis -0.006 0.6 0.04",
             "0.0000 body %.4f %.4f 0.05" % (p["tx"], p["ty"]), "0.0000 apex %.4f %.4f 0.05" % (p["tipx"], p["tipy"]),
             "0.0000 mouth %.4f %.4f 0.05" % (p["ap"], p["pro"]), "0.0000 spread %.4f 0 0.05" % p["lw"],
             "0.0000 larynx %.4f 0 0.05" % p["drop"], "0.0000 angle %.4f 0 0.05" % p["jaw"], "0.4500 mark stop 0"]
    d = tempfile.mkdtemp(prefix="v3_", dir="/tmp/es")
    try:
        open(d + "/s.score", "w").write("\n".join(lines) + "\n")
        subprocess.run([SPEAK, d + "/s.score", d + "/s.wav"], env=ENV, capture_output=True, timeout=60)
        r, x = wavfile.read(d + "/s.wav")
    finally:
        shutil.rmtree(d, ignore_errors=True)
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


PROBE = os.path.expanduser("~/ducklips/tract/build/geometry_probe")
KEYS = ["jaw", "tx", "ty", "tipx", "tipy", "ap", "pro", "lw", "drop"]
LOW = [0.5, 2.2, 5.5, 5.5, 6.0, 0.1, 0.0, 1.5, 0.0]
HIGH = [20.0, 5.6, 8.6, 7.8, 9.8, 2.5, 1.2, 5.0, 1.2]  # (jaw 14 until 2026-10-18)
story = {}
for l in open(os.path.expanduser("~/ducklips/tract/data/story_2008_areas.txt")):
    f = l.split()
    if f and not f[0].startswith("#"): story[f[0]] = numpy.array([float(v) for v in f[2:]])
POS = (numpy.arange(16) + 0.5) / 16
MAXFEV = int(os.environ.get("MAXFEV", "260"))


def read(path):
    out = {}
    for l in open(path):
        f = l.split("#")[0].split()
        if len(f) > 1 and "=" in f[1] and not f[0].startswith(("rule", "anat")):
            out[f[0].lower()] = {k: float(x) for k, x in (y.split("=") for y in f[1:])}
    return out


def shape_error(v, p):
    if v not in story:
        return 0.0  # (no MRI shape for it: formants alone)
    line = " ".join("%s=%.5f" % (k, p[k]) for k in KEYS) + " " + os.environ.get("PROBE_EXTRA", "")  # e.g. anat.vtl=1
    lines = subprocess.run([PROBE], input=line + "\n", capture_output=True, text=True).stdout.split("\n")
    length = float(lines[0].split()[1]); arcs, areas = [], []
    for l in lines[1:]:
        if l.startswith("cut"):
            f = l.split(); arcs.append(float(f[1])); areas.append(float(f[2]))
    A = story[v]; target = A[(POS * len(A)).astype(int)]
    model = numpy.interp(POS * length, arcs, areas)
    return float(numpy.sqrt(numpy.mean(numpy.log((model + 0.1) / (target + 0.1)) ** 2)))


def cost(v, x):
    p = dict(zip(KEYS, numpy.clip(x, LOW, HIGH)))
    f = formants(p); b = BDL[v]
    err = float(numpy.sqrt(sum((12 * numpy.log2(f[i] / b[i])) ** 2 for i in range(3))))
    shape = shape_error(v, p)
    return err + W * shape, err, shape, f


# Anatomy the formants alone do not enforce: in a vowel the tip rests behind
# the lower incisors (lingual surface about (7.5, 7.3), floor (7.1, 6.6) in
# the jaw's frame), and unrounded vowels do not push the lips forward. (Free,
# the tuner lifted the tip toward the ridge in AH 9.0, IY 8.4, IH 8.0 and
# protruded AA's lips 0.69 cm: the owner heard "a bad lisp", a speech
# impediment.)
UNROUNDED = {"iy", "ih", "ey", "eh", "ae", "ah", "aa", "ax"}
# The jaw opens for open vowels (degrees). Fitted to bdl's formants with
# the jaw free, /eh/ came out at 1.05, /ae/ 6.5 and /iy/ 1.9: the tongue
# bunched under a nearly shut jaw filled the mouth, and the owner heard
# alfa, golf, mike, quebec and yankee "misshapen" (2026-10-18; the NATO
# reader's /ae/ F1 about 1000 Hz against our 650). The owner, on refitting
# with these minimums: "yes, do that."
JAW_MIN = {"aa": 14.0, "ae": 13.0, "ao": 11.0, "ah": 9.0, "eh": 8.0, "ey": 6.0, "ow": 6.0, "ih": 5.0, "uh": 5.0,
           "ax": 4.0, "iy": 3.0, "uw": 2.0}


def bounds_for(v):
    low = list(LOW); high = list(HIGH)
    low[KEYS.index("tipx")], high[KEYS.index("tipx")] = 6.7, 7.4
    low[KEYS.index("tipy")], high[KEYS.index("tipy")] = 6.6, 7.4
    if v in JAW_MIN:
        low[KEYS.index("jaw")] = JAW_MIN[v]
    if v in UNROUNDED:
        high[KEYS.index("pro")] = 0.15
    if v in ("uw", "ow"):
        # The owner (2026-10-17): "ooo does not need to pull the tongue
        # back." The lips round it (protruded, a small opening, the jaw
        # mostly closed); the body stays at or forward of neutral. (Free,
        # the tuner held /ow/ with the body back at x 2.75 and the jaw 12.7
        # degrees open.)
        low[KEYS.index("tx")] = 3.6
        high[KEYS.index("jaw")] = 8.0
        low[KEYS.index("pro")] = 0.6
        high[KEYS.index("ap")] = 0.6
    if v == "ax":
        # A relaxed schwa: jaw a little open, lips about 1 cm apart. (Free,
        # the formants were met with the jaw shut and the lips 2.4 cm open.)
        low[KEYS.index("jaw")], high[KEYS.index("jaw")] = max(3.0, JAW_MIN["ax"]), 8.0
        low[KEYS.index("ap")], high[KEYS.index("ap")] = 0.6, 1.3
    return low, high


def tune(v):
    low, high = bounds_for(v)
    x0 = numpy.clip(numpy.array([starts[v][k] for k in KEYS]), low, high)
    c0 = cost(v, x0)
    best = {"c": c0[0], "x": x0, "info": c0}

    def objective(x):
        x = numpy.clip(x, low, high)
        c = cost(v, x)
        if c[0] < best["c"]:
            best.update(c=c[0], x=x, info=c)
        return c[0]
    minimize(objective, x0, method="Powell", bounds=list(zip(low, high)), options={"maxfev": MAXFEV, "xtol": 0.02, "ftol": 0.01})
    return v, c0, best


starts = read(start_path)
with ThreadPoolExecutor(8) as pool:
    results = list(pool.map(tune, want))
tuned = {v: b for v, c0, b in results}
with open(outpath, "w") as out:
    for l in open(current_path):
        f = l.split("#")[0].split()
        name = f[0].lower() if f else ""
        if name in tuned:
            b = tuned[name]; c, err, shape, fm = b["info"]
            out.write("%s %s  # cost %.2f (formants %.2f st, shape rms %.3f) F %d %d %d\n" % (
                f[0], " ".join("%s=%.3f" % (k, x) for k, x in zip(KEYS, b["x"])), c, err, shape, fm[0], fm[1], fm[2]))
        else:
            out.write(l)
for v, c0, b in results:
    print("%-3s start cost %.2f (%.2f st, shape %.3f)  ->  %.2f (%.2f st, shape %.3f) F %s" % (
        v, c0[0], c0[1], c0[2], b["info"][0], b["info"][1], b["info"][2], b["info"][3].astype(int)))
