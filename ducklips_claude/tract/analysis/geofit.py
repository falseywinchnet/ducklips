"""Fit the midsagittal model to Story's (2008) MRI area functions: one set of
area rules (alpha, beta per region) shared by all vowels, and one posture
per vowel. Residuals: log area at Story's section centres, and the tract
length. Writes the fitted rules and postures.
usage: geofit.py OUT.txt"""
import sys, subprocess, numpy
from scipy.optimize import least_squares

PROBE = "/Users/joshuahkuttenkuler/ducklips/tract/build/geometry_probe"
story = {}
for l in open("/Users/joshuahkuttenkuler/ducklips/tract/data/story_2008_areas.txt"):
    f = l.split()
    if f and not f[0].startswith("#"):
        story[f[0]] = (float(f[1]), numpy.array([float(v) for v in f[2:]]))
VOWELS = sorted(story)
REGIONS = ["larynx", "pharynx", "velar", "palatal", "alveolar", "teeth"]
NSHARED = 2 * len(REGIONS) + 2   # rules, then anat.scale and anat.root
KEYS = ["jaw", "tx", "ty", "tipx", "tipy", "ap", "pro", "lw", "drop"]
LOW = [0.0, 2.0, 4.5, 5.0, 6.0, 0.05, 0.0, 1.5, 0.0]
HIGH = [20.0, 6.0, 8.2, 8.0, 9.0, 3.0, 1.5, 5.0, 1.5]
START = {"iy": [4, 4.9, 7.0, 6.9, 7.4, 1.0, 0, 4.0, 0], "ih": [5, 4.6, 6.9, 6.9, 7.3, 1.0, 0, 4.0, 0],
         "ey": [6, 4.6, 6.7, 6.9, 7.2, 1.2, 0, 4.0, 0], "eh": [8, 4.3, 6.5, 6.9, 7.1, 1.4, 0, 4.0, 0],
         "ae": [12, 4.0, 6.0, 6.9, 7.0, 1.8, 0, 4.2, 0], "ah": [8, 3.6, 6.2, 6.9, 7.0, 1.2, 0, 3.6, 0],
         "aa": [12, 3.0, 5.6, 6.9, 7.0, 1.8, 0, 3.6, 0], "ao": [10, 3.0, 5.9, 6.9, 7.0, 1.2, 0.4, 3.0, 0.3],
         "ow": [7, 3.1, 6.4, 6.9, 7.0, 0.8, 0.7, 2.6, 0.6], "uh": [6, 3.3, 6.8, 6.9, 7.0, 0.8, 0.6, 2.6, 0.8],
         "uw": [5, 3.2, 7.0, 6.9, 7.0, 0.3, 1.0, 2.0, 1.0]}
RULE_START = [1.4, 1.4, 0.84, 1.4, 0.6, 1.6, 0.53, 2.0, 0.62, 1.5, 2.7, 1.0, 0.97, 0.5]
probe = subprocess.Popen([PROBE], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)


def shape(posture, rules):
    keys = " ".join("%s=%.5f" % (k, v) for k, v in zip(KEYS, posture))
    keys += " " + " ".join("rule.%s=%.5f,%.5f" % (r, rules[2 * i], rules[2 * i + 1]) for i, r in enumerate(REGIONS))
    keys += " anat.scale=%.5f anat.root=%.5f" % (rules[-2], rules[-1])
    probe.stdin.write(keys + "\n"); probe.stdin.flush()
    head = probe.stdout.readline().split()
    arcs, areas = [], []
    while True:
        l = probe.stdout.readline()
        if l.startswith("end"): break
        f = l.split(); arcs.append(float(f[1])); areas.append(float(f[2]))
    return float(head[1]), numpy.array(arcs), numpy.array(areas)


def residuals(p):
    rules = p[:NSHARED]; out = []
    for i, v in enumerate(VOWELS):
        posture = p[NSHARED + 9 * i:NSHARED + 9 + 9 * i]
        L, A = story[v]
        length, arcs, areas = shape(posture, rules)
        centres = (numpy.arange(len(A)) + 0.5) * L
        model = numpy.interp(centres, arcs, areas, left=areas[0], right=areas[-1]) if len(arcs) else numpy.zeros(len(A))
        out.extend(numpy.log((model + 0.1) / (A + 0.1)))
        out.append(1.5 * (length - L * len(A)))
    return numpy.array(out)


import os
if os.path.exists("/tmp/es/geofit1.txt"):
    for l in open("/tmp/es/geofit1.txt"):
        f = l.split("#")[0].split()
        if f and f[0] in START:
            START[f[0]] = [float(x.split("=")[1]) for x in f[1:]]
p0 = numpy.array(RULE_START + sum((START[v] for v in VOWELS), []), float)
low = numpy.array([0.3, 0.9] * len(REGIONS) + [0.85, 0.0] + LOW * len(VOWELS)); high = numpy.array([5.0, 2.2] * len(REGIONS) + [1.15, 1.0] + HIGH * len(VOWELS))
r0 = residuals(p0)
print("start rms log error %.3f" % numpy.sqrt(numpy.mean(r0 ** 2)), flush=True)
fit = least_squares(residuals, p0, bounds=(low, high), diff_step=1e-3, max_nfev=400, verbose=0)
r = residuals(fit.x)
print("fitted rms log error %.3f" % numpy.sqrt(numpy.mean(r ** 2)))
with open(sys.argv[1], "w") as out:
    for i, reg in enumerate(REGIONS):
        out.write("rule.%s=%.4f,%.4f\n" % (reg, fit.x[2 * i], fit.x[2 * i + 1]))
    out.write("anat.scale=%.4f\nanat.root=%.4f\n" % (fit.x[NSHARED - 2], fit.x[NSHARED - 1]))
    for i, v in enumerate(VOWELS):
        posture = fit.x[NSHARED + 9 * i:NSHARED + 9 + 9 * i]
        L, A = story[v]
        res = r[i * (len(A) + 1):(i + 1) * (len(A) + 1)]
        out.write("%s %s  # rms %.3f length %+.2f cm\n" % (v, " ".join("%s=%.3f" % (k, x) for k, x in zip(KEYS, posture)),
                                                       numpy.sqrt(numpy.mean(res[:-1] ** 2)), res[-1] / 1.5))
print(open(sys.argv[1]).read())
