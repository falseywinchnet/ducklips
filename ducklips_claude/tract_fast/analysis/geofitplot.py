"""Plot every vowel of a geofit result file against Story's MRI areas.
usage: geofitplot.py FIT.txt OUT.png"""
import sys, subprocess, numpy, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PROBE = "/Users/joshuahkuttenkuler/ducklips/tract/build/geometry_probe"
story = {}
for l in open("/Users/joshuahkuttenkuler/ducklips/tract/data/story_2008_areas.txt"):
    f = l.split()
    if f and not f[0].startswith("#"): story[f[0]] = (float(f[1]), numpy.array([float(v) for v in f[2:]]))
lines = [l.split("#")[0].strip() for l in open(sys.argv[1]) if l.strip()]
rules = " ".join(l for l in lines if l.startswith("rule") or l.startswith("anat"))
vowels = [l for l in lines if not (l.startswith("rule") or l.startswith("anat"))]
fig, ax = plt.subplots(2, len(vowels), figsize=(3.2 * len(vowels), 6.5), squeeze=False)
for k, l in enumerate(vowels):
    v, keys = l.split(" ", 1)
    out = subprocess.run([PROBE, "outline"], input=keys + " " + rules + "\n", capture_output=True, text=True).stdout.split("\n")
    hull = numpy.array([[float(a) for a in x.split()[1:]] for x in out if x.startswith("hull")])
    organ = numpy.array([[float(a) for a in x.split()[1:]] for x in out if x.startswith("organ")])
    cuts = [x.split() for x in out if x.startswith("cut")]
    a = ax[0, k]; a.plot(hull[:, 0], hull[:, 1], "k-"); a.plot(numpy.r_[organ[:, 0], organ[0, 0]], numpy.r_[organ[:, 1], organ[0, 1]], "r-")
    a.plot([float(c[5]) for c in cuts], [float(c[6]) for c in cuts], "b.", ms=2); a.set_aspect("equal"); a.set_xlim(-1, 11); a.set_ylim(-2, 11.5)
    a.set_title(v + " L %.1f" % float(out[0].split()[1]), fontsize=9); a.set_xticks([]); a.set_yticks([])
    b = ax[1, k]; b.step([float(c[1]) for c in cuts], [float(c[2]) for c in cuts], where="mid")
    L, A = story[v]; b.step((numpy.arange(len(A)) + 0.5) * L, A, where="mid"); b.set_ylim(0, 7); b.set_xlim(0, 21)
plt.tight_layout(); plt.savefig(sys.argv[2], dpi=55)
