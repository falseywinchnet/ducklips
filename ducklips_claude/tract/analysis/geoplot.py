"""Draw midsagittal shapes and their area functions from geometry_probe outline output.
usage: geoplot.py OUT.png "label|key=value ..." ...   (MRI=vowel overlays Story's area function)"""
import sys, subprocess, numpy, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PROBE = sys.argv[1].rsplit("/", 1)[0] + "/geometry_probe" if False else "/Users/joshuahkuttenkuler/ducklips/tract/build/geometry_probe"
story = {}
for l in open("/Users/joshuahkuttenkuler/ducklips/tract/data/story_2008_areas.txt"):
    f = l.split()
    if f and not f[0].startswith("#"): story[f[0]] = (float(f[1]), numpy.array([float(v) for v in f[2:]]))
specs = sys.argv[2:]
fig, ax = plt.subplots(2, len(specs), figsize=(4.5 * len(specs), 8), squeeze=False)
for k, spec in enumerate(specs):
    label, keys = spec.split("|")
    out = subprocess.run([PROBE, "outline"], input=keys + "\n", capture_output=True, text=True).stdout.split("\n")
    hull = numpy.array([[float(a) for a in l.split()[1:]] for l in out if l.startswith("hull")])
    organ = numpy.array([[float(a) for a in l.split()[1:]] for l in out if l.startswith("organ")])
    cuts = [l.split() for l in out if l.startswith("cut")]
    head = out[0]
    a = ax[0, k]; a.plot(hull[:, 0], hull[:, 1], "k-"); a.plot(numpy.r_[organ[:, 0], organ[0, 0]], numpy.r_[organ[:, 1], organ[0, 1]], "r-")
    a.plot([float(c[5]) for c in cuts], [float(c[6]) for c in cuts], "b.", ms=3); a.set_aspect("equal"); a.set_xlim(-1, 11); a.set_ylim(-2, 11.5); a.set_title(label + "\n" + head, fontsize=8)
    b = ax[1, k]; b.step([float(c[1]) for c in cuts], [float(c[2]) for c in cuts], where="mid", label="model")
    if label in story:
        L, A = story[label]; b.step((numpy.arange(len(A)) + 0.5) * L, A, where="mid", label="Story 2008")
    b.set_ylim(0, 8); b.legend(fontsize=7); b.set_xlabel("cm from glottis")
plt.tight_layout(); plt.savefig(sys.argv[1], dpi=60)
