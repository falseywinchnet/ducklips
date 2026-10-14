"""Plot heard F0 against the plan, subglottal pressure, ear correction and
velum through a stretch of a traced render.
usage: ptrace.py WAV TRACE SCORE FROM TO OUT.png"""
import sys, re, numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/tmp")
from vw import load, pitch_at
wav, tr, score, a, b, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
x = load(wav)
ts = numpy.arange(a, b, 0.005)
f0 = numpy.array([pitch_at(x, int((t - 0.015) * 16000))[0] for t in ts])
rows = [l.split() for l in open(tr)]
t = numpy.array([float(r[0]) for r in rows]); plan = numpy.array([float(r[3]) for r in rows]); psub = numpy.array([float(r[6]) for r in rows])
ear = numpy.array([float(r[r.index("ear") + 1]) for r in rows]); vel = numpy.array([float(r[r.index("velum") + 1]) for r in rows])
area = numpy.array([float(r[5]) for r in rows])
m = (t >= a) & (t < b)
fig, ax = plt.subplots(4, 1, figsize=(14, 10), sharex=True)
ax[0].plot(ts[f0 > 0], f0[f0 > 0], ".", ms=3, label="heard F0"); ax[0].plot(t[m], plan[m], label="plan"); ax[0].set_ylim(60, 220); ax[0].legend()
ax[1].plot(t[m], psub[m], label="Psub (Pa)"); ax[1].legend()
ax[2].plot(t[m], ear[m], label="ear correction"); ax[2].plot(t[m], vel[m], label="velum"); ax[2].legend()
ax[3].plot(t[m], area[m], label="glottal area col"); ax[3].legend()
for l in open(score):
    mm = re.match(r"# ([A-Z]+)-?\d? ([\d.]+) ([\d.]+)", l)
    if mm and mm.group(1) != "phones":
        s = float(mm.group(2))
        if a <= s < b:
            for k in ax: k.axvline(s, color="0.8", lw=0.5)
            ax[0].text(s + 0.005, 205, mm.group(1), fontsize=8)
plt.tight_layout(); plt.savefig(out, dpi=80)
