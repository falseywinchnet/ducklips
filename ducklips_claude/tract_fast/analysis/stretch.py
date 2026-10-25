"""Stretch the movement times of chosen parts in a score (tongue body, jaw,
opening or rounding lips, pitch; closures keep their speed), centred on each
event's original midpoint: a move starting at t over m seconds becomes one
starting at t - (k - 1) m / 2 over k m. Pitch takes its own factor.
usage: stretch.py IN.score OUT.score ARTIC_FACTOR PITCH_FACTOR"""
import sys
inp, out, ka, kp = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
ARTIC = {"tongue", "jaw"}  # plus lips moves that do not close (target over 0.5)
lines = []
for l in open(inp):
    f = l.split()
    if len(f) == 5 and not l.startswith("#") and (f[1] in ARTIC or f[1] == "pitch" or (f[1] == "lips" and float(f[2]) > 0.5)):
        k = kp if f[1] == "pitch" else ka
        t, m = float(f[0]), float(f[4])
        f[0] = "%.4f" % max(0.0, t - (k - 1) * m / 2); f[4] = "%.3f" % (k * m)
        l = " ".join(f) + "\n"
    lines.append(l)
head = [l for l in lines if l.startswith("#") or len(l.split()) != 5]
body = sorted([l for l in lines if not (l.startswith("#") or len(l.split()) != 5)], key=lambda l: float(l.split()[0]))
open(out, "w").write("".join(body + head))
