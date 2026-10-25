"""Formants at the end of vowels before a velar (k g ng) and before a labial
(p b m): LPC F2, F3 over the vowel's last 25 ms against its middle; bdl
(lab: end times) vs ours (score). Velar pinch: F3-F2 closes at the end.
usage: pinch.py bdl LABDIR WAVDIR N | ours DIR"""
import sys, os, re, glob, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, lfilter
from scipy.linalg import solve_toeplitz
VOW = {"aa", "ae", "ah", "ao", "aw", "ay", "eh", "er", "ey", "ih", "iy", "ow", "oy", "uh", "uw", "ax"}
def load(p):
    r, x = wavfile.read(p); x = x.astype(float)
    x = resample_poly(x, 10000, r); return lfilter([1, -0.63], [1], x)
def formants(y):
    if len(y) < 200: return None
    f = y[:250] * numpy.hamming(min(250, len(y)))[:len(y[:250])]; c = numpy.correlate(f, f, "full")[len(f) - 1:len(f) + 12]
    if c[0] <= 0: return None
    a = solve_toeplitz(c[:12], -c[1:13]); r = numpy.roots(numpy.r_[1, a]); r = r[numpy.imag(r) > 0]
    fr = numpy.angle(r) * 10000 / (2 * numpy.pi); bw = -numpy.log(numpy.abs(r)) * 10000 / numpy.pi
    q = sorted(x for x, b in zip(fr, bw) if x > 200 and b < 700)
    return q[:3] if len(q) >= 3 else None
def items_bdl(labdir, wavdir, n):
    for lab in sorted(glob.glob(labdir + "/*.lab"))[:n]:
        rows = [l.split() for l in open(lab) if len(l.split()) == 3]; prev = 0; ph = []
        for r in rows: ph.append((r[2], prev, float(r[0]))); prev = float(r[0])
        yield os.path.join(wavdir, os.path.basename(lab)[:-4] + ".wav"), ph
def items_ours(d):
    for sc in sorted(glob.glob(d + "/*.score")):
        w = sc[:-6] + ("_ours.wav" if os.path.exists(sc[:-6] + "_ours.wav") else ".wav")
        if not os.path.exists(w): continue
        ph = [(m.group(1).rstrip("-0123456789").lower(), float(m.group(2)), float(m.group(2)) + float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sc)) if m]
        yield w, ph
src = items_bdl(sys.argv[2], sys.argv[3], int(sys.argv[4])) if sys.argv[1] == "bdl" else items_ours(sys.argv[2])
res = {"velar": [], "labial": [], "alveolar": []}
for w, ph in src:
    x = load(w)
    for (v, s, e), (c, _, _) in zip(ph[:-1], ph[1:]):
        if v not in VOW or e - s < 0.06: continue
        kind = "velar" if c in ("k", "g", "ng") else "labial" if c in ("p", "b", "m") else "alveolar" if c in ("t", "d", "n") else None
        if not kind: continue
        mid = formants(x[int((s + 0.4 * (e - s)) * 10000):int((s + 0.4 * (e - s)) * 10000) + 250])
        end = formants(x[int((e - 0.027) * 10000):int((e - 0.002) * 10000)])
        if mid and end: res[kind].append([end[1] - mid[1], end[2] - mid[2], (end[2] - end[1]) - (mid[2] - mid[1])])
for k, v in res.items():
    v = numpy.array(v)
    print("%-6s %-8s n %3d  dF2 %+5.0f  dF3 %+5.0f  d(F3-F2) %+5.0f Hz" % (sys.argv[1], k, len(v), *numpy.median(v, axis=0)))
