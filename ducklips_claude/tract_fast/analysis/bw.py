"""Median LPC F1-F3 per vowel (middle 40%), bdl (a0001-a0060) against ours
(voice/out/arctic renders). usage: f123.py"""
import sys, glob, os, re, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, lfilter
from scipy.linalg import solve_toeplitz
R = os.path.expanduser("~/ducklips/voice")
V = ["iy", "ih", "eh", "ae", "aa", "ah", "ao", "ow", "uh", "uw", "er", "ey", "ay"]
def load(p):
    r, x = wavfile.read(p); return lfilter([1, -0.63], [1], resample_poly(x.astype(float), 10000, r))
def lpc3(y):
    out = []
    for i in range(0, len(y) - 250, 100):
        f = y[i:i + 250] * numpy.hamming(250); c = numpy.correlate(f, f, "full")[249:262]
        if c[0] <= 0: continue
        a = solve_toeplitz(c[:12], -c[1:13]); rt = numpy.roots(numpy.r_[1, a]); rt = rt[numpy.imag(rt) > 0]
        fr = numpy.angle(rt) * 10000 / (2 * numpy.pi); bw = -numpy.log(numpy.abs(rt)) * 10000 / numpy.pi
        q = sorted((a_, b_) for a_, b_ in zip(fr, bw) if a_ > 150 and b_ < 900)[:4]
        if len(q) == 4: out.append([x for pair in q for x in pair])
    return out
def collect(items):
    got = {}
    for wav, ph in items:
        x = load(wav)
        for n, s, e in ph:
            if n in V and e - s > 0.06:
                got.setdefault(n, []).extend(lpc3(x[int((s + 0.3 * (e - s)) * 10000):int((e - 0.3 * (e - s)) * 10000)]))
    return {k: numpy.median(v, axis=0) for k, v in got.items() if v}
bdl = []
for lab in sorted(glob.glob(R + "/ref/arctic_bdl/lab/*.lab")):
    rows = [l.split() for l in open(lab) if len(l.split()) == 3]; prev = 0; ph = []
    for r in rows: ph.append((r[2], prev, float(r[0]))); prev = float(r[0])
    bdl.append((R + "/ref/arctic_bdl/wav/" + os.path.basename(lab)[:-4] + ".wav", ph))
ours = []
for sc in sorted(glob.glob(R + "/out/arctic/arctic_a00[0-5]?.score")):
    ph = [(m.group(1).rstrip("-0123456789").lower(), float(m.group(2)), float(m.group(2)) + float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sc)) if m]
    ours.append((sc[:-6] + "_ours.wav", ph))
b, o = collect(bdl), collect(ours)
print("vowel   ours F1 B1 F2 B2 F3 B3 F4 B4          bdl")
for v in V:
    if v in b and v in o:
        print("%-3s  %s    %s" % (v, " ".join("%4.0f" % x for x in o[v]), " ".join("%4.0f" % x for x in b[v])))
mb = numpy.median([b[v] for v in V if v in b and v in o], axis=0); mo = numpy.median([o[v] for v in V if v in b and v in o], axis=0)
print("median B1-B4 ours %s  bdl %s" % (" ".join("%4.0f" % x for x in mo[1::2]), " ".join("%4.0f" % x for x in mb[1::2])))
