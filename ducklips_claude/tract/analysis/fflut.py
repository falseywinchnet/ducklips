"""Formant flutter in vowels, ours against bdl: LPC F1-F3 every 5 ms (25 ms
frames) over each vowel's middle 60%; the median |second difference| (Hz)
of each track, and the share of frames whose second difference exceeds
100 / 150 / 200 Hz (F1 / F2 / F3): wobble faster than articulation.
usage: fflut.py DIR SUFFIX"""
import sys, glob, os, re, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, lfilter
from scipy.linalg import solve_toeplitz
R = os.path.expanduser("~/ducklips/voice")
VOW = {"aa", "ae", "ah", "ao", "aw", "ay", "eh", "ey", "ih", "iy", "ow", "uw"}
def track(y):
    y = lfilter([1, -0.63], [1], y); out = []
    for i in range(0, len(y) - 250, 50):
        f = y[i:i + 250] * numpy.hamming(250); c = numpy.correlate(f, f, "full")[249:262]
        if c[0] <= 0: out.append(None); continue
        a = solve_toeplitz(c[:12], -c[1:13]); rt = numpy.roots(numpy.r_[1, a]); rt = rt[numpy.imag(rt) > 0]
        fr = numpy.angle(rt) * 10000 / (2 * numpy.pi); bw = -numpy.log(numpy.abs(rt)) * 10000 / numpy.pi
        q = sorted(a_ for a_, b_ in zip(fr, bw) if a_ > 200 and b_ < 500)[:3]
        out.append(q if len(q) == 3 else None)
    return out
def flutter(items):
    d2 = []
    for wav, ph in items:
        r, x = wavfile.read(wav); x = resample_poly(x.astype(float), 10000, r)
        for n, s, e in ph:
            if n in VOW and e - s > 0.09:
                t = track(x[int((s + 0.2 * (e - s)) * 10000):int((e - 0.2 * (e - s)) * 10000)])
                for a, b, c in zip(t, t[1:], t[2:]):
                    if a and b and c: d2.append([abs(c[k] - 2 * b[k] + a[k]) for k in range(3)])
    d2 = numpy.array(d2)
    return numpy.median(d2, axis=0), [(d2[:, k] > th).mean() for k, th in enumerate((100, 150, 200))], len(d2)
bdl = []
for lab in sorted(glob.glob(R + "/ref/arctic_bdl/lab/*.lab"))[:60]:
    rows = [l.split() for l in open(lab) if len(l.split()) == 3]; prev = 0; ph = []
    for r in rows: ph.append((r[2], prev, float(r[0]))); prev = float(r[0])
    bdl.append((R + "/ref/arctic_bdl/wav/" + os.path.basename(lab)[:-4] + ".wav", ph))
ours = []
for sc in sorted(glob.glob(sys.argv[1] + "/arctic_a00[0-5]?.score")):
    ph = [(m.group(1).rstrip("-0123456789").lower(), float(m.group(2)), float(m.group(2)) + float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sc)) if m]
    ours.append((sc[:-6] + "_" + sys.argv[2] + ".wav", ph))
for name, items in (("bdl", bdl), ("ours", ours)):
    m, share, n = flutter(items)
    print("%-5s median |d2| F1 %4.0f F2 %4.0f F3 %4.0f Hz   big jumps %.2f %.2f %.2f  (n %d)" % (name, *m, *share, n))
