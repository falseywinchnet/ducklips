"""Pitch steps at nasal boundaries: F0 every 5 ms (autocorrelation, 40 ms
frames, parabolic peak); at each vowel-nasal or nasal-vowel boundary the
largest |F0 change| between frames 10 ms apart within +-25 ms, in percent;
and the same statistic in vowel middles. usage: f0seam.py DIR SUFFIX [x]"""
import sys, glob, os, re, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly
R = os.path.expanduser("~/ducklips/voice")
VOW = {"aa","ae","ah","ao","aw","ay","eh","er","ey","ih","iy","ow","oy","uh","uw","ax"}
NAS = {"m","n","ng"}
def f0track(x, r):
    hop = int(0.005 * r); n = int(0.04 * r); out = []
    for s in range(0, len(x) - n, hop):
        seg = x[s:s + n] - x[s:s + n].mean()
        if numpy.dot(seg, seg) < 1e-6: out.append(0.0); continue
        c = numpy.correlate(seg, seg, "full")[n - 1:]; lo, hi = int(r / 250), int(r / 65)
        k = lo + numpy.argmax(c[lo:hi])
        if c[k] < 0.4 * c[0] or k <= lo or k >= hi - 1: out.append(0.0); continue
        a, b, g = c[k - 1], c[k], c[k + 1]; den = a - 2 * b + g
        out.append(r / (k + (0.5 * (a - g) / den if den != 0 else 0)))
    return numpy.array(out)
def run(items):
    at, mid = [], []
    for wav, ph in items:
        r, x = wavfile.read(wav); x = x.astype(float)
        if x.ndim > 1: x = x[:, 0]
        x = resample_poly(x, 16000, r); f = f0track(x, 16000)
        t = numpy.arange(len(f)) * 0.005 + 0.02
        ok = (f[2:] > 0) & (f[:-2] > 0)
        ch = numpy.where(ok, numpy.abs(f[2:] - f[:-2]) / numpy.maximum(f[:-2], 1) * 100, numpy.nan); tc = t[1:-1]
        for (a, s1, e1), (b, s2, e2) in zip(ph, ph[1:]):
            if (a in VOW and b in NAS) or (a in NAS and b in VOW):
                m = (numpy.abs(tc - e1) < 0.025) & ~numpy.isnan(ch)
                if m.sum() >= 3: at.append(numpy.nanmax(ch[m]))
        for n, s, e in ph:
            if n in VOW and e - s > 0.1:
                m = (tc > s + 0.3 * (e - s)) & (tc < e - 0.3 * (e - s)) & ~numpy.isnan(ch)
                if m.sum() >= 3: mid.append(numpy.nanmax(ch[m]))
    return numpy.median(at), numpy.percentile(at, 90), numpy.median(mid), len(at)
bdl = []
for lab in sorted(glob.glob(R + "/ref/arctic_bdl/lab/*.lab"))[:60]:
    rows = [l.split() for l in open(lab) if len(l.split()) == 3]; prev = 0; ph = []
    for r_ in rows: ph.append((r_[2], prev, float(r_[0]))); prev = float(r_[0])
    bdl.append((R + "/ref/arctic_bdl/wav/" + os.path.basename(lab)[:-4] + ".wav", ph))
ours = []
for sc in sorted(glob.glob(sys.argv[1] + "/arctic_a00[0-5]?.score")):
    ph = [(m.group(1).rstrip("-0123456789").lower(), float(m.group(2)), float(m.group(2)) + float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sc)) if m]
    ours.append((sc[:-6] + "_" + sys.argv[2] + ".wav", ph))
if len(sys.argv) < 4:
    print("bdl   F0 step at nasal boundaries median %.2f%%  p90 %.2f%%   in vowel middles %.2f%%  (n %d)" % run(bdl))
print("%-5s F0 step at nasal boundaries median %.2f%%  p90 %.2f%%   in vowel middles %.2f%%  (n %d)" % (sys.argv[2], *run(ours)))
