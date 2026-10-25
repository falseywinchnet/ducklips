"""Seams at nasal boundaries: log spectra (0.1-5 kHz, 30 ms frames, 5 ms
hop); spectral change = mean |dB difference| between frames 10 ms apart.
At each vowel-nasal or nasal-vowel boundary, the largest change within
+-25 ms, re the median change in the middles of vowels of that utterance.
Ours against bdl. usage: seams.py DIR SUFFIX"""
import sys, glob, os, re, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly
R = os.path.expanduser("~/ducklips/voice")
VOW = {"aa","ae","ah","ao","aw","ay","eh","er","ey","ih","iy","ow","oy","uh","uw","ax"}
NAS = {"m","n","ng"}
def spectra(x, r):
    n = int(0.03 * r); hop = int(0.005 * r); win = numpy.hanning(n); f = numpy.fft.rfftfreq(n, 1 / r); band = (f > 100) & (f < 5000)
    S = [20 * numpy.log10(numpy.abs(numpy.fft.rfft(x[s:s + n] * win))[band] + 1e-6) for s in range(0, len(x) - n, hop)]
    return numpy.array(S)
def run(items):
    ratios = []
    for wav, ph in items:
        r, x = wavfile.read(wav); x = x.astype(float)
        if x.ndim > 1: x = x[:, 0]
        x = resample_poly(x, 16000, r); r = 16000
        S = spectra(x, r)
        if len(S) < 5: continue
        ch = numpy.abs(S[2:] - S[:-2]).mean(axis=1); tc = (numpy.arange(len(ch)) + 1) * 0.005 + 0.015
        inner = []
        for n, s, e in ph:
            if n in VOW and e - s > 0.08:
                m = (tc > s + 0.3 * (e - s)) & (tc < e - 0.3 * (e - s)); inner.extend(ch[m])
        if not inner: continue
        base = numpy.median(inner)
        for (a, s1, e1), (b, s2, e2) in zip(ph, ph[1:]):
            if (a in VOW and b in NAS) or (a in NAS and b in VOW):
                m = numpy.abs(tc - e1) < 0.025
                if m.any(): ratios.append(ch[m].max() / base)
    return numpy.median(ratios), numpy.percentile(ratios, 90), len(ratios)
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
    print("bdl   seam ratio median %.2f  p90 %.2f  (n %d)" % run(bdl))
print("ours  seam ratio median %.2f  p90 %.2f  (n %d)  [%s]" % (*run(ours), sys.argv[2]))
