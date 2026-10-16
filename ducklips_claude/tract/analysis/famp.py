"""Formant amplitudes per vowel, ours against bdl: medians of H1-A1, A1-A2,
A2-A3, A3-A4 (dB; Ak the spectrum's peak within 100 Hz of LPC formant k)
and LPC bandwidths B1-B4 (Hz, 25 ms frames; biased, compared like for like).
usage: famp.py DIR SUFFIX   (ours: DIR/arctic_a00NN.score + _SUFFIX.wav)"""
import sys, glob, os, re, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, lfilter
from scipy.linalg import solve_toeplitz
R = os.path.expanduser("~/ducklips/voice")
V = ["iy", "ih", "eh", "ae", "aa", "ah", "ao", "ey", "ay"]
def load(p):
    r, x = wavfile.read(p); x = x.astype(float)
    return resample_poly(x, 10000, r), resample_poly(x, 16000, r)
def formants(y):
    y = lfilter([1, -0.63], [1], y); out = []
    for i in range(0, len(y) - 250, 100):
        f = y[i:i + 250] * numpy.hamming(250); c = numpy.correlate(f, f, "full")[249:263]
        if c[0] <= 0: continue
        a = solve_toeplitz(c[:12], -c[1:13]); rt = numpy.roots(numpy.r_[1, a]); rt = rt[numpy.imag(rt) > 0]
        fr = numpy.angle(rt) * 10000 / (2 * numpy.pi); bw = -numpy.log(numpy.abs(rt)) * 10000 / numpy.pi
        q = sorted((a_, b_) for a_, b_ in zip(fr, bw) if a_ > 150 and b_ < 700)[:4]
        if len(q) == 4: out.append(q)
    return out
def measure(items):
    got = {}
    for wav, ph in items:
        x10, x16 = load(wav)
        for n, s, e in ph:
            if n not in V or e - s < 0.07: continue
            a, b = s + 0.3 * (e - s), e - 0.3 * (e - s)
            fs = formants(x10[int(a * 10000):int(b * 10000)])
            if not fs: continue
            F = numpy.median(numpy.array(fs), axis=0)
            seg = x16[int(a * 16000):int(b * 16000)]
            S = 20 * numpy.log10(numpy.abs(numpy.fft.rfft(seg * numpy.hanning(len(seg)), 8192)) + 1e-9)
            fr = numpy.fft.rfftfreq(8192, 1 / 16000)
            def peak(lo, hi): return S[(fr >= lo) & (fr < hi)].max()
            A = [peak(f - 100, f + 100) for f, _ in F]
            H1 = peak(70, 190)
            got.setdefault(n, []).append([H1 - A[0], A[0] - A[1], A[1] - A[2], A[2] - A[3]] + [b_ for _, b_ in F] + [f for f, _ in F])
    return {k: numpy.median(numpy.array(v), axis=0) for k, v in got.items()}
bdl = []
for lab in sorted(glob.glob(R + "/ref/arctic_bdl/lab/*.lab")):
    rows = [l.split() for l in open(lab) if len(l.split()) == 3]; prev = 0; ph = []
    for r in rows: ph.append((r[2], prev, float(r[0]))); prev = float(r[0])
    bdl.append((R + "/ref/arctic_bdl/wav/" + os.path.basename(lab)[:-4] + ".wav", ph))
ours = []
for sc in sorted(glob.glob(sys.argv[1] + "/arctic_a00[0-5]?.score")):
    ph = [(m.group(1).rstrip("-0123456789").lower(), float(m.group(2)), float(m.group(2)) + float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sc)) if m]
    ours.append((sc[:-6] + "_" + sys.argv[2] + ".wav", ph))
b, o = measure(bdl), measure(ours)
print("vowel      H1-A1 A1-A2 A2-A3 A3-A4 |    B1   B2   B3   B4 |   F4")
for v in V:
    if v in b and v in o:
        print("%-3s ours  %+5.0f %+5.0f %+5.0f %+5.0f | %5.0f %4.0f %4.0f %4.0f | %4.0f" % (v, *o[v][:8], o[v][11]))
        print("    bdl   %+5.0f %+5.0f %+5.0f %+5.0f | %5.0f %4.0f %4.0f %4.0f | %4.0f" % (*b[v][:8], b[v][11]))
m = numpy.mean([o[v] - b[v] for v in V if v in b and v in o], axis=0)
print("mean ours-bdl %+5.1f %+5.1f %+5.1f %+5.1f | %+5.0f %+4.0f %+4.0f %+4.0f | %+4.0f" % (*m[:8], m[11]))
