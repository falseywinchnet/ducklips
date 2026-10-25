"""Release-aligned stop bursts: the release is the sharpest rise of 2-8 kHz
energy (2.5 ms frames) within 40 ms of the stop's end; the burst is the
15 ms from it. Level re the utterance's loudest vowel (same band), band
shape, and the rise itself (dB in 5 ms). Ours vs bdl over /tmp/es/corpus."""
import glob, os, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, butter, sosfilt
VOW = set("AA AE AH AO AW AX AY EH ER EY IH IY OW OY UH UW".split())
bands = [(500, 1500), (1500, 3000), (3000, 5000), (5000, 7800)]
def collect(path, spans, rows):
    r, x = wavfile.read(path); x = x.astype(float); x = x[:, 0] if x.ndim > 1 else x
    x = resample_poly(x, 16000, r); r = 16000
    hi = sosfilt(butter(4, 2000, "high", fs=r, output="sos"), x)
    top = max([10 * numpy.log10(numpy.mean(x[int(a * r):int(b * r)]**2) + 1e-9) for n, a, b in spans if n in VOW and b - a > 0.04] or [0])
    for n, a, b in spans:
        if n not in ("P", "T", "K"): continue
        lo_i = int((b - 0.04) * r); hi_i = int((b + 0.04) * r)
        if lo_i < 80 or hi_i + 400 > len(x): continue
        e = numpy.array([10 * numpy.log10(numpy.mean(hi[i:i + 40]**2) + 1e-9) for i in range(lo_i - 80, hi_i, 40)])
        rise = e[2:] - e[:-2]
        k = int(numpy.argmax(rise)); at = lo_i - 80 + (k + 2) * 40
        s = x[at:at + 240]
        S = numpy.abs(numpy.fft.rfft(s * numpy.hanning(len(s)), 2048))**2; f = numpy.fft.rfftfreq(2048, 1 / r)
        v = numpy.array([10 * numpy.log10(S[(f >= lo) & (f < hh)].sum() + 1e-9) for lo, hh in bands])
        rows.setdefault(n, []).append(numpy.concatenate([v - v.max(), [10 * numpy.log10(numpy.mean(s**2) + 1e-9) - top, rise[k]]]))
ours = {}; bdl = {}
for score in sorted(glob.glob("/tmp/es/corpus/*.score")):
    name = os.path.basename(score)[:-6]
    sp = [(l.split()[1].split("-")[0].rstrip("012"), float(l.split()[2]), float(l.split()[2]) + float(l.split()[3])) for l in open(score) if l.startswith("# ") and len(l.split()) == 5 and l.split()[1] != "phones:"]
    collect(score[:-6] + ".wav", sp, ours)
    prev = 0.0; sp = []
    for l in open(os.path.expanduser("~/ducklips/voice/ref/arctic_bdl_full/lab/%s.lab" % name)):
        p = l.split()
        if len(p) == 3 and p[0][0].isdigit(): sp.append((p[2].upper(), prev, float(p[0]))); prev = float(p[0])
    collect(os.path.expanduser("~/ducklips/voice/ref/arctic_bdl_full/wav/%s.wav" % name), sp, bdl)
for k in ("P", "T", "K"):
    for label, d in (("bdl", bdl), ("ours", ours)):
        m = numpy.median(numpy.array(d[k]), 0)
        print("%s %-4s n=%2d  shape %s  burst %+.0f dB  rise %.0f dB/5ms" % (k, label, len(d[k]), " ".join("%+.0f" % v for v in m[:4]), m[4], m[5]))
# Fricative shapes (middle 60%), dB re strongest band.
fb2 = [(500, 1500), (1500, 3000), (3000, 5000), (5000, 7800)]
def fshape(path, spans, rows):
    r, x = wavfile.read(path); x = x.astype(float); x = x[:, 0] if x.ndim > 1 else x
    x = resample_poly(x, 16000, r); r = 16000
    for n, a, b in spans:
        if n in ("S", "SH", "F", "Z") and b - a > 0.04:
            d = b - a; s = x[int((a + 0.2 * d) * r):int((b - 0.2 * d) * r)]
            S = numpy.abs(numpy.fft.rfft(s * numpy.hanning(len(s)), 4096))**2; f = numpy.fft.rfftfreq(4096, 1 / r)
            v = numpy.array([10 * numpy.log10(S[(f >= lo) & (f < hh)].sum() + 1e-9) for lo, hh in fb2]); rows.setdefault(n, []).append(v - v.max())
fo = {}; fb = {}
for score in sorted(glob.glob("/tmp/es/corpus/*.score")):
    name = os.path.basename(score)[:-6]
    sp = [(l.split()[1].split("-")[0].rstrip("012"), float(l.split()[2]), float(l.split()[2]) + float(l.split()[3])) for l in open(score) if l.startswith("# ") and len(l.split()) == 5 and l.split()[1] != "phones:"]
    fshape(score[:-6] + ".wav", sp, fo)
    prev = 0.0; sp = []
    for l in open(os.path.expanduser("~/ducklips/voice/ref/arctic_bdl_full/lab/%s.lab" % name)):
        p = l.split()
        if len(p) == 3 and p[0][0].isdigit(): sp.append((p[2].upper(), prev, float(p[0]))); prev = float(p[0])
    fshape(os.path.expanduser("~/ducklips/voice/ref/arctic_bdl_full/wav/%s.wav" % name), sp, fb)
for k in ("S", "SH", "F", "Z"):
    if k in fo and k in fb:
        print("%-2s shape ours %s | bdl %s" % (k, " ".join("%+.0f" % v for v in numpy.median(fo[k], 0)), " ".join("%+.0f" % v for v in numpy.median(fb[k], 0))))
