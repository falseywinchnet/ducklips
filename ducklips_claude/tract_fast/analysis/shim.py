"""Local shimmer in vowels, ours against bdl: each vowel's middle 60%,
low-passed at 1 kHz; period from autocorrelation; the peak-to-peak
amplitude of each successive period; local shimmer = mean |A(i+1)-A(i)| /
mean A (%), and in dB. usage: shim.py DIR SUFFIX"""
import sys, glob, os, re, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, butter, sosfiltfilt
R = os.path.expanduser("~/ducklips/voice")
VOW = {"aa", "ae", "ah", "ao", "aw", "ay", "eh", "ey", "ih", "iy", "ow", "uw", "er"}
SOS = butter(4, 1000, fs=16000, output="sos")
def shimmer(items):
    loc, db = [], []
    for wav, ph in items:
        r, x = wavfile.read(wav); x = resample_poly(x.astype(float), 16000, r); y = sosfiltfilt(SOS, x)
        for n, s, e in ph:
            if n not in VOW or e - s < 0.08: continue
            seg = y[int((s + 0.2 * (e - s)) * 16000):int((e - 0.2 * (e - s)) * 16000)]
            c = numpy.correlate(seg, seg, "full")[len(seg) - 1:]
            lo, hi = 16000 // 250, 16000 // 70
            if len(c) <= hi: continue
            T = lo + numpy.argmax(c[lo:hi])
            if c[T] < 0.5 * c[0]: continue
            A = numpy.array([seg[i:i + T].max() - seg[i:i + T].min() for i in range(0, len(seg) - T, T)])
            if len(A) < 4: continue
            loc.append(numpy.mean(numpy.abs(numpy.diff(A))) / numpy.mean(A) * 100)
            db.append(numpy.mean(numpy.abs(20 * numpy.log10(A[1:] / A[:-1]))))
    return numpy.median(loc), numpy.median(db), len(loc)
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
    print("%-5s local shimmer %.1f %%  %.2f dB  (n %d)" % (name, *shimmer(items)))
