"""Band levels of nasal murmurs (middle 50%) re the loudest vowel of the
utterance, dB, in 0-0.5, 0.5-1, 1-2, 2-4, 4-8 kHz; bdl (lab: end times)
against ours (score '# NAME stress start dur'). Medians over all n/m/ng.
usage: nasalband.py bdl LABDIR WAVDIR N | ours DIR (score+wav pairs) [VAR...]"""
import sys, os, re, glob, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly
B = [(0, 500), (500, 1000), (1000, 2000), (2000, 4000), (4000, 8000)]
VOW = {"aa", "ae", "ah", "ao", "aw", "ay", "eh", "er", "ey", "ih", "iy", "ow", "oy", "uh", "uw", "ax"}
def load(p):
    r, x = wavfile.read(p); x = x.astype(float)
    return resample_poly(x, 16000, r) if r != 16000 else x
def bands(x, s, e):
    seg = x[int((s + 0.25 * (e - s)) * 16000):int((e - 0.25 * (e - s)) * 16000)]
    if len(seg) < 160: return None
    S = numpy.abs(numpy.fft.rfft(seg * numpy.hanning(len(seg)))) ** 2 / len(seg); f = numpy.fft.rfftfreq(len(seg), 1 / 16000)
    return numpy.array([10 * numpy.log10(S[(f >= a) & (f < b)].sum() + 1e-12) for a, b in B])
def run(items):
    out = []
    for wav, phones in items:
        x = load(wav); ref = -1e9
        for n, s, e in phones:
            if n in VOW and e - s > 0.05:
                b = bands(x, s, e)
                if b is not None: ref = max(ref, 10 * numpy.log10(numpy.sum(10 ** (b / 10))))
        for n, s, e in phones:
            if n in ("n", "m", "ng") and e - s > 0.04:
                b = bands(x, s, e)
                if b is not None: out.append(b - ref)
    return numpy.median(out, axis=0), len(out)
if sys.argv[1] == "bdl":
    items = []
    for lab in sorted(glob.glob(sys.argv[2] + "/*.lab"))[:int(sys.argv[4])]:
        rows = [l.split() for l in open(lab) if len(l.split()) == 3]; prev = 0; ph = []
        for r in rows: ph.append((r[2], prev, float(r[0]))); prev = float(r[0])
        items.append((os.path.join(sys.argv[3], os.path.basename(lab)[:-4] + ".wav"), ph))
else:
    items = []
    for sc in sorted(glob.glob(sys.argv[2] + "/*.score")):
        w = sc[:-6] + ".wav" if os.path.exists(sc[:-6] + ".wav") else sc[:-6] + "_ours.wav"
        if not os.path.exists(w): continue
        ph = [(m.group(1).rstrip("-0123456789").lower(), float(m.group(2)), float(m.group(2)) + float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sc)) if m]
        items.append((w, ph))
m, n = run(items)
print("%-10s n %3d  %s" % (sys.argv[1] if sys.argv[1] == "bdl" else sys.argv[2].split("/")[-1], n, "  ".join("%5.1f" % v for v in m)))
