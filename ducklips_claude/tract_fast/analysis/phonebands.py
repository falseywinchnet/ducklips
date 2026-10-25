"""Every phone's spectrum against bdl's: median level in 7 bands (0-0.5,
0.5-1, 1-2, 2-3, 3-4, 4-6, 6-8 kHz) of each phone's middle 50%, dB re the
utterance's loudest vowel; ours (voice/out/arctic renders, a0001-a0060)
minus bdl's, sorted by the mean absolute difference.
usage: phonebands.py"""
import os, re, glob, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly
R = os.path.expanduser("~/ducklips/voice")
B = [(0, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, 4000), (4000, 6000), (6000, 8000)]
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
    got = {}
    for wav, phones in items:
        x = load(wav); ref = -1e9
        for n, s, e in phones:
            if n in VOW and e - s > 0.05:
                b = bands(x, s, e)
                if b is not None: ref = max(ref, 10 * numpy.log10(numpy.sum(10 ** (b / 10))))
        for n, s, e in phones:
            if e - s > 0.03 and n not in ("pau", "h#", "ssil"):
                b = bands(x, s, e)
                if b is not None: got.setdefault(n, []).append(b - ref)
    return {k: (numpy.median(v, axis=0), len(v)) for k, v in got.items()}
bdl = []
for lab in sorted(glob.glob(R + "/ref/arctic_bdl/lab/*.lab")):
    rows = [l.split() for l in open(lab) if len(l.split()) == 3]; prev = 0; ph = []
    for r in rows: ph.append((r[2], prev, float(r[0]))); prev = float(r[0])
    bdl.append((R + "/ref/arctic_bdl/wav/" + os.path.basename(lab)[:-4] + ".wav", ph))
ours = []
for sc in sorted(glob.glob(R + "/out/arctic/arctic_a00[0-5]?.score")):
    ph = [(m.group(1).rstrip("-0123456789").lower(), float(m.group(2)), float(m.group(2)) + float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sc)) if m]
    ours.append((sc[:-6] + "_ours.wav", ph))
b, o = run(bdl), run(ours)
rows = []
for k in b:
    if k in o and b[k][1] >= 8 and o[k][1] >= 8:
        d = o[k][0] - b[k][0]; rows.append((numpy.mean(numpy.abs(d)), k, d, b[k][0]))
print("phone  mean|d|   ours - bdl (dB): 0-.5 .5-1 1-2 2-3 3-4 4-6 6-8     bdl levels")
for m, k, d, lv in sorted(rows, reverse=True):
    print("%-4s %6.1f   %s     %s" % (k, m, " ".join("%+5.0f" % v for v in d), " ".join("%4.0f" % v for v in lv)))
