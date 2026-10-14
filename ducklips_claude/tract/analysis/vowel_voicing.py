"""Voicing in the vowels of a spoken score: the share of 10 ms frames inside
the score's vowels (phone list in the score's comments) that are periodic
(normalised autocorrelation >= 0.5 at a lag of 60-250 Hz), their median F0,
the share of frames whose F0 is within 25% of the score's planned pitch, and
per-vowel detail. usage: vowel_voicing.py SCORE WAV [--detail]"""
import sys, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly
score, wav = sys.argv[1], sys.argv[2]
phones = []
for line in open(score):
    p = line.split()
    if len(p) == 5 and p[0] == "#" and p[1][-1] in "012" and "-" not in p[1]:
        phones.append((p[1], float(p[2]), float(p[3]), float(p[4])))
r, x = wavfile.read(wav); x = resample_poly(x.astype(float), 16000, r); r = 16000
total = voiced = near = 0; f0s = []; rows = []
for name, start, dur, plan in phones:
    v = n = 0; fs = []
    for t in numpy.arange(start + 0.01, start + dur - 0.03, 0.01):
        s = x[int(t * r):int((t + 0.03) * r)]; s = s - s.mean()
        if len(s) < 100: continue
        c = numpy.correlate(s, s, "full")[len(s) - 1:]; lo, hi = int(r / 250), int(r / 60)
        k = lo + numpy.argmax(c[lo:hi]); n += 1
        if c[0] > 0 and c[k] >= 0.35 * c[0]:
            v += 1; f = r / k; fs.append(f)
            if plan > 0 and abs(f / plan - 1) < 0.25: near += 1
    total += n; voiced += v; f0s += fs
    rows.append("%s@%.2f %d/%d %s" % (name, start, v, n, "%.0f(%.0f)" % (numpy.median(fs), plan) if fs else "-"))
print("vowel frames voiced %.0f%% (%d/%d), median F0 %.0f, within 25%% of plan %.0f%%" % (
    100 * voiced / max(total, 1), voiced, total, numpy.median(f0s) if f0s else 0, 100 * near / max(total, 1)))
if "--detail" in sys.argv: print("  " + "  ".join(rows))
