"""Vowel-to-vowel variation, ours against bdl: per vowel instance (middle
50%), level dB re the utterance's loudest vowel, and tilt (2-5 kHz re 0-1
kHz, dB). Prints the spread (std) of level and of tilt within utterances,
and the level-tilt correlation (louder vowels brighter in real speech).
usage: vlev.py DIR SUFFIX"""
import sys, glob, os, re, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly
R = os.path.expanduser("~/ducklips/voice")
VOW = {"aa", "ae", "ah", "ao", "aw", "ay", "eh", "er", "ey", "ih", "iy", "ow", "oy", "uh", "uw", "ax"}
def stats(items):
    sl, st, ll, tt = [], [], [], []
    for wav, ph in items:
        r, x = wavfile.read(wav); x = resample_poly(x.astype(float), 16000, r)
        lev, tilt = [], []
        for n, s, e in ph:
            if n in VOW and e - s > 0.05:
                seg = x[int((s + 0.25 * (e - s)) * 16000):int((e - 0.25 * (e - s)) * 16000)]
                S = numpy.abs(numpy.fft.rfft(seg * numpy.hanning(len(seg)))) ** 2; f = numpy.fft.rfftfreq(len(seg), 1 / 16000)
                lev.append(10 * numpy.log10(S.sum() / len(seg) + 1e-12))
                tilt.append(10 * numpy.log10(S[(f >= 2000) & (f < 5000)].sum() / S[f < 1000].sum()))
        if len(lev) > 3:
            lev = numpy.array(lev) - max(lev); tilt = numpy.array(tilt)
            sl.append(numpy.std(lev)); st.append(numpy.std(tilt)); ll.extend(lev - lev.mean()); tt.extend(tilt - tilt.mean())
    return numpy.mean(sl), numpy.mean(st), numpy.corrcoef(ll, tt)[0, 1]
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
    print("%-5s level spread %.1f dB  tilt spread %.1f dB  level-tilt r %+.2f" % (name, *stats(items)))
