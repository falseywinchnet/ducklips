"""Sentence prosody: per file, F0 every 10 ms where voiced (normalized
autocorrelation >= 0.6, 40 ms window), in semitones re the file's median.
Prints medians over files of: F0 standard deviation (st), 5-95% range (st),
speaking rate (voiced frames per second of speech span), and the median
F0 (Hz); also the share of voiced time.
usage: prosody.py LABEL GLOB [LABEL GLOB ...]"""
import sys, glob, numpy
sys.path.insert(0, "/tmp")
from vw import load


def track(x):
    out = []
    for i in range(0, len(x) - 640 - 270, 160):
        a = x[i:i + 640] - x[i:i + 640].mean()
        ea = numpy.dot(a, a)
        if ea <= 0:
            out.append(0.0); continue
        c = numpy.correlate(x[i:i + 640 + 267] - x[i:i + 640 + 267].mean(), a, "valid")[::-1] if False else None
        best, lag = 0.0, 0
        full = numpy.correlate(a, a, "full")[639:]
        lo, hi = 53, 267
        k = lo + int(numpy.argmax(full[lo:hi]))
        r = full[k] / full[0]
        out.append(16000.0 / k if r >= 0.6 else 0.0)
    return numpy.array(out)


for label, pattern in zip(sys.argv[1::2], sys.argv[2::2]):
    sds, ranges, meds, voiced = [], [], [], []
    for p in sorted(glob.glob(pattern)):
        f = track(load(p))
        v = f[f > 0]
        if len(v) < 20:
            continue
        st = 12 * numpy.log2(v / numpy.median(v))
        sds.append(st.std()); ranges.append(numpy.percentile(st, 95) - numpy.percentile(st, 5)); meds.append(numpy.median(v))
        idx = numpy.where(f > 0)[0]
        voiced.append(len(v) / (idx[-1] - idx[0] + 1))
    print("%-8s files %d  F0 sd %.2f st  range %.1f st  median %.0f Hz  voiced share %.2f" % (
        label, len(sds), numpy.median(sds), numpy.median(ranges), numpy.median(meds), numpy.median(voiced)))
