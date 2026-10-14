"""F0 statistics (autocorrelation, 10 ms) for ours vs ARCTIC across utterances."""
import os
import sys
import wave

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca


def track(path):
    data, rate = ca.read_wave(path)
    size = int(0.04 * rate)
    peak = numpy.abs(data).max()
    values = []
    for start in range(0, len(data) - size, int(0.01 * rate)):
        seg = data[start:start + size]
        if numpy.sqrt(numpy.mean(seg ** 2)) < peak * 0.05:
            values.append(numpy.nan)
            continue
        seg = seg - seg.mean()
        c = numpy.correlate(seg, seg, "full")[size - 1:]
        lo, hi = int(rate / 300), int(rate / 60)
        lag = lo + int(numpy.argmax(c[lo:hi]))
        values.append(rate / lag if c[lag] / c[0] > 0.5 else numpy.nan)
    return numpy.array(values)


def stats(values):
    v = values[~numpy.isnan(values)]
    semitones = 12 * numpy.log2(v / numpy.median(v))
    steps = numpy.diff(12 * numpy.log2(values))
    steps = steps[~numpy.isnan(steps) & (numpy.abs(steps) < 3)]
    return numpy.median(v), numpy.percentile(semitones, 90) - numpy.percentile(semitones, 10), numpy.mean(numpy.abs(steps)) * 100


rows = {"ref": [], "ours": []}
for number in range(1, 21):
    name = "arctic_a%04d" % number
    rows["ref"].append(stats(track(os.path.join(ca.REF, "wav", name + ".wav"))))
    rows["ours"].append(stats(track(os.path.join(ca.OUT, name + "_ours.wav"))))
for label, values in rows.items():
    values = numpy.array(values)
    print("%-4s median F0 %.0f Hz   range p10-p90 %.1f semitones   movement %.1f semitones/s" % (label, *numpy.median(values, axis=0)))
