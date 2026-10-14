"""Long-term average spectrum (voiced/active frames), ours vs human, all 60."""
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca

edges = [80, 160, 315, 630, 1250, 2500, 5000, 8000]
totals = {"ref": numpy.zeros(len(edges) - 1), "ours": numpy.zeros(len(edges) - 1)}
for number in range(1, 61):
    name = "arctic_a%04d" % number
    for label, path in (("ref", os.path.join(ca.REF, "wav", name + ".wav")), ("ours", os.path.join(ca.OUT, name + "_ours.wav"))):
        data, rate = ca.read_wave(path)
        size = 1024
        frames = [data[i:i + size] for i in range(0, len(data) - size, size // 2)]
        energy = numpy.array([numpy.mean(f ** 2) for f in frames])
        active = [f for f, e in zip(frames, energy) if e > energy.max() * 1e-3]
        spectrum = numpy.mean([numpy.abs(numpy.fft.rfft(f * numpy.hanning(size))) ** 2 for f in active], axis=0)
        f = numpy.fft.rfftfreq(size, 1 / rate)
        totals[label] += [spectrum[(f >= a) & (f < b)].sum() for a, b in zip(edges[:-1], edges[1:])]
for label, total in totals.items():
    db = 10 * numpy.log10(total / total.max())
    print("%-4s " % label + "  ".join("%d-%d:%+.1f" % (a, b, d) for a, b, d in zip(edges[:-1], edges[1:], db)))
