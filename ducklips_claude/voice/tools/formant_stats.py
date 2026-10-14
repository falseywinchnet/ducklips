"""F1/F2 movement statistics (LPC, 10 ms) for ours vs ARCTIC."""
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
from analyze import lpc


def formant_track(path):
    data, rate = ca.read_wave(path)
    size = int(0.025 * rate)
    peak = numpy.abs(data).max()
    track = []
    for start in range(0, len(data) - size, int(0.01 * rate)):
        seg = data[start:start + size]
        if numpy.sqrt(numpy.mean(seg ** 2)) < peak * 0.08:
            track.append((numpy.nan, numpy.nan))
            continue
        seg = numpy.append(seg[0], seg[1:] - 0.94 * seg[:-1]) * numpy.hamming(size)
        roots = numpy.roots(lpc(seg, 14))
        found = []
        for root in roots:
            if root.imag <= 0:
                continue
            f = numpy.angle(root) * rate / (2 * numpy.pi)
            b = -numpy.log(abs(root)) * rate / numpy.pi
            if 200 < f < 3500 and b < 400:
                found.append(f)
        found.sort()
        track.append((found[0], found[1]) if len(found) >= 2 else (numpy.nan, numpy.nan))
    return numpy.array(track)


def movement(track):
    d = numpy.abs(numpy.diff(track, axis=0)) * 100  # Hz per second
    out = []
    for k in range(2):
        v = d[:, k][~numpy.isnan(d[:, k])]
        v = v[v < 8000]
        out.append(numpy.median(v))
    return out


rows = {"ref": [], "ours": []}
for number in range(1, 21):
    name = "arctic_a%04d" % number
    rows["ref"].append(movement(formant_track(os.path.join(ca.REF, "wav", name + ".wav"))))
    rows["ours"].append(movement(formant_track(os.path.join(ca.OUT, name + "_ours.wav"))))
for label, values in rows.items():
    print("%-4s median |dF1/dt| %.0f Hz/s   |dF2/dt| %.0f Hz/s" % (label, *numpy.median(numpy.array(values), axis=0)))
