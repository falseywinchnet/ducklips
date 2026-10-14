"""Median F1-F3 of each vowel in ARCTIC bdl (middle 40% of each vowel,
LPC order 14 at 16 kHz, candidates with bandwidth < 400 Hz)."""
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
from analyze import lpc


def formants(seg, rate):
    seg = numpy.append(seg[0], seg[1:] - 0.94 * seg[:-1]) * numpy.hamming(len(seg))
    found = []
    for root in numpy.roots(lpc(seg, 14)):
        if root.imag <= 0:
            continue
        f = numpy.angle(root) * rate / (2 * numpy.pi)
        b = -numpy.log(abs(root)) * rate / numpy.pi
        if 180 < f < 4000 and b < 400:
            found.append(f)
    found.sort()
    return found[:3] if len(found) >= 3 else None


acc = {}
for number in range(1, 61):
    name = "arctic_a%04d" % number
    data, rate = ca.read_wave(os.path.join(ca.REF, "wav", name + ".wav"))
    for phone, start, end in ca.reference_phones(name):
        if phone not in ca.VOWELS or end - start < 0.05:
            continue
        frames = []
        for t in numpy.linspace(start + 0.3 * (end - start), end - 0.3 * (end - start), 4):
            seg = data[int((t - 0.0125) * rate):int((t + 0.0125) * rate)]
            f = formants(seg, rate)
            if f:
                frames.append(f)
        if frames:
            acc.setdefault(phone, []).append(numpy.median(numpy.array(frames), axis=0))
for phone, values in sorted(acc.items()):
    m = numpy.median(numpy.array(values), axis=0)
    print("%-3s n=%3d  F1 %4.0f  F2 %4.0f  F3 %4.0f" % (phone, len(values), *m))
