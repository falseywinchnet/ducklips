"""Median F1-F3 per vowel in running speech: ours vs human (ARCTIC a0001-60)."""
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
from invert import lpc_formants

def survey(side):
    acc = {}
    for number in range(1, 61):
        name = "arctic_a%04d" % number
        if side == "ref":
            data, rate = ca.read_wave(os.path.join(ca.REF, "wav", name + ".wav"))
            phones = ca.reference_phones(name)
        else:
            data, rate = ca.read_wave(os.path.join(ca.OUT, name + "_ours.wav"))
            phones = ca.our_phones(os.path.join(ca.OUT, name + ".score"))
        for phone, start, end in phones:
            if phone not in ca.VOWELS or end - start < 0.04:
                continue
            values = []
            for t in numpy.linspace(start + 0.35 * (end - start), end - 0.35 * (end - start), 3):
                seg = data[int((t - 0.0125) * rate):int((t + 0.0125) * rate)]
                if len(seg) == int(0.025 * rate):
                    f = lpc_formants(seg, rate)
                    if f:
                        values.append(f)
            if values:
                acc.setdefault(phone, []).append(numpy.median(values, axis=0))
    return {p: numpy.median(v, axis=0) for p, v in acc.items() if len(v) >= 5}

ref, ours = survey("ref"), survey("ours")
print("vowel   human F1 F2 F3        ours F1 F2 F3")
for p in sorted(ref):
    if p in ours:
        print("%-4s  %5.0f %5.0f %5.0f    %5.0f %5.0f %5.0f" % (p, *ref[p], *ours[p]))
