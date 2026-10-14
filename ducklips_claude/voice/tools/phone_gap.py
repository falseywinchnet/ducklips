"""Per-phone fingerprint gap: our centroid minus the human centroid, by time
slice and band, for the phones the test confuses most.
usage: phone_gap.py PHONE [PHONE ...]"""
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
import phone_test as pt

names = ["arctic_a%04d" % n for n in range(1, 61)]
ref = pt.collect(names, "ref", ["cq"])
ours = pt.collect(names, "ours", ["cq"])
bands = len(pt.BANDS)
centres = [ca.GRID[a:b].mean() for a, b in pt.BANDS]
slices = ["-30ms", "start", "25%", "50%", "75%", "end", "+30ms"]
shown = list(range(0, bands, 3))
for phone in sys.argv[1:]:
    r = numpy.array([s[1] for s in ref if s[0] == phone]).reshape(-1, 7, bands)
    o = numpy.array([s[1] for s in ours if s[0] == phone]).reshape(-1, 7, bands)
    if len(r) == 0 or len(o) == 0:
        continue
    gap = o.mean(axis=0) - r.mean(axis=0)
    print("\n/%s/ (human n=%d, ours n=%d) ours-human dB; columns Hz" % (phone, len(r), len(o)))
    print("%7s" % "" + "".join("%6d" % centres[i] for i in shown))
    for k, label in enumerate(slices):
        print("%7s" % label + "".join("%+6.0f" % gap[k, i] for i in shown))
    print("%7s" % "human" + "".join("%6.0f" % r.mean(axis=0)[3, i] for i in shown) + "   (mid slice, dB re loud)")
