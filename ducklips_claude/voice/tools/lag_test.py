"""Does our acoustic timing lag the planned phone times? Shift our phone
boundaries by a constant and rerun the identification test."""
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
import phone_test as pt

views = ["cq", "causal"]
train = pt.collect(["arctic_a%04d" % n for n in range(1, 41)], "ref", views)
known = {s[0] for s in train}
names = ["arctic_a%04d" % n for n in range(1, 61)]
fields = {}
for name in names:
    prefix = os.path.join(ca.OUT, name + "_ours")
    fields[name] = [pt.load_levels(prefix, view) for view in views]
for shift in (-0.04, -0.02, 0.0, 0.02, 0.04, 0.06):
    test = []
    for name in names:
        for phone, start, end in ca.our_phones(os.path.join(ca.OUT, name + ".score")):
            phone = {"dx": "t", "ax": "ah"}.get(phone, phone)
            test.append((phone, pt.fingerprint(fields[name], start + shift, end + shift), name))
    results = [r for r in pt.classify(train, test, 5) if r[0] in known]
    print("shift %+.2f s: accuracy %.1f%%" % (shift, 100.0 * numpy.mean([a == b for a, b in results])), flush=True)
