"""Ours vs ARCTIC around several occurrences of one phone, in three views.
usage: phone_view.py PHONE OUT.png [max_count]"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca

phone = sys.argv[1]
count = int(sys.argv[3]) if len(sys.argv) > 3 else 4
found = []
for number in range(1, 21):
    name = "arctic_a%04d" % number
    ref = ca.reference_phones(name)
    ours = ca.our_phones(os.path.join(ca.OUT, name + ".score"))
    for i, j in ca.align(ref, ours):
        if ref[i][0] == phone:
            found.append((name, ref[i], ours[j]))
found = found[:count]
views = ["stft", "cq", "causal"]
figure, axes = plt.subplots(len(views), 2 * len(found), figsize=(3.2 * 2 * len(found), 8), squeeze=False)
for k, (name, rp, op) in enumerate(found):
    for side, (label, item) in enumerate((("ref", rp), ("ours", op))):
        for row, view in enumerate(views):
            prefix = os.path.join(ca.OUT, name + "_" + label)
            field = numpy.load("%s_%s_enhanced.npy" % (prefix, view))
            hz = numpy.load("%s_%s_enhanced_hz.npy" % (prefix, view))[0]
            seconds = numpy.load("%s_seconds.npy" % prefix)[0]
            mask = (seconds > item[1] - 0.12) & (seconds < item[2] + 0.12)
            level = 20 * numpy.log10(field[:, mask] + 1e-7)
            top = numpy.percentile(20 * numpy.log10(field + 1e-7), 99.5)
            axis = axes[row][2 * k + side]
            axis.pcolormesh(seconds[mask] - item[1], hz, level, vmin=top - 60, vmax=top, cmap="magma", shading="auto")
            axis.axvline(0, color="c", lw=0.6)
            axis.axvline(item[2] - item[1], color="c", lw=0.6)
            axis.set_yscale("log")
            axis.set_ylim(100, 7600)
            axis.set_title("%s %s %s /%s/" % (name[-4:], label, view, phone), fontsize=7)
figure.tight_layout()
figure.savefig(sys.argv[2], dpi=70)
