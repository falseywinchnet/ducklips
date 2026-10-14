"""Glottal pulse shape from a vowel_probe trace (area flow ps out per sample).
usage: glottal.py TRACE start_s dur_s"""
import sys

import numpy

d = numpy.loadtxt(sys.argv[1])
rate = 48000
seg = d[int(float(sys.argv[2]) * rate):int((float(sys.argv[2]) + float(sys.argv[3])) * rate)]
area, flow = seg[:, 0], seg[:, 1]
low, high = area.min(), area.max()
open_mask = area > low + 0.05 * (high - low)
# cycles from rising edges of the open mask
edges = numpy.where(numpy.diff(open_mask.astype(int)) == 1)[0]
oq, sq, closing = [], [], []
for start, end in zip(edges[:-1], edges[1:]):
    cycle = flow[start:end]
    mask = open_mask[start:end]
    opened = numpy.where(mask)[0]
    if len(opened) < 4:
        continue
    peak = opened[0] + numpy.argmax(cycle[opened])
    close = opened[-1]
    oq.append(len(opened) / (end - start))
    sq.append((peak - opened[0]) / max(close - peak, 1))
    derivative = numpy.diff(cycle)
    closing.append(-derivative.min() * rate)
print("cycles %d  F0 %.1f  OQ %.2f  speed quotient %.2f  peak closing rate %.3g cm3/s2  flow mean %.0f peak %.0f min %.0f  area max %.3f" % (
    len(oq), rate / numpy.mean(numpy.diff(edges)), numpy.mean(oq), numpy.mean(sq), numpy.mean(closing),
    flow.mean(), flow.max(), flow.min(), high))
