"""Cycle-to-cycle jitter and shimmer from a vowel_probe trace.
usage: jitter.py TRACE start dur"""
import sys

import numpy

d = numpy.loadtxt(sys.argv[1])
rate = 48000
seg = d[int(float(sys.argv[2]) * rate):int((float(sys.argv[2]) + float(sys.argv[3])) * rate)]
area = seg[:, 0]
out = seg[:, 3]
level = area.min() + 0.3 * (area.max() - area.min())
rising = numpy.where((area[:-1] < level) & (area[1:] >= level))[0]
# sub-sample crossing times
times = rising + (level - area[rising]) / (area[rising + 1] - area[rising])
periods = numpy.diff(times) / rate
amplitudes = numpy.array([numpy.abs(out[int(a):int(b)]).max() for a, b in zip(times[:-1], times[1:])])
jitter = numpy.mean(numpy.abs(numpy.diff(periods))) / numpy.mean(periods) * 100
shimmer = numpy.mean(numpy.abs(numpy.diff(amplitudes))) / numpy.mean(amplitudes) * 100
print("F0 %.1f Hz  jitter %.2f%%  shimmer %.2f%%  (human modal: jitter 0.2-1%%, shimmer 2-5%%)" % (1 / periods.mean(), jitter, shimmer))
