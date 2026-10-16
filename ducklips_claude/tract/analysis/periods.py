"""Glottal periods to a fraction of a microsecond: the cord flow from
speak WAVEDUMP (one line per 96 kHz step), each rise through half its
range timed by linear interpolation. usage: periods.py DUMP.txt"""
import sys, numpy
t, f = [], []
for line in open(sys.argv[1]):
    if line.startswith("wave"):
        parts = line.split()
        t.append(float(parts[1])); f.append(float(parts[-1]))
t = numpy.array(t); f = numpy.array(f)
level = 0.5 * (f.max() + f.min())
rise = numpy.where((f[:-1] < level) & (f[1:] >= level))[0]
times = t[rise] + (level - f[rise]) / (f[rise + 1] - f[rise]) * (t[rise + 1] - t[rise])
T = numpy.diff(times) * 1e6
rap = numpy.mean(numpy.abs(T[1:-1] - (T[:-2] + T[1:-1] + T[2:]) / 3))
print("cycles %d  F0 %.2f Hz  period %.1f us  mean |dT| %.2f us  RAP %.2f us (%.4f%%)  sd %.2f us" % (
    len(T), 1e6 / T.mean(), T.mean(), numpy.mean(numpy.abs(numpy.diff(T))), rap, 100 * rap / T.mean(), T.std()))
