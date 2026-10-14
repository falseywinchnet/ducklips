"""Formants of an area function (lines "areas a0 a1 ..." on stdin, glottis
first, 0.3646 cm sections): lossy plane-wave chain matrix with piston
radiation, closed glottis, peaks of the volume-velocity transfer."""
import sys

import numpy

SECTION = 0.3646
C = 35000.0
RHO_C = 1.14e-3 * C


def transfer(areas, hz):
    omega = 2 * numpy.pi * hz
    lip = max(areas[-1], 0.02)
    radius = numpy.sqrt(lip / numpy.pi)
    ka = omega / C * radius
    pressure = RHO_C / lip * (0.5 * ka * ka + 1j * 8 * ka / (3 * numpy.pi))
    flow = numpy.ones_like(hz, dtype=complex)
    k = omega / C - 1j * (1.2e-5 * numpy.sqrt(hz) + 0.25 / hz)
    for area in areas[::-1]:
        z = RHO_C / max(area, 0.005)
        c, s = numpy.cos(k * SECTION), numpy.sin(k * SECTION)
        pressure, flow = c * pressure + 1j * z * s * flow, 1j * s / z * pressure + c * flow
    return -20 * numpy.log10(numpy.abs(flow))


def formants(areas):
    hz = numpy.arange(150.0, 5000.0, 5.0)
    level = transfer(numpy.array(areas), hz)
    peaks = [hz[i] for i in range(1, len(hz) - 1) if level[i] > level[i - 1] and level[i] >= level[i + 1]]
    return peaks[:4]


for line in sys.stdin:
    if line.startswith("areas"):
        values = [float(v) for v in line.split()[1:]]
        print(" ".join("%5.0f" % f for f in formants(values)), "  min %.2f at %d" % (min(values), int(numpy.argmin(values))))
