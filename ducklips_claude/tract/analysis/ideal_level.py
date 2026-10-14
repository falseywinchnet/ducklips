"""Ideal output level of a held vowel from the sheet's own cross-sections
(sheet_areas rows above the cord), against another: a plane-wave chain with
modest losses (area_formants.py's), piston radiation, a glottal flow
spectrum falling 12 dB/octave above 150 Hz, harmonics of 120 Hz to 5 kHz,
radiated pressure ~ f U_lips. Prints F1-F3, the level, and per-band levels.
usage: ideal_level.py AREAS.txt [AREAS.txt...]"""
import sys, numpy

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
    return 1.0 / numpy.abs(flow)


def main():
    for path in sys.argv[1:]:
        rows = [l.split() for l in open(path)]
        start = max(i for i, r in enumerate(rows) if r[2] == "05_Vocal_cord") + 1
        areas = numpy.array([float(r[4]) for r in rows[start:]])
        hz = numpy.arange(100.0, 5000.0, 5.0)
        gain = transfer(areas, hz)
        db = 20 * numpy.log10(gain)
        peaks = [hz[i] for i in range(1, len(hz) - 1) if db[i] > db[i - 1] and db[i] >= db[i + 1]][:3]
        harmonics = numpy.arange(120.0, 5000.0, 120.0)
        g = transfer(areas, harmonics)
        source = 1.0 / (1.0 + (harmonics / 150.0) ** 2)
        p = (source * g * harmonics) ** 2
        bands = [(0, 300), (300, 600), (600, 1000), (1000, 2000), (2000, 5000)]
        print("%-20s nodes %d length %.1f cm  F %s  level %.1f dB  bands %s" % (
            path.split("/")[-1], len(areas), len(areas) * SECTION, " ".join("%4.0f" % f for f in peaks),
            10 * numpy.log10(p.sum()),
            " ".join("%.0f" % (10 * numpy.log10(p[(harmonics >= a) & (harmonics < b)].sum() + 1e-30)) for a, b in bands)))


if __name__ == "__main__":
    main()
