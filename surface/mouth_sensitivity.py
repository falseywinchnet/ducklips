"""What the mouth's shape does to its formants, and what shape a heard
formant pattern asks for.

For one area function (glottis to lips) this solves the plane-wave tube
(closed glottis, radiating lips, small wall loss) for F1-F3 and the standing
wave at each. A section where the wave's kinetic energy exceeds its
potential energy is a place where narrowing lowers that formant; where
potential energy dominates, narrowing raises it (Fant's sensitivity
functions):

    dF/F = sum_i S_i dA_i / A_i,   S_i = (KE_i - PE_i) / (KE + PE total)

Given target formants (for example, measured from a human recording), it
then finds the smoothest change of log area that moves the model's formants
there (minimum curvature, Newton steps on the sensitivity functions). The
answer is one shape among many that give the same three formants, the
smoothest; it says which regions must narrow or widen, not the exact
anatomy.

usage: python3 mouth_sensitivity.py AREAS.txt TIME [F1 F2 F3]
  AREAS.txt: rows of "time area_1 ... area_N" (cm^2), from surface_say --areas
"""
import sys

import numpy

DENSITY = 1.14
SOUND_SPEED = 350.0
LENGTH = 0.175
FIRST_SHAPED = 8  # sections below this are the fixed epilarynx


def read_row(path, time):
    best = None
    best_distance = 1.0e9
    for line in open(path):
        values = line.split("|")[0].split()
        distance = abs(float(values[0]) - time)
        if distance < best_distance:
            best_distance = distance
            best = values
    areas = numpy.zeros(len(best) - 1)
    for index in range(1, len(best)):
        areas[index - 1] = float(best[index]) * 1.0e-4
    return float(best[0]), areas


def hertz_list(values):
    text = ""
    for index in range(len(values)):
        if index > 0:
            text += " "
        text += "%.0f" % values[index]
    return text


def standing_wave(areas, frequency):
    """Pressure and volume velocity at every section (lips: U = 1), and the
    volume velocity left at the closed glottis."""
    count = len(areas)
    section = LENGTH / count
    omega = 2.0 * numpy.pi * frequency
    wavenumber = omega / SOUND_SPEED * (1.0 - 0.004j)
    lip = areas[-1]
    radius = numpy.sqrt(lip / numpy.pi)
    ka = omega / SOUND_SPEED * radius
    radiation = DENSITY * SOUND_SPEED / lip * (0.5 * ka * ka + 1j * 8.0 * ka / (3.0 * numpy.pi))
    pressure = numpy.zeros(count, dtype=complex)
    flow = numpy.zeros(count, dtype=complex)
    p = radiation * 1.0
    u = 1.0 + 0.0j
    for index in range(count - 1, -1, -1):
        impedance = DENSITY * SOUND_SPEED / areas[index]
        c = numpy.cos(wavenumber * section)
        s = numpy.sin(wavenumber * section)
        p_in = c * p + 1j * impedance * s * u
        u_in = 1j * s / impedance * p + c * u
        pressure[index] = 0.5 * (p + p_in)
        flow[index] = 0.5 * (u + u_in)
        p = p_in
        u = u_in
    return pressure, flow, u


def formants(areas, count=3):
    grid = numpy.arange(80.0, 5000.0, 2.0)
    magnitude = numpy.zeros(len(grid))
    for index in range(len(grid)):
        pressure, flow, glottal = standing_wave(areas, grid[index])
        magnitude[index] = abs(glottal)
    found = []
    for index in range(1, len(grid) - 1):
        if magnitude[index] < magnitude[index - 1] and magnitude[index] <= magnitude[index + 1]:
            found.append(grid[index])
            if len(found) == count:
                break
    return numpy.array(found)


def sensitivity(areas, frequency):
    pressure, flow, glottal = standing_wave(areas, frequency)
    section = LENGTH / len(areas)
    kinetic = DENSITY * section * numpy.abs(flow) ** 2 / areas
    potential = areas * section * numpy.abs(pressure) ** 2 / (DENSITY * SOUND_SPEED ** 2)
    total = kinetic.sum() + potential.sum()
    result = (kinetic - potential) / total
    return result


def invert(areas, targets, steps=8):
    """Smoothest log-area change (shaped sections only) that brings F1-F3
    to the targets."""
    count = len(areas)
    log_area = numpy.log(areas)
    shaped = numpy.arange(FIRST_SHAPED, count)
    curvature = numpy.zeros((len(shaped) - 2, len(shaped)))
    for row in range(len(shaped) - 2):
        curvature[row, row] = 1.0
        curvature[row, row + 1] = -2.0
        curvature[row, row + 2] = 1.0
    regular = curvature.T @ curvature + 1.0e-3 * numpy.eye(len(shaped))
    current = formants(numpy.exp(log_area))
    for step in range(steps):
        current = formants(numpy.exp(log_area))
        if len(current) < 3:
            break
        rows = numpy.zeros((3, len(shaped)))
        for index in range(3):
            rows[index] = sensitivity(numpy.exp(log_area), current[index])[shaped]
        error = numpy.log(targets / current[:3])
        # minimise x' R x subject to rows x = error
        inverse = numpy.linalg.solve(regular, rows.T)
        change = inverse @ numpy.linalg.solve(rows @ inverse, error)
        log_area[shaped] += 0.7 * change
    return numpy.exp(log_area), formants(numpy.exp(log_area))


def region_of(index, count):
    from_lips = (count - index - 0.5) * LENGTH / count * 100.0
    if from_lips < 1.0:
        return "lips"
    if from_lips < 2.5:
        return "teeth/alveolar"
    if from_lips < 5.5:
        return "hard palate"
    if from_lips < 8.0:
        return "soft palate"
    if from_lips < 11.0:
        return "upper pharynx"
    return "lower pharynx"


def describe(areas, label):
    count = len(areas)
    names = ["lower pharynx", "upper pharynx", "soft palate", "hard palate", "teeth/alveolar", "lips"]
    sums = {}
    counts = {}
    minimum = {}
    for name in names:
        sums[name] = 0.0
        counts[name] = 0
        minimum[name] = 1.0e9
    for index in range(FIRST_SHAPED, count):
        name = region_of(index, count)
        sums[name] += areas[index]
        counts[name] += 1
        minimum[name] = min(minimum[name], areas[index])
    parts = []
    for name in names:
        if counts[name] > 0:
            parts.append("%s %.2f (min %.2f)" % (name, sums[name] / counts[name] * 1.0e4, minimum[name] * 1.0e4))
    print("  %-8s mean area cm^2: %s" % (label, ", ".join(parts)))


def main():
    path = sys.argv[1]
    time, areas = read_row(path, float(sys.argv[2]))
    count = len(areas)
    found = formants(areas)
    print("t = %.3f s: model formants %s Hz" % (time, hertz_list(found)))
    describe(areas, "model")
    print("  sensitivity: + means narrowing here RAISES the formant, - lowers it")
    print("  %-16s %6s %6s %6s" % ("region", "F1", "F2", "F3"))
    names = ["lower pharynx", "upper pharynx", "soft palate", "hard palate", "teeth/alveolar", "lips"]
    curves = []
    for index in range(len(found)):
        curves.append(sensitivity(areas, found[index]))
    for name in names:
        values = []
        for index in range(len(found)):
            total = 0.0
            for section in range(FIRST_SHAPED, count):
                if region_of(section, count) == name:
                    total += curves[index][section]
            # narrowing (dA < 0) changes F by -S: report the effect of narrowing
            values.append(-total)
        print("  %-16s %+6.2f %+6.2f %+6.2f" % (name, values[0], values[1], values[2]))
    if len(sys.argv) >= 6:
        targets = numpy.array([float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])])
        shape, reached = invert(areas, targets)
        print("  target %s Hz -> reached %s Hz" % (hertz_list(targets), hertz_list(reached)))
        describe(shape, "needed")
        ratio = shape / areas
        print("  change by region (area ratio needed/model, mean and extreme):")
        for name in names:
            selected = []
            for section in range(FIRST_SHAPED, count):
                if region_of(section, count) == name:
                    selected.append(ratio[section])
            selected = numpy.array(selected)
            extreme = selected[numpy.argmax(numpy.abs(numpy.log(selected)))]
            print("    %-16s x%.2f (extreme x%.2f)" % (name, numpy.exp(numpy.log(selected).mean()), extreme))


main()
