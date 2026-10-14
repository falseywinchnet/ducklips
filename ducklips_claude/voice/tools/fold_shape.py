"""How the folds move, judged by the shape of each cycle rather than by
spectra: sustained vowels from build/vowel_probe across pitch and adduction,
read from the simulator's own glottal area and flow (TRACE).

Per case, medians over cycles of the steady part (0.6-1.4 s):
  OQ       share of the cycle the glottis is open (area above the chink
           floor + 15% of the swing, so a rebound is not a cycle)
  SQa      area opening time / closing time (human folds open slowly and
           close fast: > 1)
  SQf      the same for the flow, which lags the area (tract inertance)
  snap     peak flow-decline rate / peak AC flow, 1/ms (how abruptly the
           flow is cut at closure: the main excitation)
  rebound  largest area reopening after closure, re the cycle's swing
  jitter   mean |period difference| / mean period, %
  flow     mean and AC peak, cm^3/s

Fold material overrides come from the environment as for vowel_probe
(F_LM, F_UK, THICK1, ...), so candidate folds can be compared by eye.
usage: fold_shape.py"""
import os
import subprocess

import numpy
from scipy.signal import find_peaks

RATE = 48000.0
# (label, tongue place, area, lips, pitch Hz, lungs Pa, rest gap cm)
CASES = (
    ("a 100 modal", 7.0, 0.724, 5.0, 100.0, 750.0, -0.008),
    ("a 125 modal", 7.0, 0.724, 5.0, 125.0, 800.0, -0.008),
    ("a 155 modal", 7.0, 0.724, 5.0, 155.0, 900.0, -0.008),
    ("i 125 modal", 13.75, 0.179, 4.9, 125.0, 800.0, -0.008),
    ("a 125 loose", 7.0, 0.724, 5.0, 125.0, 800.0, 0.01),
    ("a 125 press", 7.0, 0.724, 5.0, 125.0, 800.0, -0.02),
)
PROBE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "build", "vowel_probe")


def cycle_shapes(area, flow):
    floor = numpy.percentile(area, 2)
    swing = numpy.percentile(area, 98) - floor
    if swing < 0.005:
        return None
    opened = area > floor + 0.15 * swing
    starts = numpy.where(opened[1:] & ~opened[:-1])[0] + 1
    rows = []
    for first, last in zip(starts[:-1], starts[1:]):
        a = area[first:last]
        f = flow[first:last]
        open_part = numpy.where(opened[first:last])[0]
        if len(open_part) < 4:
            continue
        close = open_part[-1]
        peak_area = int(numpy.argmax(a[:close + 1]))
        peak_flow = int(numpy.argmax(f))
        # The flow's own open phase, bounded on each side of its peak where
        # it falls to 5% of its swing (the window ends where the next cycle's
        # area opens, and the lagging flow may already be rising there).
        flow_floor = f.min()
        low = f <= flow_floor + 0.05 * (f.max() - flow_floor)
        before = numpy.where(low[:peak_flow])[0]
        later = numpy.where(low[peak_flow:])[0]
        flow_open = before[-1] if len(before) else 0
        flow_close = peak_flow + later[0] if len(later) else len(f) - 1
        derivative = numpy.diff(f) * RATE / 1000.0  # cm^3/s per ms
        ac = f.max() - flow_floor
        # A rebound is a second hump in the area after closure, not the
        # next cycle's rise at the window's end.
        humps, properties = find_peaks(a[close:], prominence=0.01 * swing)
        rebound = properties["prominences"].max() / swing if len(humps) else 0.0
        rows.append((len(open_part) / len(a), max(peak_area, 1) / max(close - peak_area, 1),
                     max(peak_flow - flow_open, 1) / max(flow_close - peak_flow, 1),
                     -derivative.min() / max(ac, 1e-9), rebound, (last - first) / RATE, f.mean(), ac))
    return numpy.array(rows) if rows else None


def run_case(case):
    label, place, area, lips, pitch, lungs, gap = case
    trace = "/tmp/fold_shape_trace.txt"
    environment = dict(os.environ, TRACE=trace)
    subprocess.run([PROBE, "/tmp/fold_shape.wav", str(place), str(area), str(lips), str(pitch), str(lungs), str(gap),
                    "1.5"], env=environment, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    data = numpy.loadtxt(trace)[int(0.6 * RATE):int(1.4 * RATE)]
    if not numpy.isfinite(data).all():
        return label, None
    return label, cycle_shapes(data[:, 0], data[:, 1])


def main():
    print("case          F0   OQ   SQa  SQf  snap/ms rebound jitter%  mean  AC cm3/s")
    for case in CASES:
        label, rows = run_case(case)
        if rows is None or len(rows) < 4:
            print("%-12s  no steady vibration" % label)
            continue
        periods = rows[:, 5]
        jitter = 100.0 * numpy.mean(numpy.abs(numpy.diff(periods))) / periods.mean()
        median = numpy.median(rows, axis=0)
        print("%-12s %4.0f %5.2f %5.2f %4.2f %6.2f %7.2f %7.1f %5.0f %5.0f" % (
            label, 1.0 / median[5], median[0], median[1], median[2], median[3], median[4], jitter, median[6],
            median[7]))


if __name__ == "__main__":
    main()
