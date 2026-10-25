"""Mouth postures from sound: the reversible core's second half.

A posture (the articulatory model's jaw, tongue body, tip, lips, larynx)
gives the tube (geometry_probe with the VocalTractLab-style tongue, cut into
the engine's nodes). The sound at the lips, peeled back through a candidate
tube (reverse_core.peel), gives the flow at the glottis that tube implies.
Only the right tube leaves that flow still while the folds are shut: a
wrong one leaves its uncancelled resonances ringing through every closed
phase (choosing among the eleven Story 2008 vowels by that alone recovered
11 of 11; by the peeled flow's energy, 6 of 11). From a hinted posture, the
search moves the organs until the closed phases fall still.

Known-answer test: a posture is rendered forward (sustain, AREAFILE), its
sound alone is carried back to a posture, and the two are compared.
usage: posture_recover.py TRUE_POSTURE HINT_POSTURE   (key=value lists)"""
import os
import subprocess
import sys
import tempfile

import numpy
from scipy.optimize import differential_evolution, minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reverse_core

ROOT = os.path.expanduser("~/ducklips/tract")
PROBE = os.path.join(ROOT, "build_dev", "geometry_probe")
SUSTAIN = os.path.join(ROOT, "build_dev", "sustain")
NODE = reverse_core.SOUND_SPEED / reverse_core.STEP_RATE
VISCOUS = 1.0   # the engine's viscous loss
# The pitch ear's gain in every render (sustain PITCHFB; the engine's
# default 0.08): off. It corrects the folds' tension once per closure, the
# closures timed in whole steps, so the equation error jumped as a
# hypothesis moved (the stalled joints: along /aa m aa/'s bad column the
# jump vanished with it off). Truths with it off: /aa b aa/ -242.3, Echo
# -245.6, /aa m aa/ -160.6 dB. Env PITCH_EAR restores a gain.
PITCH_EAR = float(os.environ.get("PITCH_EAR", "0"))
HOP = 0.02      # its per-hop low-pass (the articulatory body's)
JET_CLOSING = 0.04   # cm^2: with a seal, the reversible jet fades out below this (sustain REVJETCLOSE)
NOSE_POLE = 0.8   # the reversible nose's hop filter: a shelf (Tube::set_wave_shelf),
NOSE_ZERO = 0.79  # -0.32 dB a hop at 5 kHz, -0.47 at the half step rate
FIELDS = ("jaw", "tx", "ty", "tipx", "tipy", "ap", "pro", "lw", "drop")
# How far each field may move from the hint, and its scale for the search.
SPAN = {"jaw": 8.0, "tx": 1.2, "ty": 1.0, "tipx": 0.8, "tipy": 0.8, "ap": 1.5, "pro": 0.8, "lw": 2.0, "drop": 1.0}
FLOOR = {"ap": 0.02, "pro": 0.0, "lw": 0.5, "drop": 0.0, "jaw": 0.0}


def parse(text):
    return {k: float(v) for k, v in (item.split("=") for item in text.split())}


def profile(posture):
    """The articulatory model's cuts for a posture: arc positions, areas, and
    the tract's length (geometry_probe, the VocalTractLab-style tongue)."""
    line = " ".join("%s=%.6f" % (k, posture[k]) for k in FIELDS) + " anat.vtl=1\n"
    out = subprocess.run([PROBE], input=line, capture_output=True, text=True, cwd=ROOT).stdout.split("\n")
    length = float(out[0].split()[1])
    arcs = []
    areas = []
    for row in out[1:]:
        parts = row.split()
        if parts and parts[0] == "cut":
            arcs.append(float(parts[1]))
            areas.append(float(parts[2]))
    return numpy.array(arcs), numpy.array(areas), length


def nodes(posture, count=None):
    """The tube as the engine's nodes: whole nodes, each the narrowest of the
    cut profile (linear between cuts) over its span, so a closure is never
    stepped over and a cut sliding along moves the area continuously; and
    the length past the last whole node (cm), which the open end carries as
    air mass (Tube::set_active_length), so the length moves continuously
    too. (Rounded to whole nodes, the length stepped 0.36 cm at a time and
    the search's landscape stepped with it.) With `count`, that many nodes
    whatever the length (the profile's last area continued past its end;
    the length held within the engine's 0 to 2.5 nodes of air mass): the
    tube then moves continuously with the posture. (Free, it jumped as the
    length crossed a whole node: 0.05 spans from the true /ae/, a start
    rounding to 43 nodes where the truth has 44 never came back, +3.0 dB.)"""
    arcs, areas, length = profile(posture)
    if count is None:
        count = max(int(numpy.floor(length / NODE)), 2)
    else:
        length = min(max(length, (count + 1e-6) * NODE), (count + 2.5) * NODE)
    out = numpy.zeros(count)
    for node in range(count):
        low = node * NODE
        high = (node + 1) * NODE
        inside = (arcs > low) & (arcs < high)
        values = [numpy.interp(low, arcs, areas), numpy.interp(high, arcs, areas), numpy.interp(0.5 * (low + high), arcs, areas)]
        values.extend(areas[inside])
        out[node] = max(min(values), 1e-4)
    return out, length


def render(tube_and_length, seconds=0.3, pitch=120.0, second=None, glide=None, forced=None, keyframes=None,
           wave=None, jet=0.0, seal=0.0, velum=None, nose_dump=None, port_waves=None, folds=None, observed=None,
           jet_gate=None):
    """The engine's sound for a fixed tube as the reversible engine makes it
    (sustain: local steady linings, the air on the rest geometry, viscous
    and hop losses, no side branches or noise): the flow entering and the
    flow radiated, per step. With `second` (areas, the same count) and
    `glide` (start, duration in s), the tube glides to it (AREAFILE2).
    With `forced` (per step), the folds are driven by those returning waves
    and the first column is the flow they pass (FORCEWAVE). With
    `keyframes` (times, tubes), the tube moves through them (AREASEQ). With
    `wave`, sustain's wave file is kept there; `jet` the reversible jet's
    gain (REVJET); `seal` the area at which a node seals (SEAL); `velum`
    (times, areas) the port's keyframes, with the reversible nose; `folds`
    (times, postures as lower rest, upper rest, chink, move seconds) the
    folds' posture keyframes (FOLDSEQ); `jet_gate` (times, gates) the
    jet's gate keyframes (JETSEQ); `observed` (flows per step, NaN for
    none; gain) the forced folds observed against those flows (FORCEFLOW,
    OBSERVE)."""
    tube, length = tube_and_length
    # Each call leaves nothing behind: 36236 leftover folders filled the
    # Mini's disk during the first multi-start search.
    with tempfile.TemporaryDirectory(dir="/tmp/es") as folder:
        area_file = os.path.join(folder, "tube.txt")
        raw = os.path.join(folder, "tube.raw")
        with open(area_file, "w") as handle:
            handle.write("%.17g\n" % NODE)
            for area in tube:
                handle.write("%.17g\n" % area)
        if second is not None:
            with open(os.path.join(folder, "second.txt"), "w") as handle:
                handle.write("%.17g\n" % NODE)
                for area in second:
                    handle.write("%.17g\n" % area)
        environment = dict(os.environ, BRANCHES="", ASPIRATION="0", FRICATION="0", FLUX="0", FOLDTURB="0", CLOSING="0",
                           SKIN="0", AXIAL="0", WALLVAR="0", STATICAREA="1", VISCOUS="%g" % VISCOUS, HOPLOSS="%g" % HOP,
                           AREAFILE=area_file, RAW=raw, RAWDOUBLE="1", ACTIVELENGTH="%.17g" % length,
                           PITCHFB="%.17g" % PITCH_EAR)
        if keyframes is not None:
            with open(os.path.join(folder, "sequence.txt"), "w") as handle:
                for time, frame in zip(keyframes[0], keyframes[1]):
                    handle.write("%.17g " % time + " ".join("%.17g" % area for area in frame) + "\n")
            environment["AREASEQ"] = os.path.join(folder, "sequence.txt")
        if jet > 0.0:
            environment["REVJET"] = "%.17g" % jet
            if jet_gate is not None:
                with open(os.path.join(folder, "jet.txt"), "w") as handle:
                    for time, gate in zip(jet_gate[0], jet_gate[1]):
                        handle.write("%.17g %.17g\n" % (time, gate))
                environment["JETSEQ"] = os.path.join(folder, "jet.txt")
            if seal > 0.0:
                environment["REVJETCLOSE"] = "%.17g" % JET_CLOSING
        if seal > 0.0:
            environment["SEAL"] = "%.17g" % seal
        if velum is not None:
            # The reversible nose (REVNOSE; its hop a shelf) and the velum's
            # port over time.
            with open(os.path.join(folder, "velum.txt"), "w") as handle:
                for time, area in zip(velum[0], velum[1]):
                    handle.write("%.17g %.17g\n" % (time, area))
            environment.update(VELUMSEQ=os.path.join(folder, "velum.txt"), REVNOSE="1", NOSEHOP="%.17g" % NOSE_POLE,
                               NOSESHELF="%.17g" % NOSE_ZERO)
            if nose_dump is not None:
                environment["NOSEDUMP"] = nose_dump
            if port_waves is not None:
                environment["PORTWAVES"] = port_waves
        if folds is not None:
            # Three values a keyframe: smooth keyframes (FOLDKEYS); five
            # (with a move time): triggers (FOLDSEQ).
            with open(os.path.join(folder, "folds.txt"), "w") as handle:
                for time, posture in zip(folds[0], folds[1]):
                    handle.write("%.17g " % time + " ".join("%.17g" % value for value in posture) + "\n")
            smooth = len(folds[1]) > 0 and len(folds[1][0]) == 3
            environment["FOLDKEYS" if smooth else "FOLDSEQ"] = os.path.join(folder, "folds.txt")
        if forced is not None:
            numpy.asarray(forced, numpy.float64).tofile(os.path.join(folder, "forced.raw"))
            environment["FORCEWAVE"] = os.path.join(folder, "forced.raw")
            if observed is not None and observed[1] > 0.0:
                numpy.asarray(observed[0], numpy.float64).tofile(os.path.join(folder, "observed.raw"))
                environment.update(FORCEFLOW=os.path.join(folder, "observed.raw"), OBSERVE="%.17g" % observed[1])
        if second is not None:
            environment["AREAFILE2"] = os.path.join(folder, "second.txt")
            environment["GLIDE"] = "%.17g,%.17g" % glide
        subprocess.run([SUSTAIN, "aa", wave if wave is not None else os.path.join(folder, "tube.wav"), "%g" % pitch,
                        "800", "%g" % seconds],
                       env=environment, capture_output=True, cwd=ROOT)
        pairs = numpy.fromfile(raw, dtype=numpy.float64).reshape(-1, 2)
    return pairs[:, 0], pairs[:, 1]


def period_of(flow):
    """The voice's period (steps) from the radiated flow's autocorrelation."""
    x = flow - flow.mean()
    spectrum = numpy.fft.rfft(x, 2 * len(x))
    correlation = numpy.fft.irfft(numpy.abs(spectrum) ** 2)[:len(x)]
    low, high = int(reverse_core.STEP_RATE / 400), int(reverse_core.STEP_RATE / 60)
    # The shortest lag within 90% of the strongest: twice the period
    # correlates as well and on 0.2 s it won (60 Hz for a 120 Hz voice).
    window = correlation[low:high]
    candidates = numpy.nonzero(window >= 0.9 * window.max())[0]
    first = candidates[0]
    stop = first
    while stop + 1 < len(window) and window[stop + 1] >= window[stop]:
        stop += 1
    return low + int(stop)


def ringing(glottal, period, share=0.35):
    """The peeled flow's ringing while the folds are shut: per cycle, the
    stillest stretch (`share` of the period with the lowest mean flow), its
    flow-derivative energy against the whole cycle's (dB)."""
    derivative = numpy.diff(glottal)
    window = int(share * period)
    kernel = numpy.ones(window) / window
    ring = 0.0
    total = 0.0
    for cycle in range(len(derivative) // period):
        flow = glottal[cycle * period:(cycle + 1) * period]
        slope = derivative[cycle * period:(cycle + 1) * period]
        start = int(numpy.argmin(numpy.convolve(flow, kernel, "valid")))
        ring += numpy.sum(slope[start:start + window] ** 2)
        total += numpy.sum(slope ** 2)
    return 10.0 * numpy.log10(ring / max(total, 1e-30))


def ripple(glottal, period, top=5000.0):
    """The peeled flow's spectral ripple: its harmonics' levels (dB) against
    a smooth tilt (a quadratic in log frequency), rms over harmonics up to
    `top` Hz. A source has a smooth envelope; a tube that is not the true
    one leaves its uncancelled formants as peaks and dips in it."""
    cycles = len(glottal) // period
    x = glottal[:cycles * period] - numpy.mean(glottal[:cycles * period])
    spectrum = numpy.abs(numpy.fft.rfft(x * numpy.hanning(len(x))))
    f0 = reverse_core.STEP_RATE / period
    harmonics = numpy.arange(1, int(top / f0) + 1)
    bins = numpy.round(harmonics * cycles).astype(int)
    level = 20.0 * numpy.log10(numpy.maximum(spectrum[bins], 1e-30))
    logf = numpy.log(harmonics * f0)
    fit = numpy.polyval(numpy.polyfit(logf, level, 2), logf)
    return float(numpy.sqrt(numpy.mean((level - fit) ** 2)))


def harmonic_levels(flow, period, top=5000.0):
    """A flow's harmonic levels (dB) up to `top` Hz over whole cycles."""
    cycles = len(flow) // period
    x = flow[:cycles * period] - numpy.mean(flow[:cycles * period])
    spectrum = numpy.abs(numpy.fft.rfft(x * numpy.hanning(len(x))))
    f0 = reverse_core.STEP_RATE / period
    bins = numpy.round(numpy.arange(1, int(top / f0) + 1) * cycles).astype(int)
    return 20.0 * numpy.log10(numpy.maximum(spectrum[bins], 1e-30))


def consistency(tube, flow_out, start, length):
    """How far the glottal flow a tube implies (the whole sound peeled back
    through it, from silence) is from the flow our folds make driving that
    same tube (the forward render), over [start, start + length): the
    mismatch's level against the driven flow's swing (dB; lower is closer).
    Only the true tube makes them the same flow. (By the flows' harmonic
    levels in dB the landscape was jagged, its weak harmonics swinging
    most; by waveform it falls steadily toward the true posture: +4 dB at
    the hint, -11 at 60% of the way, -25 at the truth with a clip peeled
    from mid-voice, the tube's state at the clip's start unknown.)"""
    areas, tract_length = tube
    implied = reverse_core.peel_tube(areas, flow_out, reverse_core.lining_constants(areas), viscous=VISCOUS, hop=HOP,
                                     end_extra=tract_length - len(areas) * NODE)[start:start + length]
    driven, _ = render(tube)
    driven = driven[start:start + length]
    return -reverse_core.snr(driven - driven.mean(), implied - implied.mean())


def to_posture(unit, hint):
    """A point of the search (each field -1..1 of its span about the hint)
    as a posture."""
    posture = dict(hint)
    for k, u in zip(FIELDS, unit):
        posture[k] = max(hint[k] + SPAN[k] * float(numpy.clip(u, -1.0, 1.0)), FLOOR.get(k, -1e9))
    return posture


PRIOR = 1.0  # dB of mismatch per unit of squared distance from the hint (the codeword's posture)


def consistency_cost(unit, hint, flow_out, start, length):
    """Consistency plus a weak pull toward the hinted posture: where two
    postures sound the same (the lips' opening against their width at one
    area), the codeword's own posture decides."""
    return consistency(nodes(to_posture(unit, hint)), flow_out, start, length) + PRIOR * float(numpy.sum(unit ** 2))


def recover(flow_out, hint, start, length, workers=10):
    """The posture whose tube our folds drive to the same glottal flow the
    sound implies, searched over each field's span about the hint (a global
    search: the landscape has many local minima, and a local one from the
    hint stopped at 4.3 dB). Returns the posture, its cost, evaluations."""
    result = differential_evolution(consistency_cost, [(-1.0, 1.0)] * len(FIELDS),
                                    args=(hint, flow_out, start, length), popsize=12, maxiter=40,
                                    tol=1e-3, seed=1, workers=workers, updating="deferred", polish=False)
    # Then down into the basin: the true posture's well is deep and narrow
    # (-73 dB at it), coarser than a population samples.
    polished = minimize(consistency_cost, result.x, args=(hint, flow_out, start, length), method="Powell",
                        options={"xtol": 1e-4, "ftol": 1e-3, "maxfev": 3000})
    return to_posture(polished.x, hint), polished.fun, result.nfev + polished.nfev


def compare(name, true, found, tube):
    for k in FIELDS:
        print("  %-5s true %7.3f  %s %7.3f" % (k, true[k], name, found[k]))
    found_tube = nodes(found)
    size = min(len(tube[0]), len(found_tube[0]))
    print("  tube: true %.2f cm, %s %.2f cm; area rms difference %.3f cm^2 (true mean %.2f)" % (
        tube[1], name, found_tube[1], numpy.sqrt(numpy.mean((tube[0][:size] - found_tube[0][:size]) ** 2)),
        tube[0].mean()))
    return found_tube


def main():
    true = parse(sys.argv[1])
    hint = parse(sys.argv[2])
    tube = nodes(true)
    flow_in, flow_out = render(tube)
    start = int(0.1 * reverse_core.STEP_RATE)
    period = period_of(flow_out[start:])
    length = 15 * period
    print("true tube %.2f cm; period %d steps (%.1f Hz)" % (tube[1], period, reverse_core.STEP_RATE / period))
    print("consistency (mismatch, dB): true %.1f, hint %.1f" % (
        consistency(tube, flow_out, start, length), consistency(nodes(hint), flow_out, start, length)))
    found, value, evaluations = recover(flow_out, hint, start, length)
    print("recovered (%d evaluations, cost %.2f; mismatch %.1f dB):" % (
        evaluations, value, consistency(nodes(found), flow_out, start, length)))
    found_tube = compare("recovered", true, found, tube)
    # Forward again from the recovered posture: does the sound come back?
    _, again = render(found_tube)
    clip = flow_out[start:start + length]
    again = again[start:start + length]
    print("round trip: the recovered posture's sound against the true one's: %.1f dB (waveform), %.2f dB rms harmonic levels" % (
        reverse_core.snr(clip, again),
        numpy.sqrt(numpy.mean((harmonic_levels(clip, period) - harmonic_levels(again, period)) ** 2))))


if __name__ == "__main__":
    main()
