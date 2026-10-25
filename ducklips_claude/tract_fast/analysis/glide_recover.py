"""A moving mouth from sound: two tubes and the glide between them.

The engine's frame motion (each node's rest radius following its target,
the air the moving walls push, the linings recompiled, the viscous law
gliding) depends only on the targets, so for a hypothesis (tube A, tube B,
the glide's start and duration) reverse_core.moving_geometry gives the
geometry of every step and reverse_core.peel_moving carries the sound back
through it to the glottal flow it implies. The residual, as for one tube, is
the flow our folds make driving that moving tube (sustain GLIDE) less the
implied flow; the true hypothesis makes them one flow.

Stages: tube A from the stretch before the glide (tube_recover.recover,
hinted by the first phone's posture), tube B likewise from the stretch after
it, the timing on a grid, then Levenberg-Marquardt over everything.

Known-answer test: profile postures for two phones, rendered gliding; the
sound alone carried back; the recovered hypothesis rendered again.
usage: glide_recover.py FIRST SECOND   (phone names in the bdl5 profile)"""
import os
import sys
from multiprocessing import Pool

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import posture_recover
import reverse_core
import tube_recover

STEP_RATE = reverse_core.STEP_RATE
DATA = os.path.join(posture_recover.ROOT, "data")
SECONDS = 0.35          # the known answer's length
TRUE_GLIDE = (0.12, 0.10)
FIRST_WINDOW = (0.03, 0.115)     # before the glide (s)
SECOND_WINDOW = (0.235, 0.295)   # after it, within a 0.3 s static render
WHOLE_WINDOW = (0.03, 0.34)


def table(path):
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        name, rest = line.split(None, 1)
        out[name] = posture_recover.parse(rest.split("#")[0])
    return out


def resampled(areas, count):
    return numpy.interp(numpy.linspace(0.0, 1.0, count), numpy.linspace(0.0, 1.0, len(areas)), areas)


class GlideShape:
    """A hypothesis's parameters: tube A's and tube B's node log areas (the
    same count), the glide's start and duration (s); the tract's length
    fixed."""

    def __init__(self, count, length):
        self.count = count
        self.length = length

    def split(self, parameters):
        first = numpy.exp(numpy.clip(parameters[:self.count], numpy.log(0.02), numpy.log(30.0)))
        second = numpy.exp(numpy.clip(parameters[self.count:2 * self.count], numpy.log(0.02), numpy.log(30.0)))
        start = float(numpy.clip(parameters[-2], 0.0, 0.3))
        duration = float(numpy.clip(parameters[-1], 0.005, 0.3))
        return first, second, start, duration

    def residual(self, parameters, target):
        first, second, start, duration = self.split(parameters)
        return glide_residual(first, second, self.length, start, duration, target)


def glide_residual(first, second, length, start, duration, target):
    """For a gliding hypothesis, over the target's window and scale: our
    folds driven by the tract the sound implies (the equation error, as
    tube_recover.residual_of_equation; the first node's area followed step
    by step) less the implied flow; or, for a target without `equation`,
    our folds running free less the implied flow."""
    schedule = reverse_core.glide_schedule(first, second, len(target.flow_out), start, duration)
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    implied, pressure = reverse_core.peel_moving_glottis(geometry, target.flow_out, hop=posture_recover.HOP,
                                                         end_extra=length - len(first) * posture_recover.NODE)
    if target.equation:
        returning = 0.5 * (pressure - implied * reverse_core.RHO_C / geometry[:, 0, 0])
        driven, _ = posture_recover.render((first, length), target.seconds, second=second, glide=(start, duration),
                                           forced=returning)
    else:
        driven, _ = posture_recover.render((first, length), target.seconds, second=second, glide=(start, duration))
    size = min(len(driven), len(implied))
    difference = driven[:size] - implied[:size]
    if target.filter is not None:
        difference = tube_recover.sosfiltfilt(target.filter, difference)
    return tube_recover.fixed_clip(difference, target)


def level(residual, reference):
    """A residual's level against the reference flow's swing (dB)."""
    return 10.0 * numpy.log10(numpy.sum(residual ** 2) / numpy.sum((reference - reference.mean()) ** 2))


def window(bounds):
    return int(bounds[0] * STEP_RATE), int((bounds[1] - bounds[0]) * STEP_RATE)


def timing_cost(task):
    first, second, length, start, duration, target = task
    return float(numpy.sum(glide_residual(first, second, length, start, duration, target) ** 2))


def timing_grid(pool, first, second, length, starts, durations, target):
    """The best (start, duration) of a grid; returns them and the grid's
    shape for the log."""
    tasks = []
    for glide_start in starts:
        for duration in durations:
            if duration > 0.004:
                tasks.append((first, second, length, float(glide_start), float(duration), target))
    costs = numpy.array(pool.map(timing_cost, tasks))
    best = tasks[int(numpy.argmin(costs))]
    print("timing grid (%d, band %.0f Hz): start %.4f duration %.4f" % (len(tasks), target.band, best[3], best[4]))
    sys.stdout.flush()
    return best[3], best[4], len(tasks), target.band


def main():
    profile = table(os.path.join(DATA, "profile_bdl5.organs"))
    hints = table(os.path.join(DATA, "organs_vtl_tongue.txt"))
    first_name, second_name = sys.argv[1], sys.argv[2]
    true_first, length = posture_recover.nodes(profile[first_name])
    count = len(true_first)
    # The second tube as the posture model makes it at the first's count
    # (the moving tube keeps one count and one length).
    true_second = posture_recover.nodes(profile[second_name], count)[0]
    flow_in, flow_out = posture_recover.render((true_first, length), SECONDS, second=true_second, glide=TRUE_GLIDE)
    start, size = window(WHOLE_WINDOW)
    whole = tube_recover.Target(flow_out, start, size, 0.0, 1.0, True)
    reference = flow_in[start:start + size]
    truth = numpy.concatenate([numpy.log(true_first), numpy.log(true_second), TRUE_GLIDE])
    shape = GlideShape(count, length)
    print("%s -> %s glide %.3f s + %.3f s, %d nodes, %.2f cm; the truth's residual %.1f dB" % (
        first_name, second_name, TRUE_GLIDE[0], TRUE_GLIDE[1], count, length,
        level(shape.residual(truth, whole), reference)))
    sys.stdout.flush()
    with Pool(10) as pool:
        # Tube A, from the stretch before the glide.
        begin, span = window(FIRST_WINDOW)
        value, found_first = tube_recover.recover(pool, hints[first_name], flow_out, begin, span, (true_first, length),
                                                  count)
        first = resampled(found_first[0], count)
        print("A: %.1f dB, tube rms %.3f" % (value, tube_recover.tube_rms((first, length), (true_first, length))))
        # Tube B, from the stretch after it (the peel from the start through
        # B alone is wrong before the glide; its error has died by then).
        begin, span = window(SECOND_WINDOW)
        value, found_second = tube_recover.recover(pool, hints[second_name], flow_out, begin, span,
                                                   (true_second, length), count)
        second = resampled(found_second[0], count)
        print("B: %.1f dB, tube rms %.3f" % (value, tube_recover.tube_rms((second, length), (true_second, length))))
        sys.stdout.flush()
        # The timing, coarse to fine: 1 ms off the true start costs the
        # residual 111 dB of its depth (-134 -> -22.8 dB), so a 10 ms grid
        # is read through the low band first, then 1 ms, then 0.2 ms.
        low = tube_recover.Target(flow_out, start, size, 1500.0, 1.0, True)
        best = timing_grid(pool, first, second, length, numpy.arange(0.06, 0.22, 0.01),
                           numpy.array([0.03, 0.06, 0.09, 0.12, 0.15, 0.18]), low)
        best = timing_grid(pool, first, second, length, best[0] + numpy.arange(-0.006, 0.0065, 0.002),
                           best[1] + numpy.arange(-0.015, 0.016, 0.005), whole)
        best = timing_grid(pool, first, second, length, best[0] + numpy.arange(-0.0012, 0.0013, 0.0004),
                           best[1] + numpy.arange(-0.003, 0.0031, 0.001), whole)
        parameters = numpy.concatenate([numpy.log(first), numpy.log(second), [best[0], best[1]]])
        print("before the joint search: %.1f dB" % level(shape.residual(parameters, whole), reference))
        sys.stdout.flush()
        for band in (1500.0, 0.0, 0.0):
            parameters, used = tube_recover.stage(pool, shape, parameters,
                                                  tube_recover.Target(flow_out, start, size, band, 1.0, True), 30)
            found = shape.split(parameters)
            print("joint, band %5.0f Hz: %d evaluations; %.1f dB; A rms %.3f, B rms %.3f; start %.4f duration %.4f" % (
                band, used, level(shape.residual(parameters, whole), reference),
                tube_recover.tube_rms((found[0], length), (true_first, length)),
                tube_recover.tube_rms((found[1], length), (true_second, length)), found[2], found[3]))
            sys.stdout.flush()
    first, second, glide_start, duration = shape.split(parameters)
    _, again = posture_recover.render((first, length), SECONDS, second=second, glide=(glide_start, duration))
    print("round trip: the recovered glide's sound against the true one's: %.1f dB (whole window)" % reverse_core.snr(
        flow_out[start:start + size], again[start:start + size]))


if __name__ == "__main__":
    main()
