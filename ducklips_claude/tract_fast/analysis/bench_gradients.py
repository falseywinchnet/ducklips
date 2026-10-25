"""Measurements of the derivative toolchain (gradients.py, ad/): primal
equivalence with the reference pipeline, the exact Jacobian against finite
differences, and the cost of each.

The word: bdl5 profile vowels held and moved on sequence_recover's plan (the
known-answer test's truth), its sound rendered by sustain; the hypothesis
near the truth (log areas and times perturbed, seeded), where a search's
Jacobians are taken.

usage: bench_gradients.py check|jacobian|sweep|timing|descend|reference [PHONES]   (default AA,IY,UW)"""
import os
import resource
import sys
import time

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gradients
import posture_recover
import reverse_core
import sequence_recover
import tube_recover

DATA = os.path.join(posture_recover.ROOT, "data")
AREA_NOISE = 0.02     # rms of the hypothesis's log-area offsets from the truth
TIME_NOISE = 0.0003   # s, rms of its times' offsets


def word(phones_text):
    """The truth, its sound, the target and a hypothesis near the truth."""
    profile = sequence_recover.table(os.path.join(DATA, "profile_bdl5.organs"))
    names = phones_text.split(",")
    first, length = posture_recover.nodes(profile[names[0]])
    count = len(first)
    tubes = numpy.array([posture_recover.nodes(profile[name], count)[0] for name in names])
    arrivals, leavings = sequence_recover.plan(len(names))
    seconds = leavings[-1] + 0.05
    times, frames = sequence_recover.keyframes(tubes, arrivals, leavings)
    flow_in, flow_out, _ = gradients.sustain((tubes[0], length), seconds, keyframes=(times, frames))
    start, size = sequence_recover.window(0.02, seconds - 0.01)
    target = tube_recover.Target(flow_out, start, size, 0.0, 1.0, True)
    target.reference = flow_in[start:start + size]   # the glottal flow, sequence_recover.level's reference
    shape = sequence_recover.SequenceShape(len(names), count, length)
    truth = shape.join(tubes, arrivals, leavings)
    generator = numpy.random.default_rng(11)
    hypothesis = truth.copy()
    hypothesis[:len(names) * count] += AREA_NOISE * generator.standard_normal(len(names) * count)
    hypothesis[len(names) * count:] += TIME_NOISE * generator.standard_normal(shape.time_count)
    return shape, target, truth, hypothesis


def reference_residual(shape, parameters, target, keep=None):
    """sequence_recover.sequence_residual (no seal, nose or jet) with the
    sustain render (gradients.sustain: posture_recover.render's run)."""
    tubes, arrivals, leavings = shape.split(parameters)
    times, frames = sequence_recover.keyframes(tubes, arrivals, leavings)
    schedule = reverse_core.keyframe_schedule(times, frames, len(target.flow_out))
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    implied, pressure = reverse_core.peel_moving_glottis(geometry, target.flow_out, hop=posture_recover.HOP,
                                                         end_extra=shape.length - frames.shape[1] * posture_recover.NODE)
    forced = 0.5 * (pressure - implied * reverse_core.RHO_C / geometry[:, 0, 0])
    driven, _, dump = gradients.sustain((frames[0], shape.length), target.seconds, keyframes=(times, frames),
                                        forced=forced, geodump=keep)
    size = min(len(driven), len(implied))
    difference = driven[:size] - implied[:size]
    residual = tube_recover.fixed_clip(difference, target)
    if keep is not None:
        keep.update(implied=implied, pressure=pressure, forced=forced, driven=driven, geometry=geometry, dump=dump)
    return residual


def deviation(a, b):
    a = numpy.asarray(a, float)
    b = numpy.asarray(b, float)
    size = min(len(a), len(b))
    a = a[:size]
    b = b[:size]
    scale = max(numpy.max(numpy.abs(b)), 1e-300)
    identical = int(numpy.sum(a == b))
    return "max |diff| %.3g (relative to max %.3g), identical %d of %d" % (
        numpy.max(numpy.abs(a - b)), numpy.max(numpy.abs(a - b)) / scale, identical, size)


def check(phones_text):
    """The native copies' primal outputs against the reference's, stage by stage."""
    shape, target, truth, hypothesis = word(phones_text)
    model = gradients.SequenceJacobian(shape, target)
    keep = {}
    began = time.perf_counter()
    reference = reference_residual(shape, hypothesis, target, keep)
    reference_seconds = time.perf_counter() - began
    count, steps, schedule, end_extra, n = model._layout(hypothesis)
    print("word %s: %d nodes, %d sound steps, %d render steps, %d parameters" % (
        phones_text, n, count, steps, len(hypothesis)))
    native = gradients.native_geometry(schedule)
    geometry = keep["geometry"]
    for column, name in ((0, "area"), (1, "pushed air"), (2, "link resistance")):
        print("geometry %s, native vs reverse_core.moving_geometry: %s" % (
            name, deviation(native[:count, :, column].ravel(), geometry[:, :, column].ravel())))
    dump = keep["dump"]
    for column, name in ((0, "area"), (1, "pushed air"), (2, "link resistance")):
        print("geometry %s, reverse_core vs sustain GEODUMP: %s" % (
            name, deviation(geometry[:, :, column].ravel(), dump[:count, :, column].ravel())))
    implied, pressure, area0 = gradients.peel(schedule, count, target.flow_out, end_extra)
    print("peel implied flow, native vs numba: %s" % deviation(implied, keep["implied"]))
    print("peel first-node pressure, native vs numba: %s" % deviation(pressure, keep["pressure"]))
    # The folds alone: sustain's own forced waves and first-node area.
    waves = numpy.zeros(steps)
    shared = min(count, steps)
    waves[:shared] = keep["forced"][:shared]
    folds_sustain_area = gradients.forced_folds(waves, dump[:steps, 0, 0], n)
    print("forced folds (sustain's A0), native vs sustain FORCEWAVE: %s" % deviation(folds_sustain_area,
                                                                                    keep["driven"]))
    folds_own_area = gradients.forced_folds(waves, area0[:steps], n)
    print("forced folds (native A0), native vs sustain FORCEWAVE: %s" % deviation(folds_own_area, keep["driven"]))
    began = time.perf_counter()
    residual = model.residual(hypothesis)
    native_seconds = time.perf_counter() - began
    print("residual, native vs reference: %s; levels %.2f / %.2f dB" % (
        deviation(residual, reference), sequence_recover.level(residual, target.flow_out[target.start:][:target.length]),
        sequence_recover.level(reference, target.flow_out[target.start:][:target.length])))
    truth_residual = model.residual(truth)
    print("at the truth: native residual rms %.3g (sound rms %.3g)" % (
        numpy.sqrt(numpy.mean(truth_residual ** 2)), numpy.std(target.flow_out)))
    print("one residual: reference %.3f s, native %.3f s" % (reference_seconds, native_seconds))


def jacobian(phones_text, block=gradients.BLOCK):
    """The exact Jacobian against central differences of the same residual."""
    shape, target, truth, hypothesis = word(phones_text)
    model = gradients.SequenceJacobian(shape, target, block)
    model.residual(hypothesis)
    began = time.perf_counter()
    cpu_began = time.process_time()
    exact, residual = model.jacobian(hypothesis, with_residual=True)
    wall = time.perf_counter() - began
    cpu = time.process_time() - cpu_began
    began = time.perf_counter()
    model.residual(hypothesis)
    one = time.perf_counter() - began
    print("Jacobian %d x %d: %.2f s wall, %.2f s CPU (block %d); one native residual %.3f s" % (
        exact.shape[0], exact.shape[1], wall, cpu, block, one))
    steps = shape.steps(len(hypothesis))
    phones = shape.phones
    count = shape.count
    agreement = numpy.zeros(len(hypothesis))
    scales = numpy.zeros(len(hypothesis))
    for column in range(len(hypothesis)):
        size = 1e-6 if column < phones * count else 1e-8
        up = hypothesis.copy()
        down = hypothesis.copy()
        up[column] += size
        down[column] -= size
        central = (model.residual(up) - model.residual(down)) / (2.0 * size)
        scale = max(numpy.linalg.norm(central), numpy.linalg.norm(exact[:, column]), 1e-300)
        agreement[column] = numpy.linalg.norm(central - exact[:, column]) / scale
        scales[column] = numpy.linalg.norm(exact[:, column])
    tubes = agreement[:phones * count]
    times = agreement[phones * count:]
    print("log-area columns (central differences, step 1e-6): relative error median %.2g, 90%% %.2g, max %.2g" % (
        numpy.median(tubes), numpy.percentile(tubes, 90), numpy.max(tubes)))
    print("time columns (central differences, step 1e-8 s): relative errors %s" % " ".join(
        "%.2g" % value for value in times))
    worst = numpy.argsort(agreement)[::-1][:6]
    print("worst columns: %s" % " ".join("%d:%.2g (|J| %.3g)" % (c, agreement[c], scales[c]) for c in worst))
    numpy.savez(os.path.join(gradients.SCRATCH, "agreement_%s.npz" % phones_text.replace(",", "_").lower()),
                agreement=agreement, scales=scales, hypothesis=hypothesis)


def sweep(phones_text, columns=(2, 48, 100, 138, 141)):
    """Central differences of the native residual at several steps against
    the exact columns: truncation error falls as the step squared, rounding
    grows as one over it; the exact column sits at the floor of both."""
    shape, target, truth, hypothesis = word(phones_text)
    model = gradients.SequenceJacobian(shape, target)
    exact = model.jacobian(hypothesis)
    tube_count = shape.phones * shape.count
    for column in columns:
        sizes = (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8) if column < tube_count else (1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
        errors = []
        for size in sizes:
            up = hypothesis.copy()
            down = hypothesis.copy()
            up[column] += size
            down[column] -= size
            central = (model.residual(up) - model.residual(down)) / (2.0 * size)
            errors.append(numpy.linalg.norm(central - exact[:, column]) / max(numpy.linalg.norm(exact[:, column]),
                                                                               1e-300))
        print("column %d (%s): %s" % (column, "log area" if column < tube_count else "time",
                                      "  ".join("h %.0e: %.2g" % (size, error) for size, error in zip(sizes, errors))))


def cpu_now():
    """This process's CPU seconds and its finished children's (sustain)."""
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def timing(phones_text, repeats=5):
    """CPU seconds (this process and sustain's), interleaved repeats, medians:
    the reference residual (numba peel, sustain render), the native
    residual, the exact Jacobian; and the search's forward-difference
    Jacobian's cost as (columns + 1) reference residuals."""
    shape, target, truth, hypothesis = word(phones_text)
    model = gradients.SequenceJacobian(shape, target)
    model.jacobian(hypothesis)
    reference_residual(shape, hypothesis, target)
    costs = {"reference residual": [], "native residual": [], "exact Jacobian": []}
    walls = {"reference residual": [], "native residual": [], "exact Jacobian": []}
    for repeat in range(repeats):
        for name in costs:
            began = cpu_now()
            wall = time.perf_counter()
            if name == "reference residual":
                reference_residual(shape, hypothesis, target)
            elif name == "native residual":
                model.residual(hypothesis)
            else:
                model.jacobian(hypothesis)
            costs[name].append(cpu_now() - began)
            walls[name].append(time.perf_counter() - wall)
    for name in costs:
        print("%s: %.3f s CPU, %.3f s wall (medians of %d)" % (name, numpy.median(costs[name]),
                                                              numpy.median(walls[name]), repeats))
    columns = len(hypothesis)
    reference = numpy.median(costs["reference residual"])
    exact = numpy.median(costs["exact Jacobian"])
    native = numpy.median(costs["native residual"])
    print("forward-difference Jacobian (%d columns + base): reference %.1f s CPU, native %.1f s CPU; exact %.2f s CPU "
          "= %.1f reference residuals = %.1f native residuals; %.0fx less CPU than the reference differences" % (
              columns, (columns + 1) * reference, (columns + 1) * native, exact, exact / reference, exact / native,
              (columns + 1) * reference / exact))


class Recorder:
    """A residual function that keeps each evaluation's level (dB against
    the glottal flow's swing, as sequence_recover.level)."""

    def __init__(self, model, target):
        self.model = model
        self.target = target
        self.levels = []

    def __call__(self, parameters):
        residual = self.model.residual(parameters)
        self.levels.append(sequence_recover.level(residual, self.target.reference))
        return residual


class ForwardDifferences:
    """The search's own Jacobian (tube_recover.Problem.jacobian: forward
    differences at SequenceShape.steps) of the native residual, serially."""

    def __init__(self, model):
        self.model = model

    def __call__(self, parameters):
        steps = self.model.shape.steps(len(parameters))
        base = self.model.residual(parameters)
        columns = []
        for column in range(len(parameters)):
            moved = numpy.array(parameters, float)
            moved[column] += steps[column]
            columns.append((self.model.residual(moved) - base) / steps[column])
        return numpy.stack(columns, axis=1)


def descend(phones_text, budget=30):
    """Levenberg-Marquardt (tube_recover.stage's settings) from the
    hypothesis with the exact Jacobian and with the search's forward
    differences, both on the native residual."""
    from scipy.optimize import least_squares
    shape, target, truth, hypothesis = word(phones_text)
    model = gradients.SequenceJacobian(shape, target)
    tube_count = shape.phones * shape.count
    print("start: %.1f dB; truth %.1f dB" % (sequence_recover.level(model.residual(hypothesis), target.reference),
                                             sequence_recover.level(model.residual(truth), target.reference)))
    for name, jacobian in (("exact", model.jacobian), ("forward differences", ForwardDifferences(model))):
        recorder = Recorder(model, target)
        began = time.perf_counter()
        cpu_began = time.process_time()
        result = least_squares(recorder, hypothesis, jac=jacobian, method="lm", max_nfev=budget, xtol=1e-12,
                               ftol=1e-12, gtol=1e-14)
        wall = time.perf_counter() - began
        cpu = time.process_time() - cpu_began
        found = result.x
        area_error = numpy.sqrt(numpy.mean((numpy.exp(found[:tube_count]) - numpy.exp(truth[:tube_count])) ** 2))
        time_error = numpy.max(numpy.abs(found[tube_count:] - truth[tube_count:]))
        print("%s Jacobian: %d residuals, %d Jacobians, %.0f s wall, %.0f s CPU; levels %s; final %.1f dB, "
              "areas rms %.4f cm^2 off, times at most %.4f ms off" % (
                  name, result.nfev, result.njev, wall, cpu, " ".join("%.1f" % value for value in recorder.levels),
                  recorder.levels[-1], area_error, 1e3 * time_error))
        sys.stdout.flush()


def reference(phones_text, columns=6):
    """A few columns of the reference pipeline's forward differences (sustain,
    the search's own steps) against the exact Jacobian, and their cost."""
    shape, target, truth, hypothesis = word(phones_text)
    model = gradients.SequenceJacobian(shape, target)
    exact = model.jacobian(hypothesis)
    base = reference_residual(shape, hypothesis, target)
    steps = shape.steps(len(hypothesis))
    phones = shape.phones
    count = shape.count
    chosen = [5, count + 20, 2 * count + 40, phones * count, phones * count + 2, phones * count + shape.time_count - 1]
    seconds = []
    for column in chosen[:columns]:
        for size in (steps[column], 0.01 * steps[column]):
            moved = hypothesis.copy()
            moved[column] += size
            began = time.perf_counter()
            forward = (reference_residual(shape, moved, target) - base) / size
            seconds.append(time.perf_counter() - began)
            error = numpy.linalg.norm(forward - exact[:, column]) / max(numpy.linalg.norm(exact[:, column]), 1e-300)
            print("column %d, reference forward difference step %.3g: relative error %.3g against the exact" % (
                column, size, error))
    print("reference residual: %.3f s each (median of %d)" % (numpy.median(seconds), len(seconds)))


if __name__ == "__main__":
    command = sys.argv[1]
    phones_text = sys.argv[2] if len(sys.argv) > 2 else "AA,IY,UW"
    if command == "check":
        check(phones_text)
    elif command == "jacobian":
        jacobian(phones_text, int(sys.argv[3]) if len(sys.argv) > 3 else gradients.BLOCK)
    elif command == "timing":
        timing(phones_text)
    elif command == "descend":
        descend(phones_text)
    elif command == "sweep":
        if len(sys.argv) > 3:
            sweep(phones_text, tuple(int(value) for value in sys.argv[3].split(",")))
        else:
            sweep(phones_text)
    elif command == "reference":
        reference(phones_text)
