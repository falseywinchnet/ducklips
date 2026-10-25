"""Checks and costs for the driver's shortcuts (NOTES_driver.md).

usage: driver_bench.py WORD TEST[,TEST...]   (WORD as sequence_recover's;
POOL workers for the pool tests)

Tests:
  cost      one residual of each kind the word's search uses: wall, CPU,
            and (cProfile) where it goes
  truncate  each residual kind with renders stopped at the compared
            stretch's end against the whole sound's: bitwise
  stage     a posture stage through least_squares (the reference) and
            through leastsq with the memo: parameters bitwise, evaluations
  inert     the joint's Jacobian with and without the inert columns (bend 0
            and bend 1): bitwise
  locality  a column per phone and per time: the first residual row that
            differs from the base against column_locality's resume step
  first     the first vowel's equation error on the sound cut at its
            window's end plus the cache's reach against the whole: bitwise
"""
import cProfile
import io
import os
import pstats
import sys
import time
from multiprocessing import Pool
from typing import Any, Callable, Dict, List, Tuple

import numpy
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import column_locality
import driver_ledger
import posture_recover
import sequence_recover
import tube_recover

STEP_RATE = sequence_recover.STEP_RATE


def hinted_parameters(word: Any) -> numpy.ndarray:
    """A hypothesis near the truth (the hinted times, the true tubes): the
    point the searches work around."""
    return word.shape.join(word.true_tubes, word.hinted_arrivals, word.hinted_leavings,
                           numpy.where(word.true_velum > 0.0, word.true_velum + 0.02, 0.0))


def targets(word: Any) -> Dict[str, Any]:
    """One target of each kind the search uses on this word."""
    flow_out = word.flow_out
    out: Dict[str, Any] = {"whole": word.whole}
    begin, span = sequence_recover.window(word.hinted_leavings[0] - 0.02, word.hinted_arrivals[1] + 0.03)
    out["grid exact"] = tube_recover.Target(flow_out, begin, span, 0.0, 1.0, True)
    out["grid low"] = tube_recover.Target(flow_out, begin, span, 1500.0, 1.0, True)
    start, size = sequence_recover.window(0.02, word.hinted_leavings[min(2, word.phones - 1)] - sequence_recover.CLEAR)
    out["prefix"] = tube_recover.Target(flow_out, start, size, 0.0, 1.0, True, 1.0)
    sealed = [k for k in sorted(word.stops) if k + 1 < word.phones]
    if sealed:
        index = sealed[0] + 1
        first = max(int((word.hinted_arrivals[index - 1] - 0.005) * STEP_RATE), 0)
        inside = tube_recover.Target(flow_out, first, int((word.hinted_leavings[index - 1] + 0.015) * STEP_RATE) - first,
                                     0.0, 1.0, False)
        inside.output_from = first
        out["release scan"] = inside
    return out


def timed(function: Callable[[], Any]) -> Tuple[Any, float, float]:
    wall = time.time()
    cpu = driver_ledger.process_cpu()
    value = function()
    return value, time.time() - wall, driver_ledger.process_cpu() - cpu


class ResidualCall:
    """A named zero-argument call of a word's residual (no lambdas)."""

    def __init__(self, word: Any, parameters: numpy.ndarray, target: Any) -> None:
        self.word = word
        self.parameters = parameters
        self.target = target

    def __call__(self) -> numpy.ndarray:
        if self.target.output_from is not None:
            tubes, arrivals, leavings = self.word.shape.split(self.parameters)
            return sequence_recover.word_residual(tubes, arrivals, leavings, self.word.shape.velum(self.parameters),
                                                  self.word.length, self.target)
        return self.word.shape.residual(self.parameters, self.target)


def test_cost(word: Any) -> None:
    parameters = hinted_parameters(word)
    for name, target in targets(word).items():
        call = ResidualCall(word, parameters, target)
        call()   # warm (numba, the nose's constants)
        profile = cProfile.Profile()
        profile.enable()
        _, wall, cpu = timed(call)
        profile.disable()
        text = io.StringIO()
        pstats.Stats(profile, stream=text).sort_stats("tottime").print_stats(6)
        lines = [line for line in text.getvalue().split("\n") if line.strip()][-7:]
        print("cost %-13s wall %.3f s, CPU %.3f s (engine included)" % (name, wall, cpu))
        for line in lines:
            print("    " + line[:150])
        sys.stdout.flush()


def test_truncate_offsets(word: Any) -> None:
    """truncate's check at the first transition's grid points (leaving and
    arrival moved), on the whole and the grid's targets."""
    base = hinted_parameters(word)
    chosen = targets(word)
    leaving = word.phones * word.count
    for leaving_offset in (-0.0016, -0.0002, 0.0, 0.0004, 0.003):
        for arrival_offset in (-0.0008, 0.0, 0.0012):
            parameters = base.copy()
            parameters[leaving] += leaving_offset
            parameters[leaving + 1] += arrival_offset
            for name in ("grid exact", "whole"):
                call = ResidualCall(word, parameters, chosen[name])
                sequence_recover.TRUNCATE = False
                full = call()
                sequence_recover.TRUNCATE = True
                cut = call()
                print("truncate at %+.4f %+.4f %-11s bitwise %s" % (leaving_offset, arrival_offset, name,
                                                                   numpy.array_equal(full, cut)))
                sys.stdout.flush()


def test_truncate(word: Any) -> None:
    parameters = hinted_parameters(word)
    for name, target in targets(word).items():
        call = ResidualCall(word, parameters, target)
        sequence_recover.TRUNCATE = False
        full, full_wall, _ = timed(call)
        sequence_recover.TRUNCATE = True
        cut, cut_wall, _ = timed(call)
        print("truncate %-13s bitwise %s (max |difference| %.3g); wall %.3f -> %.3f s" % (
            name, numpy.array_equal(full, cut), float(numpy.max(numpy.abs(full - cut))) if len(full) == len(cut) else -1.0,
            full_wall, cut_wall))
        sys.stdout.flush()


def test_renderer(word: Any) -> None:
    """Each residual kind through the subprocess and through the in-process
    engine (posture_recover.RENDERER): bitwise, and the time (medians of 3)."""
    parameters = hinted_parameters(word)
    for name, target in targets(word).items():
        call = ResidualCall(word, parameters, target)
        values = {}
        for renderer in ("sustain", "engine", "incremental") * 3:
            posture_recover.RENDERER = renderer
            value, wall, cpu = timed(call)
            values.setdefault(renderer, []).append((value, wall, cpu))
        posture_recover.RENDERER = "sustain"
        same = all(numpy.array_equal(values["sustain"][0][0], entry[0])
                   for entry in values["engine"] + values["incremental"])
        print("renderer %-13s bitwise %s; CPU s (medians of 3): sustain %.3f, engine %.3f, incremental %.3f "
              "(the same hypothesis again: its best case)" % (
                  name, same, float(numpy.median([e[2] for e in values["sustain"]])),
                  float(numpy.median([e[2] for e in values["engine"]])),
                  float(numpy.median([e[2] for e in values["incremental"]]))))
        sys.stdout.flush()


class CountingResidual:
    """A Problem's residual counter for least_squares (the reference path)."""

    def __init__(self, problem: Any) -> None:
        self.problem = problem

    def __call__(self, parameters: numpy.ndarray) -> numpy.ndarray:
        return self.problem.residual(parameters)


def test_stage(word: Any, pool: Any) -> None:
    hint = word.hints[word.names[0]]
    begin, span = sequence_recover.window(word.hinted_arrivals[0] + sequence_recover.SETTLE,
                                          word.hinted_leavings[0] - sequence_recover.CLEAR)
    guide = tube_recover.SpectrumTarget(word.flow_out, begin, span)
    shape = tube_recover.PostureShape(hint, tube_recover.WIDEN, word.count)
    start = numpy.full(len(posture_recover.FIELDS), 0.1)
    results: Dict[str, Tuple[numpy.ndarray, int, float]] = {}
    for exact in (False, True):
        tube_recover.EXACT_SHORTCUTS = exact
        wall = time.time()
        found, used = tube_recover.stage(pool, shape, start, guide, 15)
        results["leastsq+memo" if exact else "least_squares"] = (found, used, time.time() - wall)
    tube_recover.EXACT_SHORTCUTS = True
    reference = results["least_squares"]
    fast = results["leastsq+memo"]
    print("stage (posture, 9 parameters, 15 nfev): least_squares %d evaluations %.1f s; leastsq+memo %d evaluations "
          "%.1f s; parameters bitwise %s" % (reference[1], reference[2], fast[1], fast[2],
                                             numpy.array_equal(reference[0], fast[0])))
    sys.stdout.flush()


def test_inert(word: Any, pool: Any) -> None:
    parameters = hinted_parameters(word)
    for bend in (0.0, 1.0):
        target = tube_recover.Target(word.flow_out, word.start, word.size, 0.0, 1.0, True, bend)
        jacobians = {}
        for exact in (False, True):
            tube_recover.EXACT_SHORTCUTS = exact
            problem = tube_recover.Problem(word.shape, pool, target)
            wall = time.time()
            jacobians[exact] = (problem.jacobian(parameters), problem.evaluations, time.time() - wall)
        tube_recover.EXACT_SHORTCUTS = True
        print("inert columns, bend %.1f: %d -> %d evaluations, %.1f -> %.1f s; Jacobian bitwise %s" % (
            bend, jacobians[False][1], jacobians[True][1], jacobians[False][2], jacobians[True][2],
            numpy.array_equal(jacobians[False][0], jacobians[True][0])))
        sys.stdout.flush()


def test_resume(word: Any, pool: Any) -> None:
    """The joint's Jacobian with its columns resumed from the base's record
    (batches per worker) against every column whole: bitwise, and the time."""
    parameters = hinted_parameters(word)
    for bend in (0.0, 1.0):
        target = tube_recover.Target(word.flow_out, word.start, word.size, 0.0, 1.0, True, bend)
        jacobians = {}
        for resume in (False, True):
            sequence_recover.RESUME = resume
            problem = tube_recover.Problem(word.shape, pool, target)
            wall = time.time()
            cpu = driver_ledger.LEDGER.records.get("(top)")
            before = 0.0 if cpu is None else cpu.worker_cpu
            jacobians[resume] = (problem.jacobian(parameters), problem.evaluations, time.time() - wall,
                                 driver_ledger.LEDGER.records["(top)"].worker_cpu - before)
        sequence_recover.RESUME = True
        print("resume, bend %.1f: whole %d evaluations %.1f s wall %.1f s CPU; resumed %d evaluations %.1f s wall "
              "%.1f s CPU; Jacobian bitwise %s" % (
                  bend, jacobians[False][1], jacobians[False][2], jacobians[False][3], jacobians[True][1],
                  jacobians[True][2], jacobians[True][3], numpy.array_equal(jacobians[False][0], jacobians[True][0])))
        sys.stdout.flush()


def locality_task(task: Tuple[Any, numpy.ndarray, Any]) -> numpy.ndarray:
    shape, parameters, target = task
    return shape.residual(parameters, target)


def test_locality(word: Any, pool: Any) -> None:
    parameters = hinted_parameters(word)
    target = word.whole
    shape = word.shape
    steps = shape.steps(len(parameters))
    columns: List[int] = []
    for phone in range(word.phones):
        if phone not in word.stops:
            columns.extend([phone * word.count + 5, phone * word.count + word.count - 3])
    columns.extend(range(word.phones * word.count, word.phones * word.count + shape.time_count))
    columns.extend(range(word.phones * word.count + shape.time_count, len(parameters)))
    moved_points = []
    for column in columns:
        moved = parameters.copy()
        moved[column] += steps[column]
        moved_points.append(moved)
    locality = column_locality.ColumnLocality(shape, target)
    predicted = locality.resume_steps(parameters, moved_points)
    values = pool.map(locality_task, [(shape, point, target) for point in [parameters] + moved_points], chunksize=1)
    base = values[0]
    scale = float(numpy.max(numpy.abs(base)))
    for column, step, value in zip(columns, predicted, values[1:]):
        # fixed_clip removes the stretch's mean: a change anywhere shifts
        # every row by the mean's change. The first row whose difference
        # from the base departs from that constant shift:
        shift = value - base
        departs = numpy.nonzero(numpy.abs(shift - shift[0]) > 1e-9 * scale)[0]
        first = int(departs[0]) + target.start if len(departs) else -1
        print("locality column %3d: resume step %6d (%.4f s); first departing row step %6d (%.4f s); safe %s" % (
            column, step, step / STEP_RATE, first, first / STEP_RATE, step <= first or first < 0))
    print("(rows are mean-removed over the stretch (fixed_clip): a resumed column must carry the base's raw "
          "difference before its resume step and remove the mean after; see NOTES_driver.md)")
    sys.stdout.flush()


def test_first(word: Any) -> None:
    hint = word.hints[word.names[0]]
    begin, span = sequence_recover.window(word.hinted_arrivals[0] + sequence_recover.SETTLE,
                                          word.hinted_leavings[0] - sequence_recover.CLEAR)
    reach = begin + span + sequence_recover.FIRST_VOWEL_REACH * word.count
    for label, tube in (("truth", (word.true_tubes[0], word.length)),
                        ("hint", posture_recover.nodes(hint, word.count))):
        full = tube_recover.Target(word.flow_out, begin, span, 0.0, 1.0, True)
        cut = tube_recover.Target(numpy.array(word.flow_out[:reach]), begin, span, 0.0, 1.0, True)
        a = tube_recover.residual_for(tube, full)
        b = tube_recover.residual_for(tube, cut)
        spectrum_full = tube_recover.residual_for(tube, tube_recover.SpectrumTarget(word.flow_out, begin, span))
        spectrum_cut = tube_recover.residual_for(tube, tube_recover.SpectrumTarget(numpy.array(word.flow_out[:reach]),
                                                                                    begin, span))
        print("first vowel (%s): equation error bitwise on the sound cut at +%d steps: %s; spectrum: %s" % (
            label, reach - begin - span, numpy.array_equal(a, b), numpy.array_equal(spectrum_full, spectrum_cut)))
        for margin in (0, 10, 25, 50, 100):
            short = tube_recover.Target(numpy.array(word.flow_out[:begin + span + margin]), begin, span, 0.0, 1.0, True)
            c = tube_recover.residual_for(tube, short)
            print("    cut at +%3d steps: bitwise %s" % (margin, numpy.array_equal(a, c)))
    sys.stdout.flush()


def main() -> None:
    word = sequence_recover.build_word(sys.argv[1].split(","))
    sequence_recover.configure_targets(word)
    tests = sys.argv[2].split(",")
    print("%s: %d phones, %.2f s, %d steps" % (" ".join(word.specs), word.phones, word.seconds, len(word.flow_out)))
    if "cost" in tests:
        test_cost(word)
    if "truncate" in tests:
        test_truncate(word)
    if "offsets" in tests:
        test_truncate_offsets(word)
    if "first" in tests:
        test_first(word)
    if "renderer" in tests:
        test_renderer(word)
    if any(test in tests for test in ("stage", "inert", "locality", "resume")):
        with Pool(driver_ledger.POOL) as pool:
            if "stage" in tests:
                test_stage(word, pool)
            if "inert" in tests:
                test_inert(word, pool)
            if "locality" in tests:
                test_locality(word, pool)
            if "resume" in tests:
                test_resume(word, pool)


if __name__ == "__main__":
    main()
