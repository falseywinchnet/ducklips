"""How much the sound backs each part of the recovered timing: for each
clock, a local nudge of its time field (a raised-cosine bump one quantum
high over BUMP_TICKS ticks) at each place, and the equation error's change
it makes against the sound's scale (dB). A part whose nudge changes almost
nothing is a choice the sound does not make (the least-action step makes
it there), not a measurement.

Usage: observability.py SPECS STAGE (the saved solution of fields_word.py)
"""
import multiprocessing
import os
import sys
from typing import List, Optional, Tuple

os.environ.setdefault("RENDERER", "engine")
import numpy

import fields
import fields_word
import pyramid_word

BUMP_TICKS = 10   # ticks a nudge spans (10 ms: about a fifth of a transition)
STRIDE = 5        # ticks between nudges

_problem: Optional[fields_word.TimeProblem] = None
_vector: Optional[numpy.ndarray] = None
_base: Optional[numpy.ndarray] = None


def bump(ticks: int, centre: int) -> numpy.ndarray:
    """A raised-cosine bump of height one over BUMP_TICKS ticks at centre."""
    shape = numpy.zeros(ticks)
    for tick in range(max(centre - BUMP_TICKS // 2, 0), min(centre + BUMP_TICKS // 2, ticks)):
        phase = (tick - centre) / (BUMP_TICKS / 2.0)
        shape[tick] = 0.5 * (1.0 + numpy.cos(numpy.pi * phase))
    return shape


def nudged_residual(problem: fields_word.TimeProblem, vector: numpy.ndarray, clock: int,
                    shape: numpy.ndarray) -> numpy.ndarray:
    """The residual with one clock's time field nudged by `shape` quanta."""
    clocks = problem.times(vector)
    quantum = fields_word.time_quantum(problem.context) * pyramid_word.TICK
    clocks[clock] = clocks[clock] + quantum * shape
    context = problem.context
    rows = problem.rows_of(vector) if hasattr(problem, "rows_of") else problem.rows
    space = problem.space_of(vector) if hasattr(problem, "space_of") else problem.space
    side = clocks[1] if len(clocks) > 1 else None
    command = fields.read_postures(rows, context.count, space, fields.same_time(clocks[0], context.count), side)
    difference = pyramid_word.equation_error(context, command)[context.span]
    return difference / context.scale


def worker_start(specs: List[str], stage: str, knots: list, used: int, vector: numpy.ndarray) -> None:
    global _problem, _vector, _base
    fields_word.worker_start(specs, knots, stage, used)
    _problem = fields_word._problem
    _vector = vector
    _base = _problem.residual(vector)


def worker_nudge(task: Tuple[int, int]) -> float:
    clock, centre = task
    changed = nudged_residual(_problem, _vector, clock, bump(_problem.context.ticks, centre))
    return float(numpy.linalg.norm(changed - _base))


def main() -> None:
    specs = sys.argv[1].split(",")
    stage = sys.argv[2]
    saved = sys.argv[4] if len(sys.argv) > 4 else os.path.join(
        pyramid_word.sequence_recover.SCRATCH, "fields_%s_%s.npy" % (stage, "_".join(specs).replace("!", "")))
    # Knots as fields_word stores them: (tick, value) for the time stage,
    # (tick, node, phone below) for the body stages.
    knots = []
    if len(sys.argv) > 3:
        for knot in sys.argv[3].split(";"):
            parts = knot.split(":")
            knots.append(tuple(float(x) for x in parts) if stage == "time" else
                         (float(parts[0]), int(parts[1]), int(parts[2])))
    vector = numpy.load(saved)
    context = pyramid_word.Context(specs, ["areas"])
    fields_word.worker_start(specs, knots, stage, context.levels)
    problem = fields_word._problem
    clocks = problem.clocks()
    tasks = [(clock, centre) for clock in range(clocks) for centre in range(0, context.ticks, STRIDE)]
    with multiprocessing.Pool(fields_word.POOL, initializer=worker_start,
                              initargs=(specs, stage, knots, context.levels, vector)) as pool:
        changes = pool.map(worker_nudge, tasks)
    names = ["mouth", "folds"]
    for clock in range(clocks):
        print("%s clock: the sound's change for a one-quantum nudge (dB against the scale), every %d ms:" % (
            names[clock], STRIDE))
        values = [changes[index] for index, (c, _) in enumerate(tasks) if c == clock]
        line = []
        for position, value in enumerate(values):
            line.append("%d:%.0f" % (position * STRIDE, 20.0 * numpy.log10(max(value, 1e-300))))
        print("  " + " ".join(line), flush=True)


if __name__ == "__main__":
    main()
