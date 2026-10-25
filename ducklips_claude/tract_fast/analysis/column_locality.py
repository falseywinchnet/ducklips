"""Locality of a word's Jacobian columns: a parameter of phone k acts only
from its earliest keyframe on, so a finite-difference column's residual is
the base point's up to there and could be computed from a saved state
(a resume point) instead of from the sound's start.

This module holds the driver's side of that contract:

* first_divergence(base_inputs, moved_inputs): the first engine step whose
  inputs differ between the base point and a moved one. The inputs are the
  keyframe tracks the engine and the peels read (the tube's AREASEQ, the
  velum's VELUMSEQ, the folds' FOLDKEYS, the jet's JETSEQ), all linear
  between keyframes: a change at keyframe i alters the track from keyframe
  i - 1's time on.
* resume_step(divergence, nodes): where a column's computation may start.
  The peels read the sound ahead of the step they produce (a node count of
  steps per pass; output_residual already peels 4 node counts past its
  start), so the residual is the base's only up to the divergence less
  that reach. LOOKAHEAD_NODES is the driver's conservative bound; OPEN:
  the kernels' owners to state the exact reach.
* ResumableResidual: what a residual must offer for resumed columns, and
  FullRecompute, the stub that recomputes everything (today's behaviour),
  so the driver's Jacobian can be written against the interface now and the
  kernels' snapshots plugged in later.
* ColumnLocality: per Jacobian, the projected share of steps the columns
  could skip (MEASURED in the driver's runs, reported with the ledger).
"""
import math
from typing import Any, Dict, List, Optional, Protocol, Sequence, Tuple

import numpy

import driver_ledger
import reverse_core

STEP_RATE = reverse_core.STEP_RATE
LOOKAHEAD_NODES = 4   # node counts of steps the peels read ahead (output_residual's reach)

Track = Tuple[numpy.ndarray, numpy.ndarray]   # keyframe times (s), values (one row a keyframe)


def track_divergence(base: Optional[Track], moved: Optional[Track]) -> float:
    """The first time (s) at which two linear keyframe tracks can differ
    (math.inf: identical)."""
    if base is None and moved is None:
        return math.inf
    if base is None or moved is None:
        return 0.0
    base_times = numpy.asarray(base[0], float)
    moved_times = numpy.asarray(moved[0], float)
    base_values = numpy.asarray(base[1], float)
    moved_values = numpy.asarray(moved[1], float)
    if len(base_times) != len(moved_times) or base_values.shape != moved_values.shape:
        return 0.0
    for index in range(len(base_times)):
        same_time = base_times[index] == moved_times[index]
        same_value = numpy.array_equal(base_values[index], moved_values[index])
        if not (same_time and same_value):
            if index == 0:
                return 0.0
            held = (numpy.array_equal(base_values[index - 1], base_values[index])
                    and numpy.array_equal(moved_values[index - 1], moved_values[index]))
            if same_value and held:
                # A held stretch's end moved: the track is flat up to the
                # earlier of the two times.
                return float(min(base_times[index], moved_times[index]))
            # Linear from keyframe index - 1: the track bends from there.
            return float(min(base_times[index - 1], moved_times[index - 1]))
    return math.inf


def first_divergence(base_inputs: Dict[str, Optional[Track]], moved_inputs: Dict[str, Optional[Track]]) -> float:
    """The first time (s) any of the engine's input tracks differs."""
    earliest = math.inf
    for name in set(base_inputs) | set(moved_inputs):
        earliest = min(earliest, track_divergence(base_inputs.get(name), moved_inputs.get(name)))
    return earliest


def resume_step(divergence: float, nodes: int) -> int:
    """The step from which a column's residual must be recomputed."""
    if math.isinf(divergence):
        return -1   # no input differs: the column is the base's
    return max(int(math.floor(divergence * STEP_RATE)) - LOOKAHEAD_NODES * nodes, 0)


class ResumableResidual(Protocol):
    """A residual that can start from a saved state.

    run_with_snapshots(parameters, target, steps): the residual at
    `parameters` and, for each step in `steps`, an opaque picklable
    snapshot of everything the residual's computation carries across that
    step (the engine's fold and tube state, the peels' states, the forced
    waves so far, the observer's recorded state, the march's arrays).
    resume(parameters, target, snapshot): the residual at `parameters`
    (which differ from the snapshot's base only in inputs from
    snapshot.step on), computed from the snapshot; bitwise equal to a full
    computation."""

    def run_with_snapshots(self, parameters: numpy.ndarray, target: Any,
                           steps: Sequence[int]) -> Tuple[numpy.ndarray, Dict[int, Any]]:
        ...

    def resume(self, parameters: numpy.ndarray, target: Any, snapshot: Any) -> numpy.ndarray:
        ...


class FullRecompute:
    """The stub: no snapshots, every column from the sound's start."""

    def __init__(self, shape: Any) -> None:
        self.shape = shape

    def run_with_snapshots(self, parameters: numpy.ndarray, target: Any,
                           steps: Sequence[int]) -> Tuple[numpy.ndarray, Dict[int, Any]]:
        snapshots: Dict[int, Any] = {}
        for step in steps:
            snapshots[step] = None
        return self.shape.residual(parameters, target), snapshots

    def resume(self, parameters: numpy.ndarray, target: Any, snapshot: Any) -> numpy.ndarray:
        return self.shape.residual(parameters, target)


# Projected savings, by the ledger's top stage: [columns, steps a column
# computes from the start, steps it would compute from its resume point].
STATISTICS: Dict[str, List[float]] = {}


class ColumnLocality:
    """A stage's columns' resume steps (only for a shape that names its
    input tracks: SequenceShape.input_tracks(parameters, target))."""

    def __init__(self, shape: Any, target: Any) -> None:
        self.shape = shape
        self.target = target
        self.enabled = hasattr(shape, "input_tracks") and hasattr(target, "flow_out")

    def resume_steps(self, point: numpy.ndarray, moved_points: Sequence[numpy.ndarray]) -> List[int]:
        base = self.shape.input_tracks(point, self.target)
        nodes = self.shape.count
        steps: List[int] = []
        for moved in moved_points:
            steps.append(resume_step(first_divergence(base, self.shape.input_tracks(moved, self.target)), nodes))
        return steps

    def observe(self, point: numpy.ndarray, moved_points: Sequence[numpy.ndarray]) -> None:
        if not self.enabled:
            return
        horizon = self.target.start + self.target.length
        steps = self.resume_steps(point, moved_points)
        key = driver_ledger.LEDGER.stack[0] if driver_ledger.LEDGER.stack else "(top)"
        entry = STATISTICS.setdefault(key, [0.0, 0.0, 0.0])
        for step in steps:
            entry[0] += 1
            entry[1] += horizon
            entry[2] += 0 if step < 0 else max(horizon - step, 0)


def report(stream: Any) -> None:
    stream.write("column locality (projected; resume at the first input divergence less %d node counts):\n"
                 % LOOKAHEAD_NODES)
    for key in sorted(STATISTICS):
        columns, full, resumed = STATISTICS[key]
        stream.write("  %-30s %7d columns, steps from start %.3g, from resume %.3g: %.1f%% of the columns' steps\n" % (
            key, columns, full, resumed, 100.0 * resumed / max(full, 1.0)))
    stream.flush()
