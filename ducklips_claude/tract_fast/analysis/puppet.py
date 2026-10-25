"""The puppet: the inverse in the family of articulator motions.

Phase B of the 2026-10-24 program turn (NOTES_pyramid.md). A word is a
sequence of targets with event times (each phone's leaving and the next's
arrival); every articulator moves between targets with a half-cosine, at
its own lag and with its own extra move time (the per-node traverse: the
owner's "friction"; both in seconds, so one secant serves every unknown).
The word's unknowns are its event times, two a transition; the voice's are
the lags and extras, shared across words and held here (Phase D learns
them, voice.py). The sound is matched by the equation error, as in fields_word.py,
and the secant for the Jacobian descends by decades from the hint's
uncertainty (HINT_OFFSET, 7.3 ms) to one engine step.

Why this family and not a free time field (fields_word.TimeProblem, one
knot a tick): the zero-flow residual is a linear meter of front mismatch
(a clock error d at one tick reads about 20 log10 d dB: 1e-6 -> -109 dB,
3.9e-4 -> -57 dB; outside_in.py, step 43), so a free field has hundreds of
directions the sound barely hears and goes where the onset's roughness
sends it (Alfa: no step kept at any secant up to a tick, step 42). An event
time moves a whole transition, averaging the roughness over the move.

HYPOTHESIS under test: Alfa (seal-less, blocked at -19.5 dB) descends in
this family; Echo, Golf and Papa do not regress.
"""

import multiprocessing
import os
import sys
from typing import List, Optional, Sequence, Tuple

import numpy

import fields
import fields_word
import pyramid_word
import sequence_recover

SECANTS = (0.01, 0.001, 0.0001, 1.0 / sequence_recover.STEP_RATE)   # s; decades from the hint's error to a step
GROUPS = (("lips", 0, 2), ("blade", 2, 8), ("front", 8, 16), ("back", 16, 26), ("pharynx", 26, 10000))
# articulator groups by distance from the lips in nodes (0.37 cm a node): the
# lips' two nodes, the tongue blade's six behind them, the front and back of
# the body, the pharynx. Phase D's shared unit: one word cannot determine 46
# node lags (precisions 0.5-86 ms in the single-word run, step 51).
DURATION_START = float(numpy.mean(sequence_recover.TRANSITIONS))   # s; the plan's own mean move, every group's start


def node_progress(arrivals: numpy.ndarray, leavings: numpy.ndarray, ticks: int, tick_seconds: float) -> numpy.ndarray:
    """fields.progress for one articulator (its own arrivals and leavings),
    over the ticks at once: phone k held from arrival to leaving, each
    transition a half-cosine in the share of the transition elapsed. Where a
    stretched or lagged move overruns the next phone's hold the later
    transition wins, as in fields.progress."""
    now = numpy.arange(ticks, dtype=numpy.float64) * tick_seconds
    field = numpy.zeros(ticks)
    for k in range(len(arrivals)):
        field[now >= arrivals[k]] = float(k)
    for k in range(len(arrivals) - 1):
        moving = (now > leavings[k]) & (now < arrivals[k + 1])
        share = (now[moving] - leavings[k]) / (arrivals[k + 1] - leavings[k])
        field[moving] = k + 0.5 * (1.0 - numpy.cos(numpy.pi * share))
    return field


class PuppetProblem:
    """fields_word.TimeProblem's interface (size, steps, kinds, times, time,
    command, residual) over event times: vector = [L_0, A_1, L_1, A_2, ...,
    L_{P-2}, A_{P-1}] in seconds, the first arrival and last leaving the
    word's own. Column c (a node, the port, the folds) moves with lag
    lags[c] (s) and its move lasts extras[c] (s) longer than the shared one
    (shorter when negative; never shorter than nothing)."""

    def __init__(self, context: pyramid_word.Context, lags: Optional[numpy.ndarray] = None,
                 extras: Optional[numpy.ndarray] = None, secant: float = SECANTS[0],
                 free_lags: Optional[Sequence[int]] = None, free_extras: Optional[Sequence[int]] = None,
                 crossings: bool = False) -> None:
        self.context = context
        word = context.word
        self.rows = pyramid_word.targets(context, word.true_tubes, word.true_velum)
        self.space = fields.identity_space(context.ticks, context.count)
        self.phones = word.phones
        self.columns = self.rows.shape[1]
        self.lags = numpy.zeros(self.columns) if lags is None else numpy.asarray(lags, dtype=numpy.float64)
        self.extras = numpy.zeros(self.columns) if extras is None else numpy.asarray(extras, dtype=numpy.float64)
        # The columns whose lags, then extras, are unknowns too (after the event times).
        self.free_lags = numpy.array([] if free_lags is None else list(free_lags), dtype=numpy.int64)
        self.free_extras = numpy.array([] if free_extras is None else list(free_extras), dtype=numpy.int64)
        # With `crossings`, a transition that seals or releases a stop is
        # carried as (seal crossing of the stop's node, move duration)
        # instead of (leaving, arrival): the crossing is what the sound
        # hears sharply (step 44's largest singular directions), and the
        # (leaving, arrival) pair let the 10 ms rung slide Delta 15 ms along
        # the crossing-preserving direction (step 46). The crossing's share
        # of the move, where the node's half-cosine meets the seal, is fixed
        # by the targets: f* = acos(1 - 2 (seal - a_k) / (a_{k+1} - a_k)) / pi.
        self.crossings = crossings
        self.sealing_node = numpy.full(self.phones - 1, -1, dtype=numpy.int64)
        self.crossing_share = numpy.zeros(self.phones - 1)
        if crossings:
            for k in range(self.phones - 1):
                node = -1
                if k + 1 in word.stops and k + 1 not in word.narrowings:
                    node = word.stops[k + 1]
                elif k in word.stops and k not in word.narrowings:
                    node = word.stops[k]
                if node < 0:
                    continue
                before, after = self.rows[k, node], self.rows[k + 1, node]
                if after == before:
                    continue
                height = (word.seal - before) / (after - before)
                if not 0.0 < height < 1.0:
                    continue
                self.sealing_node[k] = node
                self.crossing_share[k] = numpy.arccos(1.0 - 2.0 * height) / numpy.pi
        self.first_arrival = float(word.true_arrivals[0])
        self.last_leaving = float(word.true_leavings[-1])
        self.secant = secant
        self.used = 1
        self.knots: List[Tuple[float, float]] = []
        self.window: Optional[numpy.ndarray] = None
        # An adjoint.ResidualCache when evaluations are to be resumed from
        # the last at the first changed tick (bit-identical to the plain
        # residual; 1.6x less CPU, NOTES_adjoint.md §9). One per process:
        # attached by worker_start and main, never pickled.
        self.cache = None
        self.steps = self.difference_steps()

    def events(self) -> int:
        return 2 * (self.phones - 1)

    def size(self) -> int:
        return self.events() + len(self.free_lags) + len(self.free_extras)

    def kinds(self) -> numpy.ndarray:
        """Every unknown a time (an event, a lag, an extra; s): scaled by its secant (kind 0)."""
        return numpy.zeros(self.size())

    def difference_steps(self) -> numpy.ndarray:
        return numpy.full(self.size(), self.secant)

    def events_of(self, vector: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]:
        arrivals = numpy.empty(self.phones)
        leavings = numpy.empty(self.phones)
        arrivals[0] = self.first_arrival
        leavings[-1] = self.last_leaving
        events = self.events()
        leavings[:-1] = vector[0:events:2]
        arrivals[1:] = vector[1:events:2]
        if self.crossings:
            # (crossing, duration) -> (leaving, arrival): the stop's node
            # crosses the seal at its own lag plus f* of its own move.
            lags = self.lags_of(vector)
            extras = self.extras_of(vector)
            for k in numpy.flatnonzero(self.sealing_node >= 0):
                node = self.sealing_node[k]
                crossing, duration = vector[2 * k], vector[2 * k + 1]
                leavings[k] = crossing - lags[node] - self.crossing_share[k] * (duration + extras[node])
                arrivals[k + 1] = leavings[k] + duration
        return arrivals, leavings

    def lags_of(self, vector: numpy.ndarray) -> numpy.ndarray:
        lags = self.lags.copy()
        lags[self.free_lags] = vector[self.events():self.events() + len(self.free_lags)]
        return lags

    def extras_of(self, vector: numpy.ndarray) -> numpy.ndarray:
        extras = self.extras.copy()
        extras[self.free_extras] = vector[self.events() + len(self.free_lags):]
        return extras

    def vector_of(self, arrivals: numpy.ndarray, leavings: numpy.ndarray,
                  lags: Optional[numpy.ndarray] = None, extras: Optional[numpy.ndarray] = None) -> numpy.ndarray:
        vector = numpy.empty(self.size())
        events = self.events()
        vector[0:events:2] = leavings[:-1]
        vector[1:events:2] = arrivals[1:]
        split = events + len(self.free_lags)
        vector[events:split] = (self.lags if lags is None else lags)[self.free_lags]
        vector[split:] = (self.extras if extras is None else extras)[self.free_extras]
        if self.crossings:
            all_lags = self.lags if lags is None else lags
            all_extras = self.extras if extras is None else extras
            for k in numpy.flatnonzero(self.sealing_node >= 0):
                node = self.sealing_node[k]
                duration = arrivals[k + 1] - leavings[k]
                vector[2 * k] = leavings[k] + all_lags[node] + self.crossing_share[k] * (duration + all_extras[node])
                vector[2 * k + 1] = duration
        return vector

    def valid(self, arrivals: numpy.ndarray, leavings: numpy.ndarray) -> bool:
        """Holds no shorter than nothing, moves longer than nothing."""
        return bool(numpy.all(leavings >= arrivals) and numpy.all(arrivals[1:] > leavings[:-1]))

    def valid_vector(self, vector: numpy.ndarray) -> bool:
        arrivals, leavings = self.events_of(vector)
        return self.valid(arrivals, leavings)

    def column_events(self, arrivals: numpy.ndarray, leavings: numpy.ndarray, lags: numpy.ndarray,
                      extras: numpy.ndarray, column: int) -> Tuple[numpy.ndarray, numpy.ndarray]:
        """Column c's own event times: each move starts lags[c] after the
        shared leaving and lasts extras[c] longer than the shared move (at
        least one engine step)."""
        own_leavings = leavings + lags[column]
        own_arrivals = numpy.empty(self.phones)
        own_arrivals[0] = arrivals[0] + lags[column]
        shortest = 1.0 / sequence_recover.STEP_RATE
        own_arrivals[1:] = own_leavings[:-1] + numpy.maximum(arrivals[1:] - leavings[:-1] + extras[column], shortest)
        return own_arrivals, own_leavings

    def time_field(self, vector: numpy.ndarray) -> numpy.ndarray:
        """The time field (ticks x columns): each column's clock."""
        arrivals, leavings = self.events_of(vector)
        lags = self.lags_of(vector)
        extras = self.extras_of(vector)
        ticks = self.context.ticks
        field = numpy.empty((ticks, self.columns))
        for column in range(self.columns):
            own_arrivals, own_leavings = self.column_events(arrivals, leavings, lags, extras, column)
            field[:, column] = node_progress(own_arrivals, own_leavings, ticks, 0.001)
        return field

    def times(self, vector: numpy.ndarray) -> List[numpy.ndarray]:
        """The mouth's clock (node 0's) and, when the word has a port or
        folds, the side's (the last column's)."""
        field = self.time_field(vector)
        clocks = [field[:, 0]]
        if self.columns > self.context.count:
            clocks.append(field[:, -1])
        return clocks

    def time(self, vector: numpy.ndarray) -> numpy.ndarray:
        return self.times(vector)[0]

    def command(self, vector: numpy.ndarray) -> numpy.ndarray:
        field = self.time_field(vector)
        count = self.context.count
        side = field[:, -1] if self.columns > count else None
        return fields.read_postures(self.rows, count, self.space, field[:, :count], side)

    def residual(self, vector: numpy.ndarray) -> numpy.ndarray:
        context = self.context
        if self.cache is not None:
            return self.cache.residual(self, vector)
        if not self.valid_vector(vector):
            return numpy.full(context.word.size, 1e3)
        difference = pyramid_word.equation_error(context, self.command(vector))[context.span]
        if len(difference) < context.word.size or not numpy.all(numpy.isfinite(difference)):
            return numpy.full(context.word.size, 1e3)
        if self.window is not None:
            difference = difference * self.window
        return difference / context.scale


def half_cosine_steps(starts: numpy.ndarray, duration: float, ticks: int, tick_seconds: float) -> numpy.ndarray:
    """One articulator's clock: the sum over its moves of half-cosine steps
    from 0 to 1, each beginning at its start and lasting `duration`. Without
    overlap it is node_progress exactly (phone k held between moves); where
    moves overlap they add, the clock continuous and monotone (with
    node_progress the later move overwrote the earlier: a jump)."""
    now = numpy.arange(ticks, dtype=numpy.float64) * tick_seconds
    field = numpy.zeros(ticks)
    for start in starts:
        share = numpy.clip((now - start) / duration, 0.0, 1.0)
        field += 0.5 * (1.0 - numpy.cos(numpy.pi * share))
    return field


def group_of(context: pyramid_word.Context, column: int) -> str:
    """A command column's articulator group (GROUPS by distance from the
    lips; the port is the velum's, the fold columns the folds')."""
    count = context.count
    if column < count:
        distance = count - 1 - column
        for name, first, last in GROUPS:
            if first <= distance < last:
                return name
    if context.has_velum and column == context.channels["velum"].start:
        return "velum"
    return "folds"


class TraverseProblem(PuppetProblem):
    """The puppet with the traverse belonging to the articulators (owner,
    2026-10-24: "do the per-articulator traverse speed"). Each transition
    keeps one free time, its midpoint m_k (s); column c moves through every
    transition in its articulator's own move time durations[c] (s), starting
    at m_k + lags[c] - durations[c] / 2. vector = [m_0 .. m_{P-2}, the free
    lags, the free durations].

    The distance law is a constant move time per articulator, whatever the
    distance: the mass-spring relation of speech kinematics, peak velocity
    proportional to displacement (CANDIDATE law; for the half-cosine, peak
    velocity pi d / (2 D)). Why: with a free duration a transition, a
    vowel-to-stop move (only the seal node moving, its completion past
    contact unheard) had a duration the sound barely fixes, and those
    durations absorbed the shared lags (the blade's lag −5.9 ms against −8
    from six and nine words, step 56)."""

    def __init__(self, context: pyramid_word.Context, lags: Optional[numpy.ndarray] = None,
                 durations: Optional[numpy.ndarray] = None, secant: float = SECANTS[0],
                 free_lags: Optional[Sequence[int]] = None, free_durations: Optional[Sequence[int]] = None) -> None:
        # Set before the base initializer, which sizes the steps through size().
        self.free_durations = numpy.array([] if free_durations is None else list(free_durations), dtype=numpy.int64)
        super().__init__(context, lags, None, secant, free_lags, None, False)
        self.durations = (numpy.full(self.columns, DURATION_START) if durations is None
                          else numpy.asarray(durations, dtype=numpy.float64))

    def events(self) -> int:
        return self.phones - 1

    def size(self) -> int:
        return self.events() + len(self.free_lags) + len(self.free_durations)

    def midpoints_of(self, vector: numpy.ndarray) -> numpy.ndarray:
        return vector[:self.events()]

    def durations_of(self, vector: numpy.ndarray) -> numpy.ndarray:
        durations = self.durations.copy()
        durations[self.free_durations] = vector[self.events() + len(self.free_lags):]
        return durations

    def vector_of_midpoints(self, midpoints: numpy.ndarray, lags: Optional[numpy.ndarray] = None,
                            durations: Optional[numpy.ndarray] = None) -> numpy.ndarray:
        vector = numpy.empty(self.size())
        events = self.events()
        vector[:events] = midpoints
        split = events + len(self.free_lags)
        vector[events:split] = (self.lags if lags is None else lags)[self.free_lags]
        vector[split:] = (self.durations if durations is None else durations)[self.free_durations]
        return vector

    def valid_vector(self, vector: numpy.ndarray) -> bool:
        """Transitions in order, every move at least one engine step."""
        shortest = 1.0 / sequence_recover.STEP_RATE
        return bool(numpy.all(numpy.diff(self.midpoints_of(vector)) > 0.0)
                    and numpy.all(self.durations_of(vector) >= shortest))

    def time_field(self, vector: numpy.ndarray) -> numpy.ndarray:
        midpoints = self.midpoints_of(vector)
        lags = self.lags_of(vector)
        durations = self.durations_of(vector)
        ticks = self.context.ticks
        field = numpy.empty((ticks, self.columns))
        for column in range(self.columns):
            starts = midpoints + lags[column] - 0.5 * durations[column]
            field[:, column] = half_cosine_steps(starts, durations[column], ticks, 0.001)
        return field


def true_midpoints(word) -> numpy.ndarray:
    """Each transition's midpoint in the plan (sequence_recover.plan)."""
    return 0.5 * (word.true_leavings[:-1] + word.true_arrivals[1:])


def hinted_midpoints(word) -> numpy.ndarray:
    """The hint for the midpoints: the plan's, each pushed HINT_OFFSET
    later or earlier, alternating, as build_word pushes the hinted events
    (its pushes cancel at a transition's midpoint, so they are applied to
    the midpoint itself)."""
    midpoints = true_midpoints(word)
    for index in range(len(midpoints)):
        sign = 1.0 if (index + 1) % 2 else -1.0
        midpoints[index] += sign * sequence_recover.HINT_OFFSET
    return midpoints


def traverse_truth_durations(context: pyramid_word.Context) -> numpy.ndarray:
    """A legible move-time pattern for the synthetic speaker, one a group:
    the blade fastest (50 ms), the folds 60, the lips 70, the front 80, the
    back 100, the velum 110, the pharynx 120 ms. Drawn for the known-answer
    test of identifiability, not as a claim about real articulators."""
    table = {"blade": 0.050, "folds": 0.060, "lips": 0.070, "front": 0.080, "back": 0.100, "velum": 0.110,
             "pharynx": 0.120}
    columns = context.count + (1 if context.has_velum else 0) + (3 if context.has_folds else 0)
    durations = numpy.empty(columns)
    for column in range(columns):
        durations[column] = table[group_of(context, column)]
    return durations


def retime_truth_traverse(context: pyramid_word.Context, lags: numpy.ndarray, durations: numpy.ndarray) -> None:
    """A truth in the traverse family: the plan's midpoints, the given
    lags and move times; the context's true command, sound and scale
    replaced."""
    problem = TraverseProblem(context, lags, durations)
    context.true_command = problem.command(problem.vector_of_midpoints(true_midpoints(context.word)))
    context.true_in, context.true_out = pyramid_word.render_both(context, context.true_command)
    context.scale = float(numpy.sqrt(numpy.mean(context.true_in[context.span] ** 2)))


def retime_truth(context: pyramid_word.Context, lags: numpy.ndarray, extras: numpy.ndarray) -> None:
    """A truth whose columns move at their own lags and extra move times
    (the pyramid truth, command_of, moves them all on one clock): the
    context's true command and sound replaced, its scale with them. With
    zero lags and extras it reproduces command_of exactly (checked in
    main_friction)."""
    problem = PuppetProblem(context, lags, extras)
    word = context.word
    context.true_command = problem.command(problem.vector_of(word.true_arrivals, word.true_leavings))
    context.true_in, context.true_out = pyramid_word.render_both(context, context.true_command)
    context.scale = float(numpy.sqrt(numpy.mean(context.true_in[context.span] ** 2)))


def traverse_worker_start(specs: Sequence[str], secant: float, cached: bool) -> None:
    """The pool's problem for main_traverse: the truth retimed, the
    traverse held at the truth, the midpoints free."""
    context = pyramid_word.Context(specs, ["areas"])
    lags = friction_truth(context)
    durations = traverse_truth_durations(context)
    retime_truth_traverse(context, lags, durations)
    fields_word._problem = TraverseProblem(context, lags, durations, secant)
    if cached:
        attach_cache(fields_word._problem)


def worker_start(specs: Sequence[str], secant: float, lags: Optional[numpy.ndarray],
                 extras: Optional[numpy.ndarray], free_lags: Optional[Sequence[int]] = None,
                 true_lags: Optional[numpy.ndarray] = None, free_extras: Optional[Sequence[int]] = None,
                 true_extras: Optional[numpy.ndarray] = None, crossings: bool = False,
                 cached: bool = False) -> None:
    """The pool's problem (fields_word.worker_columns reads fields_word._problem)."""
    context = pyramid_word.Context(specs, ["areas"])
    if true_lags is not None:
        retime_truth(context, true_lags, numpy.zeros(len(true_lags)) if true_extras is None else true_extras)
    fields_word._problem = PuppetProblem(context, lags, extras, secant, free_lags, free_extras, crossings)
    if cached:
        attach_cache(fields_word._problem)


def attach_cache(problem: PuppetProblem) -> None:
    """The process's resumable residual (adjoint imports this module, so it
    is imported here, when asked for)."""
    import adjoint
    problem.cache = adjoint.ResidualCache(problem.context)


def event_report(problem: PuppetProblem, vector: numpy.ndarray, label: str) -> None:
    """Each event against the truth, in ms: the telemetry itself."""
    word = problem.context.word
    if isinstance(problem, TraverseProblem):
        found = problem.midpoints_of(vector)
        want = true_midpoints(word)
        print("%s midpoints (ms off the truth): %s" % (label, "; ".join(
            "%s>%s %+.3f" % (word.names[k], word.names[k + 1], 1000.0 * (found[k] - want[k]))
            for k in range(problem.events()))), flush=True)
        return
    arrivals, leavings = problem.events_of(vector)
    parts = []
    for k in range(problem.phones - 1):
        parts.append("%s leaves %+.3f, %s arrives %+.3f" % (
            word.names[k], 1000.0 * (leavings[k] - word.true_leavings[k]),
            word.names[k + 1], 1000.0 * (arrivals[k + 1] - word.true_arrivals[k + 1])))
    print("%s events (ms off the truth): %s" % (label, "; ".join(parts)), flush=True)


def parameter_names(problem: PuppetProblem) -> List[str]:
    word = problem.context.word
    names = []
    if isinstance(problem, TraverseProblem):
        for k in range(problem.events()):
            names.append("%s>%s" % (word.names[k], word.names[k + 1]))
        for column in problem.free_lags:
            names.append("lag%d" % column)
        for column in problem.free_durations:
            names.append("duration%d" % column)
        return names
    for k in range(problem.phones - 1):
        if problem.sealing_node[k] >= 0:
            names.append("%s|%s" % (word.names[k], word.names[k + 1]))   # the seal crossing
            names.append("%s>%s" % (word.names[k], word.names[k + 1]))   # the move's duration
        else:
            names.append("%s-" % word.names[k])
            names.append("%s+" % word.names[k + 1])
    for column in problem.free_lags:
        names.append("lag%d" % column)
    for column in problem.free_extras:
        names.append("extra%d" % column)
    return names


def spectrum(pool, problem: PuppetProblem, vector: numpy.ndarray) -> numpy.ndarray:
    """What the sound determines of the unknowns (Phase C's telemetry): the
    Jacobian's singular values (relative), each singular direction's largest
    weights, and each unknown's precision: the move (us) that would raise
    the residual's sum of squares by the amount it has left, s^2 = |r|^2 /
    (N - p), i.e. s sqrt(diag (J'J)^-1) — the sound fits as well anywhere
    within it. Columns the sound does not hear at all (zero) are named, not
    given a precision. Returns the precisions (s; inf where unheard)."""
    current = problem.residual(vector)
    matrix = fields_word.jacobian(pool, problem, vector, current)
    names = parameter_names(problem)
    heard = numpy.sqrt(numpy.sum(matrix ** 2, axis=0)) > 0.0
    precision = numpy.full(len(vector), numpy.inf)
    if numpy.any(heard):
        left, singular, rows = numpy.linalg.svd(matrix[:, heard], full_matrices=False)
        print("singular values (relative): %s" % " ".join("%.2e" % (value / singular[0]) for value in singular),
              flush=True)
        heard_names = [name for name, flag in zip(names, heard) if flag]
        for index in range(len(singular)):
            order = numpy.argsort(-numpy.abs(rows[index]))[:3]
            print("  direction %d (sigma %.2e): %s" % (
                index, singular[index] / singular[0],
                ", ".join("%s %+.2f" % (heard_names[j], rows[index, j]) for j in order)), flush=True)
        excess = max(len(current) - int(numpy.count_nonzero(heard)), 1)
        spread = numpy.sqrt(numpy.sum(current ** 2) / excess)
        inverse_diagonal = numpy.sum((rows / singular[:, None]) ** 2, axis=0)
        precision[heard] = spread * numpy.sqrt(inverse_diagonal)
    print("precision (us, the sound fits as well within): %s" % "; ".join(
        "%s %s" % (name, "unheard" if not flag else "%.2f" % (1e6 * value))
        for name, flag, value in zip(names, heard, precision)), flush=True)
    return precision


def main() -> None:
    specs = sys.argv[1].split(",")
    crossings = "crossings" in sys.argv[2:]
    cached = "plain" not in sys.argv[2:]   # the resumable residual: bit-identical on all six words (step 50)
    context = pyramid_word.Context(specs, ["areas"])
    word = context.word
    truth = fields.progress(word.true_arrivals, word.true_leavings, context.ticks, 0.001)
    problem = PuppetProblem(context, crossings=crossings)
    if cached:
        attach_cache(problem)
    if crossings:
        print("sealing transitions carried as (crossing, duration): %s" % ", ".join(
            "%s->%s node %d at share %.3f" % (word.names[k], word.names[k + 1], problem.sealing_node[k],
                                             problem.crossing_share[k])
            for k in numpy.flatnonzero(problem.sealing_node >= 0)), flush=True)
    # Check: the truth is in the family (its half-cosine command, read through
    # node_progress, against pyramid_word.command_of).
    true_vector = problem.vector_of(word.true_arrivals, word.true_leavings)
    at_truth = pyramid_word.level_db(problem.residual(true_vector))
    print("%s: %d event times; the truth in the family: equation error %.1f dB" % (
        " ".join(specs), problem.size(), at_truth), flush=True)
    vector = problem.vector_of(word.hinted_arrivals, word.hinted_leavings)
    fields_word.report(problem, vector, pyramid_word.level_db(problem.residual(vector)), "start", truth)
    event_report(problem, vector, "start")
    for secant in SECANTS:
        problem.secant = secant
        problem.steps = problem.difference_steps()
        print("secant %.3g s" % secant, flush=True)
        with multiprocessing.Pool(fields_word.POOL, initializer=worker_start,
                                  initargs=(specs, secant, None, None, None, None, None, None, crossings,
                                            cached)) as pool:
            vector = fields_word.levenberg(pool, problem, vector, truth, None)
            event_report(problem, vector, "secant %.3g s" % secant)
            if secant == SECANTS[-1]:
                spectrum(pool, problem, vector)
    numpy.save(os.path.join(sequence_recover.SCRATCH, "puppet_%s%s.npy" % (
        "_".join(specs).replace("!", ""), "_crossings" if crossings else "")), vector)


def friction_truth(context: pyramid_word.Context) -> numpy.ndarray:
    """A legible per-column lag pattern for the synthetic speaker (Phase
    D's first experiment): the tongue tip's region (the six nodes behind
    the lips' two) leads the shared clock by 8 ms, the lips' two nodes lag
    it by 6 ms, the folds lag by 4 ms, the rest move on the shared clock."""
    count = context.count
    columns = context.count + (1 if context.has_velum else 0) + (3 if context.has_folds else 0)
    lags = numpy.zeros(columns)
    lags[count - 8:count - 2] = -0.008
    lags[count - 2:count] = 0.006
    if context.has_folds:
        lags[context.channels["folds"]] = 0.004
    return lags


def friction_truth_extras(context: pyramid_word.Context) -> numpy.ndarray:
    """The pattern's extra move times: the lips' two nodes 10 ms slower than
    the shared move, the tongue tip's six 6 ms faster, the folds 5 ms
    slower, the rest on the shared move."""
    count = context.count
    columns = context.count + (1 if context.has_velum else 0) + (3 if context.has_folds else 0)
    extras = numpy.zeros(columns)
    extras[count - 8:count - 2] = -0.006
    extras[count - 2:count] = 0.010
    if context.has_folds:
        extras[context.channels["folds"]] = 0.005
    return extras


def main_friction() -> None:
    """Can the sound tell per-column traverse? The truth retimed with
    friction_truth's lags; the inverse starts at the hint with zero lags and
    every column's lag free (the spectrum names the unheard ones)."""
    specs = sys.argv[1].split(",")
    context = pyramid_word.Context(specs, ["areas"])
    word = context.word
    # Check: zero lags and unit stretches reproduce command_of.
    plain = context.true_command.copy()
    retime_truth(context, numpy.zeros(plain.shape[1]), numpy.zeros(plain.shape[1]))
    print("retimed truth with no lags against command_of: max |diff| %.3g" % numpy.abs(context.true_command - plain).max(),
          flush=True)
    true_lags = friction_truth(context)
    retime_truth(context, true_lags, numpy.zeros(len(true_lags)))
    truth = fields.progress(word.true_arrivals, word.true_leavings, context.ticks, 0.001)
    free = list(range(plain.shape[1]))
    problem = PuppetProblem(context, free_lags=free)
    true_vector = problem.vector_of(word.true_arrivals, word.true_leavings, true_lags)
    print("%s: %d unknowns (%d events, %d lags); the lagged truth in the family: %.1f dB" % (
        " ".join(specs), problem.size(), problem.events(), len(free), pyramid_word.level_db(problem.residual(true_vector))),
        flush=True)
    vector = problem.vector_of(word.hinted_arrivals, word.hinted_leavings)
    fields_word.report(problem, vector, pyramid_word.level_db(problem.residual(vector)), "start", truth)
    event_report(problem, vector, "start")
    for secant in SECANTS:
        problem.secant = secant
        problem.steps = problem.difference_steps()
        print("secant %.3g s" % secant, flush=True)
        with multiprocessing.Pool(fields_word.POOL, initializer=worker_start,
                                  initargs=(specs, secant, None, None, free, true_lags)) as pool:
            vector = fields_word.levenberg(pool, problem, vector, truth, None)
            event_report(problem, vector, "secant %.3g s" % secant)
            lags = problem.lags_of(vector)
            print("lags recovered (ms; truth): %s" % " ".join(
                "%d:%+.2f(%+.0f)" % (column, 1000.0 * lags[column], 1000.0 * true_lags[column])
                for column in free if abs(lags[column]) > 1e-5 or true_lags[column] != 0.0), flush=True)
            if secant == SECANTS[-1]:
                spectrum(pool, problem, vector)
    numpy.save(os.path.join(sequence_recover.SCRATCH, "puppet_friction_%s.npy" % "_".join(specs).replace("!", "")),
               vector)


def main_traverse() -> None:
    """The traverse family on one word (known answer): the truth retimed
    with friction_truth's lags and traverse_truth_durations' move times,
    the traverse held at the truth, the midpoints free from the hint."""
    specs = sys.argv[1].split(",")
    cached = "plain" not in sys.argv[2:]
    context = pyramid_word.Context(specs, ["areas"])
    word = context.word
    lags = friction_truth(context)
    durations = traverse_truth_durations(context)
    retime_truth_traverse(context, lags, durations)
    problem = TraverseProblem(context, lags, durations)
    if cached:
        attach_cache(problem)
    truth = fields.progress(word.true_arrivals, word.true_leavings, context.ticks, 0.001)
    print("%s: %d midpoints; the truth in the family %.1f dB" % (
        " ".join(specs), problem.size(), pyramid_word.level_db(problem.residual(
            problem.vector_of_midpoints(true_midpoints(word))))), flush=True)
    vector = problem.vector_of_midpoints(hinted_midpoints(word))
    fields_word.report(problem, vector, pyramid_word.level_db(problem.residual(vector)), "start", truth)
    event_report(problem, vector, "start")
    for secant in SECANTS:
        problem.secant = secant
        problem.steps = problem.difference_steps()
        print("secant %.3g s" % secant, flush=True)
        with multiprocessing.Pool(fields_word.POOL, initializer=traverse_worker_start,
                                  initargs=(specs, secant, cached)) as pool:
            vector = fields_word.levenberg(pool, problem, vector, truth, None)
            event_report(problem, vector, "secant %.3g s" % secant)
            if secant == SECANTS[-1]:
                spectrum(pool, problem, vector)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[2] == "friction":
        main_friction()
    elif len(sys.argv) > 2 and sys.argv[2] == "traverse":
        main_traverse()
    else:
        main()
