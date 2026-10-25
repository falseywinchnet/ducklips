"""The puppet's equation error and Jacobian with their work shared and
checkpointed, and the analytic part of the derivative chain.

The record is NOTES_adjoint.md (the cost measurements, the derivation of
the tangent and adjoint of the step, what was built and what was not).

pyramid_word.equation_error evaluates one command from nothing: the moving
geometry, the peel (observer.implied) and three whole engine renders (one
per closure window for the free folds, one for the forced folds; MEASURED
step 45 of NOTES_pyramid.md: 0.20 of a 0.29-0.50 s residual). This module
evaluates the same numbers to the bit with:

- one engine run per residual: the free renders of the closure windows and
  the forced render are prefixes and continuations of one causal run
  (RenderCache: the engine's own snapshot and resume, engine.Run, as
  engine.Incremental uses them); each window's render stops at the window's
  end, which observer.implied never reads past;
- across the evaluations of one process (a Jacobian's columns), the last
  evaluation as the base of the next: the geometry recurrence, the peel
  (observer.implied_resume) and the engine run resume from their last
  recorded state before the first command tick that differs
  (ResidualCache). An event time changes the command only from its move on,
  so a later event's column costs its tail only.

The derivative of the command with respect to the puppet's unknowns (an
event time, a lag, an extra move time) is analytic through the half-cosine
(command_tangent); the derivative of the engine and the peel with respect
to the command is not built here (NOTES_adjoint.md: CANDIDATE, with its
cost).

Usage (on the Mini, JETBOUND=1 SEALBLEND=1):
    python3 adjoint.py profile  SPECS      # one residual and one column, profiled
    python3 adjoint.py check    SPECS      # the shared / resumed residual against pyramid_word's, bit for bit
    python3 adjoint.py columns  SPECS      # the Jacobian both ways: identity, cosines, timing
    python3 adjoint.py lags     SPECS      # the voice pattern: many lag columns in one process (SHARE=12)
    python3 adjoint.py tangent  SPECS      # command_tangent against central differences of the command
    python3 adjoint.py sweep    SPECS      # the residual's change along one event over difference scales
"""
import cProfile
import multiprocessing
import multiprocessing.pool
import os
import pstats
import resource
import sys
import time
from typing import Callable, List, Optional, Sequence, Tuple

import numpy

import engine
import fields_word
import observer
import posture_recover
import puppet
import pyramid_word
import reverse_core
import sequence_recover

TICK = pyramid_word.TICK
EVERY = 960          # engine steps between run snapshots (engine.Incremental's default: 10 ms)
# CANONICAL=1: the workers' caches in the canonical form (SharedFieldFlow: the
# folds free inside every earlier hidden window, one engine run a residual);
# 0 (default): pyramid_word.equation_error's values to the bit.
CANONICAL = os.environ.get("CANONICAL", "0") == "1"
GEOMETRY_EVERY = reverse_core.GEOMETRY_EVERY
POOL = fields_word.POOL


def cpu_seconds() -> float:
    """This process's CPU time (user and system)."""
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_utime + usage.ru_stime


def same_bytes(first: numpy.ndarray, second: numpy.ndarray) -> bool:
    """Bit-for-bit equality (NaN equals NaN)."""
    a = numpy.ascontiguousarray(first)
    b = numpy.ascontiguousarray(second)
    if a.shape != b.shape or a.dtype != b.dtype:
        return False
    return a.tobytes() == b.tobytes()


def prepare_render(context: pyramid_word.Context, command: numpy.ndarray,
                   forced: Optional[numpy.ndarray]) -> engine.Settings:
    """pyramid_word.render_both's engine settings for a command (the same
    arguments, prepared and not yet run)."""
    word = context.word
    field = pyramid_word.field_parts(context, pyramid_word.physical(context, command))
    settings = engine.prepare((word.true_tubes[0], word.length), word.seconds, seal=word.seal, jet=word.jet,
                              jet_gate=context.gate, forced=forced, field=field)
    return settings


class RenderCache:
    """One engine run kept across the renders of one process. A render whose
    inputs equal the last render's before step k (engine.Settings.divergence:
    the construction compared whole, then every control tick and every
    forced wave, bit for bit) resumes from the last snapshot at or before k
    and keeps the earlier outputs. Exact: the engine's state at a step
    depends only on the construction and the inputs before it
    (bindings/equivalence.py checks engine.Incremental's results byte for
    byte; `check` below checks this cache's residuals). Renders may stop
    short of the word's end (`end`): the run is left there and continued
    by the next render that agrees with it so far."""

    def __init__(self, every: int = EVERY) -> None:
        self.every = int(every)
        self.settings: Optional[engine.Settings] = None
        self.snapshots: List[Tuple[int, engine.Run]] = []   # (step, the run's state before that step), ascending
        self.run: Optional[engine.Run] = None               # the live run, at step self.reached
        self.raw: Optional[numpy.ndarray] = None            # (steps, 2): flow entering, flow radiated
        self.reached = 0        # steps of self.raw written by the last render
        self.reused = 0         # steps the last render took from the one before
        self.resumed_from = 0   # the step the last render resumed from
        self.renders = 0
        self.steps_run = 0      # engine steps advanced, all renders

    def render(self, settings: engine.Settings, end: Optional[int] = None) -> Tuple[numpy.ndarray, bool]:
        """The RAW pairs of steps 0..end - 1 (the word's end by default) and
        whether every output was finite."""
        steps = settings.steps()
        stop = steps if end is None else min(int(end), steps)
        resume = -1
        if self.settings is not None:
            resume = min(self.settings.divergence(settings), self.reached)
        kept: List[Tuple[int, engine.Run]] = []
        for step, run in self.snapshots:
            if step <= resume:
                kept.append((step, run))
        raw = numpy.empty((steps, 2))
        if self.run is not None and resume >= self.reached and self.reached > 0:
            # Everything computed so far agrees: the live run continues.
            start = self.reached
            run = self.run
            run.settings = settings
            raw[:start] = self.raw[:start]
        elif kept:
            start, saved = kept[-1]
            run = saved.snapshot()
            run.settings = settings
            raw[:start] = self.raw[:start]
        else:
            start = 0
            run = engine.Run(settings)
            kept = [(0, run.snapshot())]
        self.resumed_from = start
        self.reused = start
        reached = start
        finite = True
        while finite and reached < stop:
            boundary = min((reached // self.every + 1) * self.every, stop)
            before = reached
            reached, finite = run.advance_into(boundary, raw, None)
            if reached == before:
                raise RuntimeError("the engine did not advance at step %d" % reached)
            self.steps_run += reached - before
            if finite and reached == boundary and reached < steps and reached % self.every == 0:
                kept.append((reached, run.snapshot()))
        self.settings = settings
        self.snapshots = kept
        self.run = run
        self.raw = raw
        self.reached = reached
        self.renders += 1
        pairs = raw[:reached].copy()
        return pairs, finite


def closure_windows(geometry: numpy.ndarray, seal: float, hide: float, count: int) -> List[Tuple[int, int]]:
    """The hidden stretches as observer._implied_run lays them out: for each
    closure (first sealed step, last sealed step, seal node), hidden from
    max(first - node, the previous stretch's end) to min(last + 2, count)."""
    windows: List[Tuple[int, int]] = []
    start = 0
    for first, last, node in observer.closures(geometry, seal, hide):
        hidden_from = max(first - node, start)
        hidden_to = min(last + 2, count)
        windows.append((hidden_from, hidden_to))
        start = hidden_to
    return windows


class SharedFieldFlow:
    """pyramid_word.FieldFlow's interface (observer.implied's free_flow) over
    a RenderCache: the folds' flow entering under the command, forced by the
    given returning waves (NaN: free). observer.implied calls it once per
    closure window with the waves NaN from the window's start on and reads
    the window's steps only; pyramid_word.equation_error calls it once more
    with every revealed wave. Each call renders only to the step the caller
    reads (the window's end; the word's end for the last call), so the
    calls are one run: each continues the last where their inputs agree."""

    def __init__(self, context: pyramid_word.Context, command: numpy.ndarray, cache: RenderCache,
                 windows: Sequence[Tuple[int, int]], canonical: bool) -> None:
        self.context = context
        self.command = command
        self.cache = cache
        self.windows = list(windows)
        # observer._implied_run's call for closure k carries, inside the
        # earlier windows, the returning wave the folds themselves made
        # there (0.5 (p - U rho c / A0) of the recorded forward run), which
        # forces them to within rounding of their free run; the final call
        # (pyramid_word.equation_error) carries NaN there. The two renders
        # therefore differ from the first window's start and cannot be one
        # run. Canonical: every earlier window's waves NaN in every call (the
        # folds free inside a hidden stretch, as the stretch was computed),
        # so the calls continue one run; the result differs from
        # pyramid_word's by the rounding of that forcing (MEASURED in
        # NOTES_adjoint.md). False: pyramid_word's values to the bit.
        self.canonical = canonical

    def end_for(self, forced: numpy.ndarray) -> int:
        """The step the caller reads to: the end of the window whose start
        is the first step from which every wave is NaN (the word's end when
        no window starts there)."""
        count = len(forced)
        finite = numpy.isfinite(forced)
        if numpy.all(finite):
            return count
        last_finite = -1
        positions = numpy.flatnonzero(finite)
        if len(positions):
            last_finite = int(positions[-1])
        hidden_from = last_finite + 1
        for window_from, window_to in self.windows:
            if window_from == hidden_from:
                return window_to
        return count

    def __call__(self, forced: numpy.ndarray, observed: Optional[numpy.ndarray] = None,
                 seconds: Optional[float] = None) -> numpy.ndarray:
        end = self.end_for(forced)
        if self.canonical and end < len(forced):
            waves = numpy.array(forced, dtype=numpy.float64)
            for window_from, window_to in self.windows:
                if window_to <= end:
                    waves[window_from:window_to] = numpy.nan
            forced = waves
        settings = prepare_render(self.context, self.command, forced)
        pairs, finite = self.cache.render(settings, end)
        flow = numpy.zeros(len(forced))
        size = min(len(pairs), len(forced))
        flow[:size] = pairs[:size, 0]
        return flow


class ChainedGeometry:
    """reverse_core.moving_geometry resumed from the last evaluation: one
    buffer of GEODUMP rows and the recurrence's state before every
    GEOMETRY_EVERY-th step; a schedule equal to the last before step k is
    computed from the last recorded state at or before k into the buffer's
    tail (the rows before it are the last evaluation's and equal, the
    recurrence being causal). The buffer is returned: valid until the next
    resume, writable."""

    def __init__(self) -> None:
        self.rest: Optional[numpy.ndarray] = None
        self.buffer: Optional[numpy.ndarray] = None
        self.snaps: Optional[numpy.ndarray] = None
        self.counters: Optional[numpy.ndarray] = None
        self.resumed_from = 0

    def geometry(self, rest: numpy.ndarray, first_changed: int, viscous: float) -> numpy.ndarray:
        rest = numpy.ascontiguousarray(rest, dtype=numpy.float64)
        steps = rest.shape[0]
        n = rest.shape[1]
        every = GEOMETRY_EVERY
        slots = max((steps + every - 1) // every, 1)
        whole = self.buffer is None or self.buffer.shape != (steps, n, 12) or first_changed <= 0
        if whole:
            self.buffer = numpy.empty((steps, n, 12))
            self.snaps = numpy.zeros((slots, n, reverse_core.GEOMETRY_STATE))
            self.counters = numpy.zeros((slots, 2))
            state = reverse_core._geometry_start(rest)
            counters = numpy.zeros(2)
            if reverse_core.ENGINE_GLIDE_PHASE != 0:
                reverse_core._engine_glide_start(rest, float(viscous), state, counters)
            begin = 0
        else:
            index = min(int(first_changed), steps) // every
            index = min(index, slots - 1)
            begin = index * every
            state = self.snaps[index].copy()
            counters = self.counters[index].copy()
        reverse_core._geometry_steps(rest, reverse_core.material_table(n), float(viscous), self.buffer, state,
                                     counters, begin, steps, every, self.snaps, self.counters)
        self.rest = rest
        self.resumed_from = begin
        return self.buffer


def first_changed_step(previous: Optional[numpy.ndarray], field: numpy.ndarray, steps: int) -> int:
    """The first engine step whose command tick differs from the previous
    evaluation's (bit for bit, every column: the areas reach the geometry
    and the peel, the folds and the port reach the engine); `steps` when
    nothing differs or there is no previous evaluation to compare with
    (then the caller computes whole: 0)."""
    if previous is None or previous.shape != field.shape:
        return 0
    tick = reverse_core.first_difference(field, previous)
    if tick >= field.shape[0]:
        return steps
    return min(tick * TICK, steps)


class ResidualCache:
    """pyramid_word.equation_error for one word, each evaluation resumed
    from the last (geometry, peel, engine run) at the first command tick
    that differs. Holds the last evaluation's inputs and records. One per
    process."""

    def __init__(self, context: pyramid_word.Context, canonical: bool = False) -> None:
        self.context = context
        self.canonical = canonical
        self.field: Optional[numpy.ndarray] = None
        self.geometry = ChainedGeometry()
        self.renders = RenderCache()
        self.implied_record = None
        self.peel_record = None
        self.evaluations = 0
        self.first_changed = 0

    def equation_error(self, command: numpy.ndarray) -> numpy.ndarray:
        """pyramid_word.equation_error's values (HIDDEN == "sound")."""
        context = self.context
        word = context.word
        sound = context.true_out
        steps = len(sound)
        field = pyramid_word.physical(context, command)
        first_changed = first_changed_step(self.field, field, steps)
        self.first_changed = first_changed
        schedule = numpy.repeat(field[:, context.channels["areas"]], TICK, axis=0)[:steps]
        geometry = self.geometry.geometry(schedule, first_changed, posture_recover.VISCOUS)
        end_extra = word.length - context.count * posture_recover.NODE
        gate = sequence_recover.jet_gate_for(word.true_arrivals, word.true_leavings, word.fricatives, steps)
        if word.seal > 0.0:
            windows = closure_windows(geometry, word.seal, sequence_recover.HIDE_AREA, steps)
            free = SharedFieldFlow(context, command, self.renders, windows, self.canonical)
            if self.implied_record is None or first_changed <= 0:
                implied, pressure, revealed, modelled, record = observer.implied_record(
                    geometry, sound, word.seal, posture_recover.HOP, end_extra, free, hide=sequence_recover.HIDE_AREA,
                    jet=word.jet, jet_closing=posture_recover.JET_CLOSING, jet_gate=gate)
            else:
                implied, pressure, revealed, modelled, record = observer.implied_resume(
                    self.implied_record, geometry, sound, free, first_changed, gate)
            self.implied_record = record
        else:
            free = SharedFieldFlow(context, command, self.renders, [], self.canonical)
            if self.peel_record is None or first_changed <= 0:
                implied, pressure, record = reverse_core.peel_moving_record(
                    geometry, sound, hop=posture_recover.HOP, end_extra=end_extra, jet=word.jet, gate=gate)
            else:
                implied, pressure, record = reverse_core.peel_moving_resume(self.peel_record, geometry, sound,
                                                                            first_changed, gate)
            self.peel_record = record
            revealed = numpy.ones(len(implied), dtype=bool)
            modelled = None
        forced = 0.5 * (pressure - implied * reverse_core.RHO_C / geometry[:, 0, 0])
        forced[~revealed] = numpy.nan
        driven = free(forced)
        size = min(len(driven), len(implied))
        difference = driven[:size] - implied[:size]
        hidden = ~revealed[:size]
        if modelled is not None:
            difference[hidden] = modelled[:size][hidden] - sound[:size][hidden]
        if pyramid_word.ONSET_MASK and context.has_folds:
            fixed = getattr(context, "fixed_voicing_mask", None)
            mask = fixed[:size] if fixed is not None else pyramid_word.voicing_mask(context, field, size)
            difference[:size] *= mask
        self.field = field
        self.evaluations += 1
        return difference

    def residual(self, problem: puppet.PuppetProblem, vector: numpy.ndarray) -> numpy.ndarray:
        """puppet.PuppetProblem.residual's values."""
        context = problem.context
        if not problem.valid_vector(vector):
            return numpy.full(context.word.size, 1e3)
        difference = self.equation_error(problem.command(vector))[context.span]
        if len(difference) < context.word.size or not numpy.all(numpy.isfinite(difference)):
            return numpy.full(context.word.size, 1e3)
        if problem.window is not None:
            difference = difference * problem.window
        return difference / context.scale


# ---- the Jacobian's columns through the cache ----

_problem: Optional[puppet.PuppetProblem] = None
_cache: Optional[ResidualCache] = None


def worker_start(specs: Sequence[str], secant: float) -> None:
    """The pool's problem and cache (puppet.worker_start's problem)."""
    global _problem, _cache
    context = pyramid_word.Context(list(specs), ["areas"])
    _problem = puppet.PuppetProblem(context, None, None, secant)
    _cache = ResidualCache(context, CANONICAL)


def column_first_tick(problem: puppet.PuppetProblem, vector: numpy.ndarray, index: int) -> int:
    """The first command tick either of a parameter's two difference points
    changes (the column's cost is the word from here on)."""
    base_field = pyramid_word.physical(problem.context, problem.command(vector))
    step = problem.steps[index]
    first = problem.context.ticks
    for sign in (1.0, -1.0):
        moved = vector.copy()
        moved[index] += sign * step
        moved_field = pyramid_word.physical(problem.context, problem.command(moved))
        first = min(first, reverse_core.first_difference(moved_field, base_field))
    return first


def worker_columns_cached(task: Tuple[numpy.ndarray, List[int]]) -> numpy.ndarray:
    """Central-difference columns over the secant, every residual through
    the process's ResidualCache (the indices in the order given: the caller
    orders them by decreasing first changed tick, so each resumes at its own
    change)."""
    base, indices = task
    problem = _problem
    cache = _cache
    columns = numpy.empty((problem.context.word.size, len(indices)))
    for position in range(len(indices)):
        index = indices[position]
        step = problem.steps[index]
        ahead = base.copy()
        ahead[index] += step
        behind = base.copy()
        behind[index] -= step
        ahead_residual = cache.residual(problem, ahead)
        behind_residual = cache.residual(problem, behind)
        columns[:, position] = (ahead_residual - behind_residual) / (2.0 * step)
    return columns


def worker_columns_plain(task: Tuple[numpy.ndarray, List[int]]) -> numpy.ndarray:
    """fields_word.worker_columns's columns (every residual from nothing)."""
    base, indices = task
    problem = _problem
    columns = numpy.empty((problem.context.word.size, len(indices)))
    for position in range(len(indices)):
        index = indices[position]
        step = problem.steps[index]
        ahead = base.copy()
        ahead[index] += step
        behind = base.copy()
        behind[index] -= step
        columns[:, position] = (problem.residual(ahead) - problem.residual(behind)) / (2.0 * step)
    return columns


def ordered_chunks(problem: puppet.PuppetProblem, vector: numpy.ndarray, workers: int) -> List[List[int]]:
    """The parameters in decreasing order of first changed tick, dealt to
    `workers` contiguous chunks. In decreasing order each evaluation differs
    from the one before it from its own change on (the one before changed
    later), so each column resumes at its own change."""
    firsts: List[Tuple[int, int]] = []
    for index in range(len(vector)):
        firsts.append((-column_first_tick(problem, vector, index), index))
    firsts.sort()
    order: List[int] = []
    for negative_first, index in firsts:
        order.append(index)
    chunks: List[List[int]] = []
    for part in numpy.array_split(numpy.array(order, dtype=numpy.int64), workers):
        if len(part) > 0:
            chunk: List[int] = []
            for value in part:
                chunk.append(int(value))
            chunks.append(chunk)
    return chunks


def jacobian_cached(pool: multiprocessing.pool.Pool, problem: puppet.PuppetProblem, vector: numpy.ndarray,
                    workers: int) -> numpy.ndarray:
    """fields_word.jacobian's matrix through the cache."""
    chunks = ordered_chunks(problem, vector, workers)
    tasks: List[Tuple[numpy.ndarray, List[int]]] = []
    for chunk in chunks:
        tasks.append((vector, chunk))
    parts = pool.map(worker_columns_cached, tasks)
    matrix = numpy.empty((problem.context.word.size, len(vector)))
    for chunk, part in zip(chunks, parts):
        matrix[:, chunk] = part
    return matrix


def jacobian_plain(pool: multiprocessing.pool.Pool, problem: puppet.PuppetProblem, vector: numpy.ndarray,
                   workers: int) -> numpy.ndarray:
    """fields_word.jacobian's matrix as it is computed today (one parameter
    a task, every residual whole)."""
    chunks = ordered_chunks(problem, vector, workers)
    tasks: List[Tuple[numpy.ndarray, List[int]]] = []
    for chunk in chunks:
        tasks.append((vector, chunk))
    parts = pool.map(worker_columns_plain, tasks)
    matrix = numpy.empty((problem.context.word.size, len(vector)))
    for chunk, part in zip(chunks, parts):
        matrix[:, chunk] = part
    return matrix


# ---- the analytic part of the chain: d command / d theta ----

def progress_tangent(arrivals: numpy.ndarray, leavings: numpy.ndarray, ticks: int, tick_seconds: float,
                     d_arrivals: numpy.ndarray, d_leavings: numpy.ndarray) -> numpy.ndarray:
    """d/d theta of puppet.node_progress for one column, given d/d theta of
    its own arrivals and leavings: on a transition k, s = k + (1 - cos(pi
    f)) / 2 with f = (t - L_k) / (A_{k+1} - L_k), so ds/dtheta = (pi / 2)
    sin(pi f) (df/dtheta), df/dtheta = ((f - 1) dL - f dA) / (A - L); zero in
    a hold. Written over the same loops as node_progress, so where a later
    transition overruns an earlier one the later wins here too."""
    now = numpy.arange(ticks, dtype=numpy.float64) * tick_seconds
    tangent = numpy.zeros(ticks)
    for k in range(len(arrivals)):
        tangent[now >= arrivals[k]] = 0.0
    for k in range(len(arrivals) - 1):
        moving = (now > leavings[k]) & (now < arrivals[k + 1])
        duration = arrivals[k + 1] - leavings[k]
        share = (now[moving] - leavings[k]) / duration
        d_share = ((share - 1.0) * d_leavings[k] - share * d_arrivals[k + 1]) / duration
        tangent[moving] = 0.5 * numpy.pi * numpy.sin(numpy.pi * share) * d_share
    return tangent


def column_event_tangent(problem: puppet.PuppetProblem, arrivals: numpy.ndarray, leavings: numpy.ndarray,
                         lags: numpy.ndarray, extras: numpy.ndarray, column: int, index: int
                         ) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """d/d theta_index of column c's own arrivals and leavings
    (puppet.PuppetProblem.column_events): the shared leaving L_k moves
    own_leavings[k] and own_arrivals[k + 1] alike (the move's start, its
    duration fixed) unless the duration is clipped at one engine step; the
    shared arrival A_{k+1} moves own_arrivals[k + 1] when not clipped; a lag
    moves every own event; an extra moves every own arrival after the first
    when not clipped."""
    phones = problem.phones
    d_arrivals = numpy.zeros(phones)
    d_leavings = numpy.zeros(phones)
    shortest = 1.0 / sequence_recover.STEP_RATE
    events = problem.events()
    clipped = numpy.zeros(phones, dtype=bool)
    for k in range(phones - 1):
        clipped[k + 1] = arrivals[k + 1] - leavings[k] + extras[column] <= shortest
    if index < events:
        k = index // 2
        if index % 2 == 0:
            # L_k: own_leavings[k] = L_k + lag; own_arrivals[k + 1] = own_leavings[k] + max(A - L + extra, shortest).
            d_leavings[k] = 1.0
            d_arrivals[k + 1] = 1.0 if clipped[k + 1] else 0.0
        else:
            # A_{k+1}: own_arrivals[k + 1] = own_leavings[k] + (A_{k+1} - L_k + extra) when not clipped.
            d_arrivals[k + 1] = 0.0 if clipped[k + 1] else 1.0
        return d_arrivals, d_leavings
    lag_position = index - events
    if lag_position < len(problem.free_lags):
        if int(problem.free_lags[lag_position]) == column:
            d_leavings[:] = 1.0
            d_arrivals[:] = 1.0
        return d_arrivals, d_leavings
    extra_position = lag_position - len(problem.free_lags)
    if int(problem.free_extras[extra_position]) == column:
        for k in range(1, phones):
            d_arrivals[k] = 0.0 if clipped[k] else 1.0
    return d_arrivals, d_leavings


def command_tangent(problem: puppet.PuppetProblem, vector: numpy.ndarray, index: int) -> numpy.ndarray:
    """d command / d theta_index (ticks x columns), before the physical
    clip: through fields.read_postures with the identity space field, each
    column linear in its clock between targets, d command[:, c] / d s_c =
    rows[upper, c] - rows[lower, c]; the side columns (port, folds) read the
    last column's clock. The physical clip (pyramid_word.physical) is C0:
    where a channel's command lies below its floor the derivative is zero
    (physical_tangent)."""
    context = problem.context
    arrivals, leavings = problem.events_of(vector)
    lags = problem.lags_of(vector)
    extras = problem.extras_of(vector)
    ticks = context.ticks
    count = context.count
    rows = problem.rows
    phones = problem.phones
    tangent = numpy.zeros((ticks, problem.columns))
    field = problem.time_field(vector)
    for column in range(problem.columns):
        own_arrivals, own_leavings = problem.column_events(arrivals, leavings, lags, extras, column)
        d_arrivals, d_leavings = column_event_tangent(problem, arrivals, leavings, lags, extras, column, index)
        if not numpy.any(d_arrivals) and not numpy.any(d_leavings):
            continue
        d_clock = progress_tangent(own_arrivals, own_leavings, ticks, 0.001, d_arrivals, d_leavings)
        clock = field[:, column]
        clamped = numpy.clip(clock, 0.0, phones - 1.0)
        lower = numpy.minimum(numpy.floor(clamped).astype(numpy.int64), max(phones - 2, 0))
        upper = numpy.minimum(lower + 1, phones - 1)
        slope = rows[upper, column] - rows[lower, column]
        tangent[:, column] = slope * d_clock
    # The side columns read the last column's clock (fields.read_postures):
    # with a port or folds, every column past the areas moves with it.
    if problem.columns > count:
        side_column = problem.columns - 1
        own_arrivals, own_leavings = problem.column_events(arrivals, leavings, lags, extras, side_column)
        d_arrivals, d_leavings = column_event_tangent(problem, arrivals, leavings, lags, extras, side_column, index)
        d_clock = progress_tangent(own_arrivals, own_leavings, ticks, 0.001, d_arrivals, d_leavings)
        clock = field[:, side_column]
        clamped = numpy.clip(clock, 0.0, phones - 1.0)
        lower = numpy.minimum(numpy.floor(clamped).astype(numpy.int64), max(phones - 2, 0))
        upper = numpy.minimum(lower + 1, phones - 1)
        for column in range(count, problem.columns):
            slope = rows[upper, column] - rows[lower, column]
            tangent[:, column] = slope * d_clock
    return tangent


def physical_tangent(context: pyramid_word.Context, command: numpy.ndarray, tangent: numpy.ndarray) -> numpy.ndarray:
    """The tangent through pyramid_word.physical: zero where the clip holds
    (an area at or below the sealed area, a port or rest gap at or below
    zero), the tangent elsewhere. The clip is C0: at the crossing the
    derivative jumps between the two."""
    out = tangent.copy()
    areas = context.channels["areas"]
    out[:, areas] = numpy.where(command[:, areas] > sequence_recover.STOP_AREA, tangent[:, areas], 0.0)
    if context.has_velum:
        velum = context.channels["velum"]
        out[:, velum] = numpy.where(command[:, velum] > 0.0, tangent[:, velum], 0.0)
    if context.has_folds:
        folds = context.channels["folds"]
        out[:, folds] = numpy.where(command[:, folds] > 0.0, tangent[:, folds], 0.0)
    return out


# ---- measurements ----

def hinted_vector(problem: puppet.PuppetProblem) -> numpy.ndarray:
    word = problem.context.word
    return problem.vector_of(word.hinted_arrivals, word.hinted_leavings)


def timed(function: Callable, *arguments) -> Tuple[object, float, float]:
    """(result, wall s, cpu s) of one call of `function` on `arguments`."""
    cpu = cpu_seconds()
    wall = time.perf_counter()
    result = function(*arguments)
    return result, time.perf_counter() - wall, cpu_seconds() - cpu


def describe_word(context: pyramid_word.Context) -> None:
    word = context.word
    events: List[str] = []
    for k in range(word.phones - 1):
        events.append("%s-%.3f %s+%.3f" % (word.names[k], word.true_leavings[k], word.names[k + 1],
                                           word.true_arrivals[k + 1]))
    print("%s: %d steps (%.3f s), %d ticks, %d nodes, window [%d, %d), seal %.3g, jet %.3g; true events (s): %s"
          % (" ".join(word.specs), len(context.true_out), word.seconds, context.ticks, context.count, word.start,
             word.start + word.size, word.seal, word.jet, " ".join(events)), flush=True)


def main_profile(specs: Sequence[str]) -> None:
    """One residual and one column (the first event's), profiled."""
    context = pyramid_word.Context(list(specs), ["areas"])
    describe_word(context)
    problem = puppet.PuppetProblem(context)
    vector = hinted_vector(problem)
    problem.residual(vector)   # numba warm
    walls: List[float] = []
    cpus: List[float] = []
    for repeat in range(5):
        result, wall, cpu = timed(problem.residual, vector)
        walls.append(wall)
        cpus.append(cpu)
    walls.sort()
    cpus.sort()
    print("plain residual: wall %.3f s (min %.3f), cpu %.3f s (min %.3f), %.1f dB" % (
        walls[2], walls[0], cpus[2], cpus[0], pyramid_word.level_db(result)), flush=True)
    profiler = cProfile.Profile()
    profiler.enable()
    problem.residual(vector)
    profiler.disable()
    stats = pstats.Stats(profiler)
    stats.sort_stats("cumulative")
    stats.print_stats(28)
    # One column: two residuals.
    step = problem.steps[0]
    ahead = vector.copy()
    ahead[0] += step
    behind = vector.copy()
    behind[0] -= step
    cpu = cpu_seconds()
    wall = time.perf_counter()
    column = (problem.residual(ahead) - problem.residual(behind)) / (2.0 * step)
    print("one column (secant %.3g s): wall %.3f s, cpu %.3f s, |column| %.3e" % (
        step, time.perf_counter() - wall, cpu_seconds() - cpu, float(numpy.linalg.norm(column))), flush=True)


def main_check(specs: Sequence[str]) -> None:
    """The cached residual against pyramid_word's over a sequence of
    commands: the hinted start, then each event moved by each secant (the
    sequence a Jacobian's worker sees), bit for bit; with the steps each
    stage resumed from."""
    context = pyramid_word.Context(list(specs), ["areas"])
    describe_word(context)
    problem = puppet.PuppetProblem(context)
    cache = ResidualCache(context)
    canonical = ResidualCache(context, True)
    vector = hinted_vector(problem)
    problem.residual(vector)
    cache.residual(problem, vector)
    canonical.residual(problem, vector)
    points: List[Tuple[str, numpy.ndarray]] = [("start", vector)]
    for secant in puppet.SECANTS:
        for index in range(problem.size()):
            for sign in (1.0, -1.0):
                moved = vector.copy()
                moved[index] += sign * secant
                points.append(("%s %+.3g" % (puppet.parameter_names(problem)[index], sign * secant), moved))
    for sign in (1.0, -1.0):
        moved = vector.copy()
        moved[:] += sign * 0.0004
        points.append(("all %+.3g" % (sign * 0.0004), moved))
    failures = 0
    plain_total = [0.0, 0.0]
    cached_total = [0.0, 0.0]
    canonical_total = [0.0, 0.0]
    worst_relative = 0.0
    worst_db = 0.0
    for label, point in points:
        plain, plain_wall, plain_cpu = timed(problem.residual, point)
        cached, cached_wall, cached_cpu = timed(cache.residual, problem, point)
        exact_steps = cache.renders.steps_run
        cache.renders.steps_run = 0
        other, other_wall, other_cpu = timed(canonical.residual, problem, point)
        canonical_steps = canonical.renders.steps_run
        canonical.renders.steps_run = 0
        plain_total[0] += plain_wall
        plain_total[1] += plain_cpu
        cached_total[0] += cached_wall
        cached_total[1] += cached_cpu
        canonical_total[0] += other_wall
        canonical_total[1] += other_cpu
        same = same_bytes(plain, cached)
        if not same:
            failures += 1
        relative = float(numpy.linalg.norm(other - plain) / max(numpy.linalg.norm(plain), 1e-300))
        worst_relative = max(worst_relative, relative)
        worst_db = max(worst_db, abs(pyramid_word.level_db(other) - pyramid_word.level_db(plain)))
        print("%-16s exact %s (%d steps run); canonical rel. diff %.2e (%d steps run); first changed step %6d; "
              "geometry from %6d; plain %.3f/%.3f s, exact %.3f/%.3f, canonical %.3f/%.3f (wall/cpu)" % (
                  label, "same" if same else "DIFFERS max %.3e" % float(numpy.max(numpy.abs(plain - cached))),
                  exact_steps, relative, canonical_steps, cache.first_changed, cache.geometry.resumed_from,
                  plain_wall, plain_cpu, cached_wall, cached_cpu, other_wall, other_cpu), flush=True)
    print("%d points, %d differ (exact); canonical: worst relative difference %.2e, worst dB difference %.2e; "
          "totals wall/cpu: plain %.2f/%.2f s, exact %.2f/%.2f s, canonical %.2f/%.2f s" % (
              len(points), failures, worst_relative, worst_db, plain_total[0], plain_total[1], cached_total[0],
              cached_total[1], canonical_total[0], canonical_total[1]), flush=True)


def column_report(plain: numpy.ndarray, cached: numpy.ndarray, names: Sequence[str]) -> None:
    for column in range(plain.shape[1]):
        a = plain[:, column]
        b = cached[:, column]
        norm_a = float(numpy.linalg.norm(a))
        norm_b = float(numpy.linalg.norm(b))
        cosine = float(numpy.dot(a, b) / (norm_a * norm_b)) if norm_a > 0.0 and norm_b > 0.0 else float("nan")
        print("  %-8s cosine %.6f  |cached|/|plain| %.6f  max |diff| %.3e  %s" % (
            names[column], cosine, norm_b / norm_a if norm_a > 0.0 else float("nan"),
            float(numpy.max(numpy.abs(a - b))), "identical" if same_bytes(a, b) else "differ"), flush=True)


def main_columns(specs: Sequence[str]) -> None:
    """The Jacobian at the hinted start, both ways, at each secant: plain
    (fields_word.worker_columns's way) and cached; identity, cosines and
    relative norms per column; wall time with POOL workers, and the
    sequential cost of all columns in one process (the voice problem's
    pattern: many columns a worker)."""
    context = pyramid_word.Context(list(specs), ["areas"])
    describe_word(context)
    problem = puppet.PuppetProblem(context)
    vector = hinted_vector(problem)
    names = puppet.parameter_names(problem)
    workers = min(POOL, problem.size())
    print("cache mode: %s; %d workers" % ("canonical" if CANONICAL else "exact", workers), flush=True)
    for secant in puppet.SECANTS:
        problem.secant = secant
        problem.steps = problem.difference_steps()
        firsts: List[str] = []
        for index in range(problem.size()):
            firsts.append("%s:%d" % (names[index], column_first_tick(problem, vector, index)))
        print("secant %.3g s: first changed tick per column: %s" % (secant, " ".join(firsts)), flush=True)
        with multiprocessing.Pool(workers, initializer=worker_start, initargs=(specs, secant)) as pool:
            # (The first map warms numba and the engine in every worker.)
            jacobian_plain(pool, problem, vector, workers)
            plain, plain_wall, _ = timed(jacobian_plain, pool, problem, vector, workers)
            cached, cached_wall, _ = timed(jacobian_cached, pool, problem, vector, workers)
            cached_again, cached_again_wall, _ = timed(jacobian_cached, pool, problem, vector, workers)
        print("  %d workers: plain %.2f s wall, cached %.2f s (first, cache empty), %.2f s (again, cache warm)" % (
            workers, plain_wall, cached_wall, cached_again_wall), flush=True)
        column_report(plain, cached, names)
        column_report(plain, cached_again, names)
    # The sequential pattern: every column in one process, plain then cached.
    problem.secant = puppet.SECANTS[-1]
    problem.steps = problem.difference_steps()
    worker_start(specs, problem.secant)
    order = ordered_chunks(_problem, vector, 1)[0]
    worker_columns_plain((vector, order[:1]))
    plain_sequential, plain_wall, plain_cpu = timed(worker_columns_plain, (vector, order))
    cached_sequential, cached_wall, cached_cpu = timed(worker_columns_cached, (vector, order))
    print("one process, %d columns (secant %.3g s): plain %.2f s wall / %.2f cpu; cached %.2f s wall / %.2f cpu; "
          "identical %s" % (len(order), problem.secant, plain_wall, plain_cpu, cached_wall, cached_cpu,
                            same_bytes(plain_sequential, cached_sequential)), flush=True)


def main_lags(specs: Sequence[str]) -> None:
    """The voice problem's pattern on one word: every column's lag free
    (puppet.main_friction's family: events + one lag a column), the columns
    of one worker's share evaluated in one process at the finest secant,
    plain and cached; identity and time. A lag moves its column's every
    move, so its first changed tick is the first transition's."""
    context = pyramid_word.Context(list(specs), ["areas"])
    describe_word(context)
    columns = context.count + (1 if context.has_velum else 0) + (3 if context.has_folds else 0)
    free = list(range(columns))
    secant = puppet.SECANTS[-1]
    global _problem, _cache
    _problem = puppet.PuppetProblem(context, None, None, secant, free)
    _cache = ResidualCache(context, CANONICAL)
    problem = _problem
    print("cache mode: %s" % ("canonical" if CANONICAL else "exact"), flush=True)
    vector = hinted_vector(problem)
    names = puppet.parameter_names(problem)
    share = int(os.environ.get("SHARE", "12"))   # columns a worker gets (the voice: ~100 columns over POOL workers)
    order = ordered_chunks(problem, vector, 1)[0]
    chosen = order[:share]
    listed: List[str] = []
    for index in chosen:
        listed.append("%s:%d" % (names[index], column_first_tick(problem, vector, index)))
    print("%d of %d columns, decreasing first changed tick: %s" % (len(chosen), len(order), " ".join(listed)),
          flush=True)
    worker_columns_plain((vector, chosen[:1]))
    plain, plain_wall, plain_cpu = timed(worker_columns_plain, (vector, chosen))
    cached, cached_wall, cached_cpu = timed(worker_columns_cached, (vector, chosen))
    cached_again, again_wall, again_cpu = timed(worker_columns_cached, (vector, chosen))
    print("plain %.2f s wall / %.2f cpu; cached %.2f / %.2f (cache empty at the start); again %.2f / %.2f; "
          "identical %s %s" % (plain_wall, plain_cpu, cached_wall, cached_cpu, again_wall, again_cpu,
                               same_bytes(plain, cached), same_bytes(plain, cached_again)), flush=True)
    chosen_names: List[str] = []
    for index in chosen:
        chosen_names.append(names[index])
    column_report(plain, cached, chosen_names)


def main_tangent(specs: Sequence[str]) -> None:
    """command_tangent (then physical_tangent) against central differences
    of the command itself, per parameter: cosine, relative norm, and the
    ticks where they differ (a kink or the clip's crossing within the
    difference)."""
    context = pyramid_word.Context(list(specs), ["areas"])
    describe_word(context)
    free = list(range(context.count + (3 if context.has_folds else 0) + (1 if context.has_velum else 0)))
    problem = puppet.PuppetProblem(context, free_lags=free, free_extras=free)
    vector = hinted_vector(problem)
    names = puppet.parameter_names(problem)
    h = 1e-7
    base = pyramid_word.physical(context, problem.command(vector))
    worst = 0.0
    zero_columns = 0
    for index in range(problem.size()):
        ahead = vector.copy()
        ahead[index] += h
        behind = vector.copy()
        behind[index] -= h
        difference = (pyramid_word.physical(context, problem.command(ahead))
                      - pyramid_word.physical(context, problem.command(behind))) / (2.0 * h)
        analytic = physical_tangent(context, problem.command(vector), command_tangent(problem, vector, index))
        norm_difference = float(numpy.linalg.norm(difference))
        norm_analytic = float(numpy.linalg.norm(analytic))
        # A column whose targets are equal across the word moves nothing:
        # the analytic tangent is exactly zero and the central difference
        # is the rounding of (1 - s) a + s a against a, divided by 2 h.
        rounding = numpy.finfo(numpy.float64).eps * float(numpy.max(numpy.abs(base))) / h
        if norm_analytic == 0.0 and float(numpy.max(numpy.abs(difference))) <= 4.0 * rounding:
            zero_columns += 1
            continue
        cosine = float(numpy.sum(difference * analytic) / (norm_difference * norm_analytic)) \
            if norm_difference > 0.0 and norm_analytic > 0.0 else float("nan")
        gap = numpy.abs(difference - analytic)
        scale = max(numpy.max(numpy.abs(difference)), 1e-300)
        off = numpy.flatnonzero(numpy.max(gap, axis=1) > 1e-6 * scale)
        worst = max(worst, float(numpy.max(gap)) / scale)
        listed: List[str] = []
        for tick in off[:8]:
            listed.append(str(int(tick)))
        print("  %-8s |difference| %.4e |analytic| %.4e  cosine %.9f  ratio %.9f  ticks off (> 1e-6 rel) %d%s" % (
            names[index], norm_difference, norm_analytic, cosine, norm_analytic / norm_difference, len(off),
            "" if len(off) == 0 else ": %s" % " ".join(listed)), flush=True)
    print("%d columns move nothing (analytic exactly 0, difference below 4 x eps max|command| / h = %.2e); "
          "worst relative gap of the rest %.3e (difference step %.1e s)" % (zero_columns, 4.0 * rounding, worst, h),
          flush=True)


def main_sweep(specs: Sequence[str]) -> None:
    """The residual's change along one event over difference scales from
    1e-7 s to 1e-2 s at the hinted start: |r(theta + d) - r(theta)| / d
    against d, and the central difference's agreement across scales
    (where it is a derivative and where it is a secant over stairs)."""
    context = pyramid_word.Context(list(specs), ["areas"])
    describe_word(context)
    problem = puppet.PuppetProblem(context)
    cache = ResidualCache(context)
    vector = hinted_vector(problem)
    names = puppet.parameter_names(problem)
    current = cache.residual(problem, vector)
    scales = (1e-7, 3e-7, 1e-6, 3e-6, 1.0 / sequence_recover.STEP_RATE, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
    for index in range(problem.size()):
        previous: Optional[numpy.ndarray] = None
        print("%s:" % names[index], flush=True)
        for d in scales:
            ahead = vector.copy()
            ahead[index] += d
            behind = vector.copy()
            behind[index] -= d
            forward = (cache.residual(problem, ahead) - current) / d
            central = (cache.residual(problem, ahead) - cache.residual(problem, behind)) / (2.0 * d)
            agreement = float("nan")
            if previous is not None:
                agreement = float(numpy.dot(central, previous)
                                  / max(numpy.linalg.norm(central) * numpy.linalg.norm(previous), 1e-300))
            print("  d %.1e s: |forward| %.4e |central| %.4e  cosine with the central at the scale before %.4f" % (
                d, float(numpy.linalg.norm(forward)), float(numpy.linalg.norm(central)), agreement), flush=True)
            previous = central


def main() -> None:
    command = sys.argv[1]
    specs = sys.argv[2].split(",")
    if command == "profile":
        main_profile(specs)
    elif command == "check":
        main_check(specs)
    elif command == "columns":
        main_columns(specs)
    elif command == "lags":
        main_lags(specs)
    elif command == "tangent":
        main_tangent(specs)
    elif command == "sweep":
        main_sweep(specs)
    else:
        raise SystemExit("unknown command %s" % command)


if __name__ == "__main__":
    main()
