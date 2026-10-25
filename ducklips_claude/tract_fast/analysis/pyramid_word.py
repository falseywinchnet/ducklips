"""A word's known answer as a smooth control field, and its inverse from the
sound by the pyramid, coarse to fine, each step an exact least squares
(Gauss-Newton). NOTES_pyramid.md has the plan and the measurements.

The truth: each phone's targets (its tube, its fold posture, its velum
port) held from arrival to leaving, each transition a half-cosine move: a
smooth command. A stop's sealed node aims past contact (the vowel's
area there, negated: the lips press beyond touching), and the contact clips
every area at the engine's sealed area (`physical`), so the closure and the
release are the moments the smooth command crosses the seal: events come
out of the field, they are not parameters.

The inverse: the unknown channels' level grids, level by level; the start
the hint tubes (the VTL tongue table) at the hinted times, as
sequence_recover starts; the residual the radiated flow against the true
sound over the word's window.

Usage: pyramid_word.py SPECS [unknown=areas|areas,folds|areas,folds,velum]
"""
import multiprocessing
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

os.environ.setdefault("RENDERER", "engine")
import numpy
import scipy.signal

import engine
import kernels_jet
import observer
import posture_recover
import reverse_core
import pyramid
import sequence_recover

LEVELS_MAX = int(os.environ.get("LEVELS_MAX", "9"))   # levels the inverse goes up to
ITERATIONS = 8          # Gauss-Newton steps a level at most
DAMPING_START = 1e-3    # Levenberg-Marquardt's first damping: lambda in J'J + lambda diag(J'J), dimensionless
                        # (Marquardt's scaling). (Multiplied by the largest squared column norm as well, every
                        # step was over-damped by that norm: "accepted" steps that moved nothing, Papa's
                        # area and time inverses alike. Fixed 2026-10-21.)
DAMPING_LIMIT = 1e12    # damping growth without a gain that ends a level
POOL = int(os.environ.get("POOL", "8"))
EPSILON = numpy.finfo(numpy.float64).eps
TICK = 96               # engine steps a control tick
HIDDEN = os.environ.get("HIDDEN", "sound")   # a closure's stretch compared by the link's flow ("link") or the
                                            # forward sound against the sound ("sound": sealed_core's)


class Context:
    """The word, its truth and the inverse's fixed parts."""

    def __init__(self, specs: Sequence[str], unknown: Sequence[str]) -> None:
        self.word = sequence_recover.build_word(list(specs))
        word = self.word
        self.total = int(word.seconds * sequence_recover.STEP_RATE)
        self.ticks = (self.total + TICK - 1) // TICK
        self.levels = int(numpy.ceil(numpy.log2(self.ticks))) + 1
        self.stack = pyramid.Pyramid(self.ticks, self.levels)
        self.count = word.count
        self.has_velum = len(word.nasals) > 0
        self.has_folds = len(word.voiceless) > 0
        self.true_command = command_of(self, word.true_tubes, word.true_arrivals, word.true_leavings,
                                       word.true_velum)
        self.channels = channel_slices(self)
        self.unknown = [name for name in unknown if name in self.channels]
        self.frozen = sorted(set(node for k, node in word.stops.items() if k not in word.narrowings))
        self.gate = sequence_recover.jet_keyframes(word.true_arrivals, word.true_leavings, word.fricatives)
        self.true_in, self.true_out = render_both(self, self.true_command)
        self.span = slice(word.start, word.start + word.size)
        # The equation error's scale: the true flow entering's swing.
        self.scale = float(numpy.sqrt(numpy.mean(self.true_in[self.span] ** 2)))


def channel_slices(context: Context) -> Dict[str, slice]:
    """The command's columns: the nodes' areas, then the velum, then the folds."""
    slices = {"areas": slice(0, context.count)}
    column = context.count
    if context.has_velum:
        slices["velum"] = slice(column, column + 1)
        column += 1
    if context.has_folds:
        slices["folds"] = slice(column, column + 3)
    return slices


def phone_at(time: float, arrivals: numpy.ndarray, leavings: numpy.ndarray) -> int:
    """The phone whose stretch holds `time`: stretches meet halfway through
    each transition."""
    index = 0
    for k in range(1, len(arrivals)):
        if time >= 0.5 * (leavings[k - 1] + arrivals[k]):
            index = k
    return index


def targets(context: Context, tubes: numpy.ndarray, velum: numpy.ndarray) -> numpy.ndarray:
    """Each phone's target row: its tube (a stop's sealed node aiming past
    contact, the opening vowel's area there negated), its port, its folds."""
    word = context.word
    width = context.count + (1 if context.has_velum else 0) + (3 if context.has_folds else 0)
    rows = numpy.zeros((word.phones, width))
    for k in range(word.phones):
        tube = numpy.asarray(tubes[k], dtype=numpy.float64).copy()
        if k in word.stops:
            node = word.stops[k]
            vowel = numpy.asarray(tubes[sequence_recover.vowel_before(k, word.stops)], dtype=numpy.float64)
            if k not in word.narrowings:
                tube[node] = -vowel[node]
        rows[k, :context.count] = tube
        column = context.count
        if context.has_velum:
            rows[k, column] = velum[k]
            column += 1
        if context.has_folds:
            posture = word.voiceless.get(k, sequence_recover.MODAL_FOLDS)
            rows[k, column:column + 3] = posture
    return rows


def command_of(context: Context, tubes: numpy.ndarray, arrivals: numpy.ndarray, leavings: numpy.ndarray,
               velum: numpy.ndarray) -> numpy.ndarray:
    """The phones' targets held from arrival to leaving, each transition a
    half-cosine move (smooth, never past either end), read at the ticks.
    (Held targets read through one pyramid level rang: /aa iy/'s ringing
    sealed node 15 mid-vowel and the peel through it failed, the truth's own
    equation error -18.8 dB; REJECTED, NOTES_pyramid.md.)"""
    rows = targets(context, tubes, velum)
    command = numpy.empty((context.ticks, rows.shape[1]))
    for tick in range(context.ticks):
        now = tick * TICK / sequence_recover.STEP_RATE
        command[tick] = rows[0]
        for k in range(len(arrivals)):
            if now >= arrivals[k]:
                command[tick] = rows[k]
            if k + 1 < len(arrivals) and leavings[k] < now < arrivals[k + 1]:
                share = (now - leavings[k]) / (arrivals[k + 1] - leavings[k])
                command[tick] = rows[k] + 0.5 * (1.0 - numpy.cos(numpy.pi * share)) * (rows[k + 1] - rows[k])
                break
    return command


def physical(context: Context, command: numpy.ndarray) -> numpy.ndarray:
    """The contact: an area never below the sealed area, a port and a rest
    gap never below zero."""
    field = command.copy()
    areas = context.channels["areas"]
    field[:, areas] = numpy.maximum(field[:, areas], sequence_recover.STOP_AREA)
    if context.has_velum:
        velum = context.channels["velum"]
        field[:, velum] = numpy.maximum(field[:, velum], 0.0)
    if context.has_folds:
        folds = context.channels["folds"]
        field[:, folds] = numpy.maximum(field[:, folds], 0.0)
    return field


def field_parts(context: Context, field: numpy.ndarray) -> Dict[str, numpy.ndarray]:
    parts = {"areas": field[:, context.channels["areas"]]}
    if context.has_velum:
        parts["velum"] = field[:, context.channels["velum"]][:, 0]
    if context.has_folds:
        parts["folds"] = field[:, context.channels["folds"]]
    return parts


def render_both(context: Context, command: numpy.ndarray, forced: Optional[numpy.ndarray] = None,
                seconds: Optional[float] = None) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """The flow entering and the radiated flow of a command (the folds
    forced by `forced` returning waves where finite)."""
    word = context.word
    return engine.render((word.true_tubes[0], word.length), word.seconds if seconds is None else seconds,
                         seal=word.seal, jet=word.jet, jet_gate=context.gate, forced=forced,
                         field=field_parts(context, physical(context, command)))


def render_command(context: Context, command: numpy.ndarray) -> numpy.ndarray:
    """The radiated flow of a command."""
    return render_both(context, command)[1]


class FieldFlow:
    """Our folds' flow entering under a command, forced by given returning
    waves (NaN: free): observer.implied's `free`."""

    def __init__(self, context: Context, command: numpy.ndarray) -> None:
        self.context = context
        self.command = command

    def __call__(self, forced: numpy.ndarray, observed: Optional[numpy.ndarray] = None,
                 seconds: Optional[float] = None) -> numpy.ndarray:
        flow = render_both(self.context, self.command, forced, seconds)[0]
        if len(flow) < len(forced):
            flow = numpy.concatenate([flow, numpy.zeros(len(forced) - len(flow))])
        return flow[:len(forced)]


ONSET_MASK = os.environ.get("ONSETMASK", "0") == "1"
ONSET_COLLISION = os.environ.get("ONSETCOLLISION", "0") == "1"   # the mask until the first collision


def voicing_mask(context: Context, field: numpy.ndarray, size: int) -> numpy.ndarray:
    """Zero where the folds change posture (their lower rest strictly
    between the word's fold targets: abducting or adducting) and for one
    glottal period after (the true flow's own period, read off it), one
    elsewhere. Near the voicing threshold the folds' forced response turns
    with the phase of their first contacts at the scale of microseconds
    (MEASURED, Alfa's /f ah/ onset: non-monotone at ~5e-5 phones, ~3 us), so
    the fold comparison there cannot carry the clock, as a sealed closure
    cannot; the folds' clock there is a choice, not a measurement."""
    folds = field[:, context.channels["folds"]][:, 0]
    targets = numpy.unique(numpy.round(numpy.array(
        [row[context.channels["folds"]][0] for row in targets_rows(context)]), 12))
    moving = numpy.ones(len(folds), dtype=bool)
    for value in targets:
        moving &= numpy.abs(folds - value) > 1e-12
    period = posture_recover.period_of(context.true_in[context.span])
    steps_moving = numpy.repeat(moving, TICK)[:size]
    mask = numpy.ones(size)
    inside = numpy.flatnonzero(steps_moving)
    if len(inside):
        runs = numpy.split(inside, numpy.flatnonzero(numpy.diff(inside) > 1) + 1)
        for run in runs:
            end = run[-1] + 1 + period
            # Only an adduction (the folds' rest falling toward voice) starts
            # an onset; an abduction into a voiceless phone has no cycles to
            # wait for. (Extended at abductions too, the mask ran across the
            # whole voiceless /f/: 40% of Alfa's steps.)
            first_tick, last_tick = run[0] // TICK, min(run[-1] // TICK, len(folds) - 1)
            adducting = folds[last_tick] < folds[first_tick]
            if ONSET_COLLISION and adducting:
                end = max(end, first_collision(context, run[-1] + 1, period) + period)
            mask[run[0]:min(end, size)] = 0.0
    return mask


def first_collision(context: Context, begin: int, period: int, flow: Optional[numpy.ndarray] = None) -> int:
    """The first glottal cycle after `begin` whose minimum falls to the
    closed-phase level: the median of the cycle minima over the stretch
    that follows (until the folds next move or the word ends). Before it the
    oscillation is still building and its phase rides on microseconds; after
    it the limit cycle holds. (The cycles here are the true glottal flow's,
    a known-answer proxy; a real inverse reads them in the flow the sound
    implies through the start's tube.)"""
    flow = context.true_in if flow is None else flow
    limit = min(len(flow), begin + 40 * period)
    minima = []
    for start in range(begin, limit - period, period):
        minima.append((start, float(flow[start:start + period].min())))
    if not minima:
        return begin
    level = float(numpy.median([value for _, value in minima]))
    for start, value in minima:
        if value <= level:
            return start
    return begin


def implied_flow(context: Context, command: numpy.ndarray) -> numpy.ndarray:
    """The glottal flow the sound implies through a command's tube (the
    peel, no seal handling: for words without one)."""
    word = context.word
    steps = len(context.true_out)
    field = physical(context, command)
    schedule = numpy.repeat(field[:, context.channels["areas"]], TICK, axis=0)[:steps]
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    end_extra = word.length - context.count * posture_recover.NODE
    gate = sequence_recover.jet_gate_for(word.true_arrivals, word.true_leavings, word.fricatives, steps)
    implied, _ = reverse_core.peel_moving_glottis(geometry, context.true_out, hop=posture_recover.HOP,
                                                  end_extra=end_extra, jet=word.jet, gate=gate)
    return implied


def targets_rows(context: Context) -> numpy.ndarray:
    word = context.word
    return targets(context, word.true_tubes, word.true_velum)


def equation_parts(context: Context, command: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray,
                                                                      Optional[numpy.ndarray], numpy.ndarray]:
    """The equation error's parts: our folds' flow driven by the implied
    returning wave, the implied flow, where the sound reveals the inside,
    the forward sound modelled inside the hidden stretches (None without a
    seal), and the moving geometry."""
    word = context.word
    sound = context.true_out
    steps = len(sound)
    field = physical(context, command)
    schedule = numpy.repeat(field[:, context.channels["areas"]], TICK, axis=0)[:steps]
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    end_extra = word.length - context.count * posture_recover.NODE
    gate = sequence_recover.jet_gate_for(word.true_arrivals, word.true_leavings, word.fricatives, steps)
    free = FieldFlow(context, command)
    if word.seal > 0.0:
        implied, pressure, revealed, modelled = observer.implied(
            geometry, sound, word.seal, posture_recover.HOP, end_extra, free, hide=sequence_recover.HIDE_AREA,
            jet=word.jet, jet_closing=posture_recover.JET_CLOSING, jet_gate=gate)
    else:
        implied, pressure = reverse_core.peel_moving_glottis(geometry, sound, hop=posture_recover.HOP,
                                                             end_extra=end_extra, jet=word.jet, gate=gate)
        revealed = numpy.ones(len(implied), dtype=bool)
        modelled = None
    forced = 0.5 * (pressure - implied * reverse_core.RHO_C / geometry[:, 0, 0])
    forced[~revealed] = numpy.nan
    driven = free(forced)
    size = min(len(driven), len(implied))
    return driven[:size], implied[:size], revealed[:size], modelled, geometry


def equation_error(context: Context, command: numpy.ndarray) -> numpy.ndarray:
    """The equation error of a command against the true sound
    (sequence_recover.sealed_core's, for a field): the sound peeled back
    through the command's tube to the flow it implies; the folds forced by
    the returning wave it implies; their flow against the implied flow.
    Inside a seal's hidden stretch, where the sound cannot reveal the
    inside, the folds run free and the forward sound is compared with the
    sound."""
    word = context.word
    sound = context.true_out
    field = physical(context, command)
    driven, implied, revealed, modelled, geometry = equation_parts(context, command)
    size = len(driven)
    difference = driven - implied
    hidden = ~revealed[:size]
    if modelled is not None and HIDDEN == "link":
        # Inside a closure's hidden stretch: the flow through the link in
        # front of its seal, read from the sound through the command's front
        # (the seal node held sealed), against zero.
        for first, last, node in observer.closures(geometry, word.seal, sequence_recover.HIDE_AREA):
            window = numpy.zeros(size, dtype=bool)
            window[max(first - node, 0):min(last + 2, size)] = True
            window &= hidden
            held = command.copy()
            held[:, node] = sequence_recover.STOP_AREA
            flow = link_flow(context, held, node)[:size]
            difference[window] = flow[window]
    elif modelled is not None:
        difference[hidden] = modelled[:size][hidden] - sound[:size][hidden]
    if ONSET_MASK and context.has_folds:
        # Fixed during a solve (set from the start field): a mask moving with
        # each candidate's fold clock would put its edges' jumps back into
        # the cost.
        fixed = getattr(context, "fixed_voicing_mask", None)
        mask = fixed[:size] if fixed is not None else voicing_mask(context, field, size)
        difference[:size] *= mask
    return difference


class Unknowns:
    """The unknown channels' grids of levels 1..used as one vector."""

    def __init__(self, context: Context, used: int, base: numpy.ndarray) -> None:
        self.context = context
        self.used = used
        self.base = base
        # The unknown command columns: the unknown channels' less the seal
        # nodes, whose events are read off the sound and held (read_steps).
        self.indices = [column for name in context.unknown
                        for column in range(context.channels[name].start, context.channels[name].stop)
                        if not (name == "areas" and column in context.frozen)]
        self.columns = len(self.indices)

    def size(self) -> int:
        return sum(self.context.stack.sizes[:self.used]) * self.columns

    def command(self, vector: numpy.ndarray) -> numpy.ndarray:
        """The full command: the known channels the truth's, the unknown the
        base command's plus the sum of the levels (corrections to it)."""
        context = self.context
        grids = []
        offset = 0
        for size in context.stack.sizes[:self.used]:
            grids.append(vector[offset:offset + size * self.columns].reshape(size, self.columns))
            offset += size * self.columns
        unknown_field = numpy.zeros((context.ticks, self.columns))
        for read, grid in zip(context.stack.reads, grids):
            unknown_field += read @ grid
        command = context.true_command.copy()
        for name in context.unknown:
            part = context.channels[name]
            command[:, part] = self.base[:, part]
        command[:, self.indices] += unknown_field
        return command

    def vector_of(self, command: numpy.ndarray) -> numpy.ndarray:
        """The unknown channels of `command` in levels 1..used, coarse to
        fine by exact least squares."""
        context = self.context
        unknown_field = command[:, self.indices]
        grids = context.stack.coarse_to_fine(unknown_field, self.used)
        return numpy.concatenate([grid.ravel() for grid in grids[:self.used]])

    def grow(self, vector: numpy.ndarray) -> Tuple["Unknowns", numpy.ndarray]:
        """One level more: the new level's grid zero."""
        larger = Unknowns(self.context, self.used + 1, self.base)
        extra = numpy.zeros(self.context.stack.sizes[self.used] * self.columns)
        return larger, numpy.concatenate([vector, extra])


def residual(context: Context, unknowns: Unknowns, vector: numpy.ndarray) -> numpy.ndarray:
    """The equation error over the word's window against the true flow
    entering's swing. (The radiated flow against the true sound was flat:
    every step of Papa's levels 1 to 7 stayed at 0.0 dB, the voice's cycle
    phase swamping the tube; REJECTED, NOTES_pyramid.md.)"""
    difference = equation_error(context, unknowns.command(vector))[context.span]
    if len(difference) < context.word.size or not numpy.all(numpy.isfinite(difference)):
        return numpy.full(context.word.size, 1e3)
    return difference / context.scale


def level_db(values: numpy.ndarray) -> float:
    return float(10.0 * numpy.log10(numpy.mean(values ** 2) + 1e-300))


_worker: Optional[Context] = None


def worker_start(specs: Sequence[str], unknown: Sequence[str]) -> None:
    global _worker
    _worker = Context(specs, unknown)


def worker_columns(task: Tuple[int, numpy.ndarray, numpy.ndarray, List[int], numpy.ndarray]) -> numpy.ndarray:
    """Forward-difference columns of the Jacobian at `base` (step
    sqrt(eps) of each value's size, at least sqrt(eps))."""
    used, start, base, indices, base_residual = task
    unknowns = Unknowns(_worker, used, start)
    columns = numpy.empty((len(base_residual), len(indices)))
    for position, index in enumerate(indices):
        step = numpy.sqrt(EPSILON) * max(1.0, abs(base[index]))
        moved = base.copy()
        moved[index] += step
        columns[:, position] = (residual(_worker, unknowns, moved) - base_residual) / step
    return columns


def jacobian(pool, unknowns: Unknowns, vector: numpy.ndarray, base_residual: numpy.ndarray) -> numpy.ndarray:
    chunks = numpy.array_split(numpy.arange(len(vector)), POOL * 4)
    tasks = [(unknowns.used, unknowns.base, vector, list(chunk), base_residual) for chunk in chunks if len(chunk) > 0]
    return numpy.concatenate(pool.map(worker_columns, tasks), axis=1)


def levenberg(pool, context: Context, unknowns: Unknowns, vector: numpy.ndarray, current: numpy.ndarray,
              used: int) -> Tuple[numpy.ndarray, numpy.ndarray, float]:
    """Levenberg-Marquardt steps on one level: each step the exact least
    squares of the damped system [J; sqrt(damping) D] (D the columns'
    norms), the damping adapted by the gain ratio (Nielsen 1999: on a gain,
    times max(1/3, 1 - (2 ratio - 1)^3); on a loss, times nu, nu doubling).
    Ends after ITERATIONS Jacobians, or when the damping has grown past
    DAMPING_LIMIT times its start without a gain."""
    cost = level_db(current)
    for iteration in range(ITERATIONS):
        matrix = jacobian(pool, unknowns, vector, current)
        norms = numpy.sqrt(numpy.sum(matrix ** 2, axis=0))
        live = norms > 0.0
        damping = DAMPING_START
        first = damping
        growth = 2.0
        accepted = False
        while not accepted and damping <= DAMPING_LIMIT * first:
            augmented = numpy.vstack([matrix[:, live], numpy.diag(numpy.sqrt(damping) * norms[live])])
            right = numpy.concatenate([-current, numpy.zeros(int(numpy.count_nonzero(live)))])
            step = numpy.zeros(len(vector))
            step[live] = numpy.linalg.lstsq(augmented, right, rcond=None)[0]
            trial = vector + step
            trial_residual = residual(context, unknowns, trial)
            predicted = float(numpy.sum(current ** 2) - numpy.sum((current + matrix @ step) ** 2))
            actual = float(numpy.sum(current ** 2) - numpy.sum(trial_residual ** 2))
            if actual > 0.0:
                ratio = actual / predicted if predicted > 0.0 else 1.0
                vector, current = trial, trial_residual
                cost = level_db(current)
                damping *= max(1.0 / 3.0, 1.0 - (2.0 * ratio - 1.0) ** 3)
                accepted = True
            else:
                damping *= growth
                growth *= 2.0
        report(context, unknowns, vector, cost, "level %d step %d (damping %.3g)%s" % (
            used, iteration, damping,
            "" if accepted else ": no damping lowered it"))
        if not accepted:
            break
    return vector, current, cost


def seal_flow(context: Context, command: numpy.ndarray, node: int) -> numpy.ndarray:
    """The flow through a seal's node implied by the sound alone: the sound
    peeled back through the command's tube from the lips to that node only.
    While the node is sealed nothing passes it, whatever lies behind; at the
    release it jumps (MEASURED, Papa's true field: about 2e-4 through the
    closure, 1487 at the release)."""
    word = context.word
    steps = len(context.true_out)
    field = physical(context, command)
    schedule = numpy.repeat(field[:, context.channels["areas"]], TICK, axis=0)[:steps]
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    front = numpy.ascontiguousarray(geometry[:, node:, :])
    end_extra = word.length - context.count * posture_recover.NODE
    return reverse_core.peel_moving(front, context.true_out, hop=posture_recover.HOP, end_extra=end_extra)


def link_flow(context: Context, command: numpy.ndarray, node: int) -> numpy.ndarray:
    """The flow through the link in front of a seal's node (node to node +
    1) implied by the sound alone, as the engine steps that link. The sound
    is peeled through the front of the tube (node + 1 on) to node + 1's
    pressure p and inflow; the wave arriving from the link is a = (inflow /
    y + p) / 2 (the junction: inflow = y (2 a - p), y the link's admittance);
    the wave leaving into it, p - a, low-passed as the engine's hop
    (shelf = 1 - HOP, pole HOP) is s; the link's flow is y (a(t+1) - s(t))
    (the engine: a(t+1) = s(t) + flow / y; sealed, a(t+1) = s(t), the flow
    zero). MEASURED: the junction inflow so peeled matched the forward
    replica's to 1e-12 on Echo's 19-node front."""
    word = context.word
    steps = len(context.true_out)
    field = physical(context, command)
    schedule = numpy.repeat(field[:, context.channels["areas"]], TICK, axis=0)[:steps]
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    end_extra = word.length - context.count * posture_recover.NODE
    admittance = 0.5 * (geometry[:, node, 0] + geometry[:, node + 1, 0]) / reverse_core.RHO_C
    if node + 1 == context.count - 1:
        # The front is the open end's node alone: the peel's open-end stage
        # on the last two nodes gives its pressure and the wave arriving
        # from the link directly. (A one-node peel took the link's
        # admittance from the end node twice.)
        front = numpy.ascontiguousarray(geometry[:, node:, :])
        pressure = numpy.zeros(steps)
        arriving = numpy.zeros(steps)
        lining = numpy.zeros((4, 2))
        ends = numpy.zeros(2)
        recorded = numpy.zeros((2, observer.STATE_COLUMNS))
        observer._states_open_end(front, numpy.ascontiguousarray(context.true_out), end_extra, pressure, arriving,
                                  lining, ends, 0, steps, -1, recorded)
    else:
        front = numpy.ascontiguousarray(geometry[:, node + 1:, :])
        # With the word's jet (a fricative in front of the seal: before the
        # closure it is the narrowest node and its jet drives the front;
        # left out, Oscar's front peel entered the closure wrong, the error
        # decaying 2e-3 -> 1.5e-5 over ~25 ms).
        gate = sequence_recover.jet_gate_for(word.true_arrivals, word.true_leavings, word.fricatives, steps)
        if word.jet > 0.0:
            # The jet as the engine places it: the narrowest node and its area
            # over the whole tube (as the velar closes it passes under the
            # fricative's narrowing, behind the front), and its noise hashed by
            # the global link index; mapped into the front's numbering.
            none = numpy.zeros(1)
            run = reverse_core.SidePeel(front, context.true_out, posture_recover.HOP, True, end_extra, word.jet,
                                        reverse_core.JET_SEED, -1, none, none, -1, gate)
            whole_node = numpy.full(steps, -1)
            whole_area = numpy.ones(steps)
            reverse_core._narrowest_steps(numpy.ascontiguousarray(geometry), 0, whole_node, whole_area)
            run.narrowest = whole_node
            run.least = whole_area
            run.link_offset = node + 1
            run.jet_closing = posture_recover.JET_CLOSING if word.seal > 0.0 else 0.0
            run.noise = kernels_jet.noise_table(reverse_core.JET_SEED, context.count - 1, steps)
            run.advance(run.count)
            inflow, pressure = run.arriving, run.pressure
        else:
            inflow, pressure = reverse_core.peel_moving_glottis(front, context.true_out, hop=posture_recover.HOP,
                                                                end_extra=end_extra)
        arriving = 0.5 * (inflow / admittance + pressure)
    leaving = pressure - arriving
    smoothed = scipy.signal.lfilter([1.0 - posture_recover.HOP], [1.0, -posture_recover.HOP], leaving)
    flow = numpy.zeros(steps)
    flow[:-1] = admittance[:-1] * (arriving[1:] - smoothed[:-1])
    return flow


def tick_energy(flow: numpy.ndarray) -> numpy.ndarray:
    ticks = len(flow) // TICK
    return numpy.sum(flow[:ticks * TICK].reshape(ticks, TICK) ** 2, axis=1)


def read_events(context: Context, command: numpy.ndarray, arrivals: numpy.ndarray,
                leavings: numpy.ndarray) -> Dict[int, Tuple[int, int]]:
    """Each stop's closure and release ticks read off the sound: the flow
    through the link in front of its seal (link_flow), the command's seal
    node held sealed over the whole word, so the flow so read is zero
    exactly where the true link is sealed (the candidate's own opening does
    not enter it), and jumps elsewhere. Per tick its energy; the closure the
    tick it falls most against the tick before, the release the tick it
    rises most, each sought between the phone before's leaving and the phone
    after's arrival (-1 for a word's first or last phone). No threshold."""
    word = context.word
    events: Dict[int, Tuple[int, int]] = {}
    tiny = numpy.finfo(numpy.float64).tiny
    for k in sorted(word.stops):
        if k in word.narrowings:
            continue
        node = word.stops[k]
        held = command.copy()
        held[:, node] = sequence_recover.STOP_AREA
        energy = tick_energy(link_flow(context, held, node))
        with numpy.errstate(divide="ignore", over="ignore"):
            rise = numpy.log(energy[1:] + tiny) - numpy.log(energy[:-1] + tiny)
        closure = -1
        release = -1
        if k > 0:
            begin, end = int(leavings[k - 1] * 1000.0), int(arrivals[k] * 1000.0) + 1
            closure = begin + 1 + int(numpy.argmin(rise[begin:end]))
        if k + 1 < word.phones:
            begin, end = int(leavings[k] * 1000.0), int(arrivals[k + 1] * 1000.0) + 1
            release = begin + 1 + int(numpy.argmax(rise[begin:end]))
        events[k] = (closure, release)
    return events


def read_steps(context: Context, command: numpy.ndarray, arrivals: numpy.ndarray,
               leavings: numpy.ndarray) -> Dict[int, Tuple[int, int]]:
    """read_events refined to the step: within the tick read and the one
    before it, the step where the link's flow (seal node held sealed) falls
    most (the closure: the first sealed step) or rises most (the release: the
    first open step), against the step before. A link seals when its glided
    area crosses the seal, at any step."""
    word = context.word
    ticks = read_events(context, command, arrivals, leavings)
    tiny = numpy.finfo(numpy.float64).tiny
    steps: Dict[int, Tuple[int, int]] = {}
    for k, (closure, release) in ticks.items():
        node = word.stops[k]
        held = command.copy()
        held[:, node] = sequence_recover.STOP_AREA
        flow = numpy.abs(link_flow(context, held, node))
        with numpy.errstate(divide="ignore"):
            change = numpy.log(flow[1:] + tiny) - numpy.log(flow[:-1] + tiny)
        found = []
        for tick, sign in ((closure, -1.0), (release, 1.0)):
            if tick < 0:
                found.append(-1)
                continue
            begin, end = (tick - 1) * TICK, (tick + 1) * TICK
            found.append(begin + 1 + int(numpy.argmax(sign * change[begin:end])))
        steps[k] = (found[0], found[1])
    return steps


def sealed_span(context: Context, command: numpy.ndarray, node: int, around: int) -> Tuple[int, int]:
    """The engine's sealed steps of a node (its glided area at or below the
    seal) in the run holding step `around`: (first, last + 1), or (-1, -1)."""
    steps = len(context.true_out)
    field = physical(context, command)
    schedule = numpy.repeat(field[:, context.channels["areas"]], TICK, axis=0)[:steps]
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    sealed = geometry[:, node, 0] <= context.word.seal
    if not sealed[min(max(around, 0), steps - 1)]:
        return -1, -1
    first = around
    while first > 0 and sealed[first - 1]:
        first -= 1
    last = around
    while last + 1 < steps and sealed[last + 1]:
        last += 1
    return first, last + 1


def align_events(context: Context, tubes: numpy.ndarray, arrivals: numpy.ndarray, leavings: numpy.ndarray,
                 velum: numpy.ndarray, read: Dict[int, Tuple[int, int]]) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """The times moved so each stop's sealed steps begin at the closure step
    read and end at the release step read: a closure's transition (the phone
    before's leaving and the stop's arrival) or a release's (the stop's
    leaving and the next phone's arrival) shifted together, by bisection on
    the shift (the sealed span's edge moves monotonically with it)."""
    arrivals = numpy.array(arrivals, dtype=numpy.float64)
    leavings = numpy.array(leavings, dtype=numpy.float64)
    word = context.word
    for k, (closure, release) in sorted(read.items()):
        node = word.stops[k]
        middle = (closure + release) // 2 if closure >= 0 else release - 10 * 96  # a release alone: just inside its closure
        for edge, target in ((0, closure), (1, release)):
            if target < 0:
                continue
            # The shifts that put the edge on the target step form a band;
            # the middle of the band (its edges are where the residual
            # jumps a step: MEASURED, a 1e-9 move at a band's edge changed
            # Papa's residual by 0.9 of its scale).
            def edge_at(shift: float) -> int:
                trial_arrivals = arrivals.copy()
                trial_leavings = leavings.copy()
                if edge == 0:
                    trial_leavings[k - 1] += shift
                    trial_arrivals[k] += shift
                else:
                    trial_leavings[k] += shift
                    trial_arrivals[k + 1] += shift
                command = command_of(context, tubes, trial_arrivals, trial_leavings, velum)
                position = sealed_span(context, command, node, middle)[edge]
                return position if position >= 0 else -1

            def boundary(below: int) -> float:
                low, high = -0.03, 0.03
                for _ in range(40):
                    shift = 0.5 * (low + high)
                    position = edge_at(shift)
                    if position < 0 or position < below:
                        low = shift
                    else:
                        high = shift
                return high

            shift = 0.5 * (boundary(target) + boundary(target + 1))
            if edge == 0:
                leavings[k - 1] += shift
                arrivals[k] += shift
            else:
                leavings[k] += shift
                arrivals[k + 1] += shift
    return arrivals, leavings


def crossings(track: numpy.ndarray, level: float) -> List[float]:
    below = track < level
    found: List[float] = []
    for tick in numpy.nonzero(below[1:] != below[:-1])[0]:
        a, b = track[tick], track[tick + 1]
        found.append(float(tick + (level - a) / (b - a)))
    return found


def report(context: Context, unknowns: Unknowns, vector: numpy.ndarray, cost: float, label: str) -> None:
    command = unknowns.command(vector)
    field = physical(context, command)
    truth = physical(context, context.true_command)
    areas = context.channels["areas"]
    area_error = float(numpy.sqrt(numpy.mean((field[:, areas] - truth[:, areas]) ** 2)))
    events = []
    for node in sorted(set(context.word.stops.values())):
        got = crossings(field[:, node], context.word.seal)
        want = crossings(truth[:, node], context.word.seal)
        if len(got) == len(want):
            events.append("node %d off %s ms" % (node, " ".join("%+.2f" % (g - w) for g, w in zip(got, want))))
        else:
            events.append("node %d: %d crossings, true %d" % (node, len(got), len(want)))
    print("%s: sound %.1f dB; area rms %.4f cm2; %s" % (label, cost, area_error, "; ".join(events)), flush=True)


def main() -> None:
    specs = sys.argv[1].split(",")
    unknown = sys.argv[2].split(",") if len(sys.argv) > 2 else ["areas"]
    context = Context(specs, unknown)
    word = context.word
    truth = physical(context, context.true_command)
    print("%s: %d ticks, levels %s; unknown %s; seal crossings (ms) %s" % (
        " ".join(specs), context.ticks, context.stack.sizes, context.unknown,
        {node: ["%.2f" % c for c in crossings(truth[:, node], word.seal)]
         for node in sorted(set(word.stops.values()))}), flush=True)
    # The start: the hint tubes (a stop the hint of the vowel it is made
    # from, its node sealed) at the hinted times.
    hint_tubes = numpy.zeros_like(word.true_tubes)
    for k, name in enumerate(word.names):
        if k not in word.stops:
            hint_tubes[k] = posture_recover.nodes(word.hints[name], word.count)[0]
    for k in sorted(word.stops):
        hint_tubes[k] = hint_tubes[sequence_recover.vowel_before(k, word.stops)]
        hint_tubes[k, word.stops[k]] = word.narrowings.get(k, sequence_recover.STOP_AREA)
    start_velum = word.hinted_velum if context.has_velum else word.true_velum
    start = command_of(context, hint_tubes, word.hinted_arrivals, word.hinted_leavings, start_velum)
    # The stops' closures and releases read off the sound (read_steps) and
    # the start's seals placed on them (align_events).
    read = read_steps(context, start, word.hinted_arrivals, word.hinted_leavings)
    arrivals, leavings = align_events(context, hint_tubes, word.hinted_arrivals, word.hinted_leavings, start_velum,
                                      read)
    start = command_of(context, hint_tubes, arrivals, leavings, start_velum)
    print("events read (steps): %s; times aligned: arrivals %s, leavings %s" % (
        read, numpy.round(arrivals, 5), numpy.round(leavings, 5)), flush=True)
    # The field is the sum of its levels (the paper's form): the unknown
    # columns start as the start's projection onto level 1 and gain detail
    # level by level; only the seal nodes, whose events are read, keep the
    # start's full resolution. (A full-resolution base stalled: its fine
    # structure lay outside the coarse levels' reach. MEASURED: from the
    # hint offset projected onto level 5, LM went -11.3 -> -28.7 dB in 6
    # steps; from the full-resolution hint start it stalled at -10.2 dB.)
    base = start.copy()
    unknowns = Unknowns(context, 1, base)
    base[:, unknowns.indices] = 0.0
    vector = unknowns.vector_of(start)
    with multiprocessing.Pool(POOL, initializer=worker_start, initargs=(specs, unknown)) as pool:
        for used in range(1, LEVELS_MAX + 1):
            if used > 1:
                unknowns, vector = unknowns.grow(vector)
            current = residual(context, unknowns, vector)
            cost = level_db(current)
            report(context, unknowns, vector, cost, "level %d (%d unknowns) start" % (used, unknowns.size()))
            vector, current, cost = levenberg(pool, context, unknowns, vector, current, used)
    numpy.savez(os.path.join(sequence_recover.SCRATCH, "pyramid_%s.npz" % "_".join(specs).replace("!", "")),
                vector=vector, used=unknowns.used, truth=context.true_command,
                command=unknowns.command(vector))


if __name__ == "__main__":
    main()
