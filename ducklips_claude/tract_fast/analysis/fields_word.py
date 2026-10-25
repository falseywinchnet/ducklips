"""The inverse in two fields (fields.py): the word's motion as postures read
through a space field and a time field, the fields' pyramid coefficients
recovered from the sound coarse to fine, each step an exact least squares
(Levenberg-Marquardt, pyramid_word.levenberg's damping).

Stage `time`: the time field alone (one value per tick, shared by every
node; the running sum of a clipped rate, every pyramid level at once), the
postures the truth's, the space field the identity; the start the hinted
phone times. NOTES_pyramid.md has the measurements.

Usage: fields_word.py SPECS [stage]
"""
import multiprocessing
import os
import sys
from typing import List, Optional, Sequence, Tuple

os.environ.setdefault("RENDERER", "engine")
import numpy
import scipy.linalg
import scipy.stats.qmc

import fields
import posture_recover
import pyramid
import pyramid_word
import tube_recover
import sequence_recover

LEVELS_MAX = int(os.environ.get("LEVELS_MAX", "11"))
ITERATIONS = int(os.environ.get("ITERATIONS", "8"))
POOL = int(os.environ.get("POOL", "8"))
CHECK_EVENTS = os.environ.get("CHECK_EVENTS", "1") == "1"
SECANT_TICKS = [0.0]   # the graduated secant's width in ticks (0: the default quantum); set by main()
GRADUATED = os.environ.get("GRADUATED", "0") == "1"
TICK_STEPS = 96   # engine steps a control tick (pyramid_word.TICK)
EPSILON = numpy.finfo(numpy.float64).eps


class TimeProblem:
    """The time field as the running sum of a clipped rate: s(tick) = sum
    over the ticks before of max(r, 0), r (phones a tick) carried by the
    pyramid's levels 1..used. Progress may pause, never reverse: each seal
    is crossed once each way. (A time field carried directly crossed Papa's
    lip seal five times where the word has three.) The truth's rate (its
    holds zero) is exact at the finest level."""

    def __init__(self, context: pyramid_word.Context, used: int) -> None:
        self.context = context
        self.used = used
        word = context.word
        self.rows = pyramid_word.targets(context, word.true_tubes, word.true_velum)
        self.space = fields.identity_space(context.ticks, context.count)
        self.steps = self.difference_steps()
        self.knots: List[Tuple[float, float]] = []

    def clocks(self) -> int:
        """The time fields: the mouth's, and the folds' and port's own when
        the word has any (they are their own articulators: read through the
        mouth's clock, Papa stalled at -38.6 dB with 69% of it at /ah/'s
        voicing onset after the second release)."""
        return 2 if self.rows.shape[1] > self.context.count else 1

    def size(self) -> int:
        return self.clocks() * sum(self.context.stack.sizes[:self.used])

    def reads(self) -> numpy.ndarray:
        return numpy.concatenate(self.context.stack.reads[:self.used], axis=1)

    def difference_steps(self) -> numpy.ndarray:
        """Each coefficient's difference step: the move that shifts its time
        field by at most one step's quantum (time_quantum) anywhere."""
        effect = numpy.cumsum(numpy.abs(self.reads()), axis=0)
        return numpy.tile(time_quantum(self.context) / numpy.max(effect, axis=0), self.clocks())

    def times(self, vector: numpy.ndarray) -> List[numpy.ndarray]:
        reads = self.reads()
        width = reads.shape[1]
        found = []
        for clock in range(self.clocks()):
            rate = numpy.maximum(reads @ vector[clock * width:(clock + 1) * width], 0.0)
            total = numpy.concatenate([[0.0], numpy.cumsum(rate)[:-1]])
            found.append(self.through_knots(total) if clock == 0 and self.knots else total)
        if getattr(self, "fixed_folds", None) is not None and len(found) > 1:
            # The folds' clock held where the read onsets put it.
            found[1] = self.fixed_folds
        return found

    def kinds(self) -> numpy.ndarray:
        """Per coefficient, 0 for a clipped rate's (scaled by its quantum),
        1 for a smooth one's (scaled by its column)."""
        return numpy.zeros(self.size())

    def through_knots(self, total: numpy.ndarray) -> numpy.ndarray:
        """The mouth clock through its knots (the events read off the sound:
        (tick, value) where a seal is crossed): between two knots, from one
        value to the next in proportion to the clipped rate's running sum
        `total`; before the first and after the last, the running sum's own
        course. Every event held exactly. (Held to first order instead, by
        linear constraints, a step moved Papa's second release 92 steps: the
        clipping's kink in the holds is not linear.)"""
        ticks = numpy.arange(len(total), dtype=numpy.float64)
        at = numpy.interp([tick for tick, _ in self.knots], ticks, total)
        field = numpy.empty_like(total)
        first_tick, first_value = self.knots[0]
        before = ticks <= first_tick
        field[before] = first_value - (at[0] - total[before])
        for index in range(len(self.knots) - 1):
            (tick_a, value_a), (tick_b, value_b) = self.knots[index], self.knots[index + 1]
            inside = (ticks > tick_a) & (ticks <= tick_b)
            span = at[index + 1] - at[index]
            share = (total[inside] - at[index]) / span if span > 0.0 else (ticks[inside] - tick_a) / (tick_b - tick_a)
            field[inside] = value_a + (value_b - value_a) * share
        last_tick, last_value = self.knots[-1]
        after = ticks > last_tick
        field[after] = last_value + (total[after] - at[-1])
        return field

    def time(self, vector: numpy.ndarray) -> numpy.ndarray:
        return self.times(vector)[0]

    def vector_of(self, time: numpy.ndarray) -> numpy.ndarray:
        rate = numpy.concatenate([numpy.diff(time), [0.0]])
        grids = self.context.stack.coarse_to_fine(rate[:, None], self.used)
        one = numpy.concatenate([grid[:, 0] for grid in grids[:self.used]])
        return numpy.tile(one, self.clocks())

    def command(self, vector: numpy.ndarray) -> numpy.ndarray:
        clocks = self.times(vector)
        side = clocks[1] if len(clocks) > 1 else None
        return fields.read_postures(self.rows, self.context.count, self.space,
                                    fields.same_time(clocks[0], self.context.count), side)

    def residual(self, vector: numpy.ndarray) -> numpy.ndarray:
        context = self.context
        difference = pyramid_word.equation_error(context, self.command(vector))[context.span]
        if len(difference) < context.word.size or not numpy.all(numpy.isfinite(difference)):
            return numpy.full(context.word.size, 1e3)
        if getattr(self, "window", None) is not None:
            # A stage's stretch of the sound only (staged solve).
            difference = difference * self.window
        return difference / context.scale


SPACE_X_LEVELS = 3     # the space field's levels along the tube (1 + 2 + 3 cells: a constriction may slide)
SPACE_T_LEVELS = int(os.environ.get("SPACE_T_LEVELS", "6"))   # and in time (37 cells over Papa's 0.54 s);
                       # 0: no space field (the postures and clocks first: the most flexible field otherwise
                       # takes what the others cannot yet explain, MEASURED: |d| 10 nodes on the first step of
                       # a truth that has none)


def hint_tubes(context: pyramid_word.Context) -> numpy.ndarray:
    """The start's tubes: each vowel's hint (the VTL tongue table), a stop
    its vowel's with its node at the seal area (as sequence_recover)."""
    word = context.word
    tubes = numpy.zeros_like(word.true_tubes)
    for k, name in enumerate(word.names):
        if k not in word.stops:
            tubes[k] = posture_recover.nodes(word.hints[name], word.count)[0]
    for k in sorted(word.stops):
        tubes[k] = tubes[sequence_recover.vowel_before(k, word.stops)]
        tubes[k, word.stops[k]] = word.narrowings.get(k, sequence_recover.STOP_AREA)
    return tubes


class FullProblem(TimeProblem):
    """The clocks, the vowels' postures (their tubes; a stop's follows from
    its vowel, its seal node aiming past contact) and the space field
    (xi = x + d(x, t), d on SPACE_X_LEVELS levels along the tube times
    SPACE_T_LEVELS in time). The knots' values follow the postures: each
    knot is (tick, seal node, the phone below the crossing), its value where
    the current postures cross the seal there."""

    def __init__(self, context: pyramid_word.Context, used: int, tubes: numpy.ndarray) -> None:
        self.base_tubes = numpy.array(tubes, dtype=numpy.float64)
        word = context.word
        self.vowels = [k for k in range(word.phones) if k not in word.stops]
        # The coarsest levels of the tube's full pyramid (1, 2, 3 cells).
        across = pyramid.Pyramid(context.count, int(numpy.ceil(numpy.log2(context.count))) + 1)
        self.space_x = numpy.concatenate(across.reads[:SPACE_X_LEVELS], axis=1)
        self.space_t = (numpy.concatenate(context.stack.reads[:SPACE_T_LEVELS], axis=1) if SPACE_T_LEVELS > 0
                        else numpy.zeros((context.ticks, 0)))
        self.event_knots: List[Tuple[float, int, int]] = []
        TimeProblem.__init__(self, context, used)

    def rate_size(self) -> int:
        return self.clocks() * sum(self.context.stack.sizes[:self.used])

    def size(self) -> int:
        return (self.rate_size() + len(self.vowels) * self.context.count +
                self.space_x.shape[1] * self.space_t.shape[1])

    def kinds(self) -> numpy.ndarray:
        kinds = numpy.ones(self.size())
        kinds[:self.rate_size()] = 0.0
        return kinds

    def difference_steps(self) -> numpy.ndarray:
        rates = TimeProblem.difference_steps(self)
        smooth = numpy.full(self.size() - len(rates), numpy.cbrt(EPSILON))
        return numpy.concatenate([rates, smooth])

    def tubes_of(self, vector: numpy.ndarray) -> numpy.ndarray:
        word = self.context.word
        tubes = self.base_tubes.copy()
        offset = self.rate_size()
        for k in self.vowels:
            tubes[k] = vector[offset:offset + self.context.count]
            offset += self.context.count
        for k in sorted(word.stops):
            tubes[k] = tubes[sequence_recover.vowel_before(k, word.stops)]
            tubes[k, word.stops[k]] = word.narrowings.get(k, sequence_recover.STOP_AREA)
        return tubes

    def rows_of(self, vector: numpy.ndarray) -> numpy.ndarray:
        return pyramid_word.targets(self.context, self.tubes_of(vector), self.context.word.true_velum)

    def space_of(self, vector: numpy.ndarray) -> numpy.ndarray:
        """xi = x + w(t) d(x, t) m(x): w = 4 f (1 - f), f the mouth clock's
        fraction (zero at every target: a posture is the tube's shape at its
        target, the space field only how the tube travels between them); m
        zero at the seal nodes (a stop's place is the transcript's).
        (Unweighted, a still space field is a different posture: the two
        traded, Papa's sound -8.1 -> -16.1 dB while its vowel areas went
        0.93 -> 1.56 cm2 off, and xi moved the seal node and its events.)"""
        offset = self.rate_size() + len(self.vowels) * self.context.count
        grid = vector[offset:].reshape(self.space_t.shape[1], self.space_x.shape[1])
        clock = self.times(vector)[0]
        fraction = clock - numpy.floor(clock)
        weight = 4.0 * fraction * (1.0 - fraction)
        mask = numpy.ones(self.context.count)
        for node in seal_nodes(self.context):
            mask[node] = 0.0
        displacement = (self.space_t @ grid @ self.space_x.T) * weight[:, None] * mask[None, :]
        return fields.identity_space(self.context.ticks, self.context.count) + displacement

    def times(self, vector: numpy.ndarray) -> List[numpy.ndarray]:
        if self.event_knots:
            rows = self.rows_of(vector)
            self.knots = []
            for tick, node, lower in self.event_knots:
                # The crossing between phone `lower` and the next.
                value = fields.level_crossing(rows, node, lower + 1, self.context.word.seal)[0]
                self.knots.append((tick, value if numpy.isfinite(value) else lower + 0.5))
        reduced = numpy.concatenate([vector[:self.rate_size()]])
        return TimeProblem.times(self, reduced)

    def vector_of_start(self, time: numpy.ndarray) -> numpy.ndarray:
        rates = TimeProblem.vector_of(self, time)
        postures = numpy.concatenate([self.base_tubes[k] for k in self.vowels])
        space = numpy.zeros(self.space_x.shape[1] * self.space_t.shape[1])
        return numpy.concatenate([rates, postures, space])

    def command(self, vector: numpy.ndarray) -> numpy.ndarray:
        clocks = self.times(vector)
        side = clocks[1] if len(clocks) > 1 else None
        return fields.read_postures(self.rows_of(vector), self.context.count, self.space_of(vector),
                                    fields.same_time(clocks[0], self.context.count), side)


class BodyProblem(FullProblem):
    """FullProblem with each vowel's posture in the articulatory model's
    coordinates (posture_recover.FIELDS: jaw, tongue body, tongue tip, lip
    aperture, protrusion, lip width, drop), not its tube's free areas: the
    body cannot make most tubes that sound alike. (With 46 free areas a
    vowel, Papa's timing and shapes traded: sound -6.9 -> -21.3 dB while
    the timing went 0.064 -> 0.141 phones off and the vowels stayed about
    1 cm2 off.)"""

    def __init__(self, context: pyramid_word.Context, used: int, postures: dict) -> None:
        self.postures = postures
        word = context.word
        tubes = numpy.zeros_like(word.true_tubes)
        FullProblem.__init__(self, context, used, tubes)
        start = [numpy.array([postures[word.names[k]][field] for field in posture_recover.FIELDS])
                 for k in self.vowels]
        self.start_parameters = start
        self.base_tubes = self.tubes_from(numpy.concatenate(start))

    def posture_width(self) -> int:
        return len(posture_recover.FIELDS)

    def size(self) -> int:
        return (self.rate_size() + len(self.vowels) * self.posture_width() +
                self.space_x.shape[1] * self.space_t.shape[1])

    def difference_steps(self) -> numpy.ndarray:
        rates = TimeProblem.difference_steps(self)
        bodies = numpy.full(len(self.vowels) * self.posture_width(), numpy.cbrt(EPSILON))
        space = numpy.full(self.space_x.shape[1] * self.space_t.shape[1], numpy.cbrt(EPSILON))
        return numpy.concatenate([rates, bodies, space])

    def tubes_from(self, parameters: numpy.ndarray) -> numpy.ndarray:
        word = self.context.word
        tubes = numpy.zeros_like(word.true_tubes)
        width = self.posture_width()
        for index, k in enumerate(self.vowels):
            posture = dict(zip(posture_recover.FIELDS, parameters[index * width:(index + 1) * width]))
            tubes[k] = posture_recover.nodes(posture, self.context.count)[0]
        for k in sorted(word.stops):
            tubes[k] = tubes[sequence_recover.vowel_before(k, word.stops)]
            tubes[k, word.stops[k]] = word.narrowings.get(k, sequence_recover.STOP_AREA)
        return tubes

    def tubes_of(self, vector: numpy.ndarray) -> numpy.ndarray:
        offset = self.rate_size()
        return self.tubes_from(vector[offset:offset + len(self.vowels) * self.posture_width()])

    def space_of(self, vector: numpy.ndarray) -> numpy.ndarray:
        offset = self.rate_size() + len(self.vowels) * self.posture_width()
        grid = vector[offset:].reshape(self.space_t.shape[1], self.space_x.shape[1])
        clock = self.times(vector)[0]
        fraction = clock - numpy.floor(clock)
        weight = 4.0 * fraction * (1.0 - fraction)
        mask = numpy.ones(self.context.count)
        for node in seal_nodes(self.context):
            mask[node] = 0.0
        displacement = (self.space_t @ grid @ self.space_x.T) * weight[:, None] * mask[None, :]
        return fields.identity_space(self.context.ticks, self.context.count) + displacement

    def vector_of_start(self, time: numpy.ndarray) -> numpy.ndarray:
        rates = TimeProblem.vector_of(self, time)
        space = numpy.zeros(self.space_x.shape[1] * self.space_t.shape[1])
        return numpy.concatenate([rates, numpy.concatenate(self.start_parameters), space])


def time_quantum(context: pyramid_word.Context) -> float:
    """The time field's change over one engine step, on average: (phones -
    1) / steps. An event lands on a step, so the sound is a staircase in the
    time field below this; a difference over it is a secant across the
    stairs. (At sqrt(eps) the columns were blind or a step's jump over 1e-8:
    Papa's levels 2-9 gained nothing. Over a whole tick the secant blurred
    the last of it: Papa stalled at -38.6 dB, 69% of it just after the
    second release.)"""
    if SECANT_TICKS[0] > 0.0:
        # Graduated: the secant spans SECANT_TICKS[0] ticks of clock, narrowed
        # by main() from ten ticks to one to one engine step as each rung
        # stops gaining (the onset's roughness reaches past a tick while the
        # cost descends at the scale of the start's error).
        return (context.word.phones - 1.0) * SECANT_TICKS[0] / context.ticks
    if os.environ.get("QUANTUM") == "tick":
        # Over a tick instead: a secant wide enough to average the voicing
        # onset's ripples (~5e-5 phones; the step's quantum, 5.8e-5, samples
        # one ripple: Alfa's Jacobian disagreed with its function, cosines
        # 0.54 / -0.48 / -0.39).
        return (context.word.phones - 1.0) / context.ticks
    return (context.word.phones - 1.0) / (context.ticks * pyramid_word.TICK)


_problem: Optional[TimeProblem] = None


def worker_start(specs: Sequence[str], knots: Optional[list] = None, stage: str = "time", used: int = 1,
                 window: Optional[numpy.ndarray] = None, voicing: Optional[numpy.ndarray] = None,
                 fixed_folds: Optional[numpy.ndarray] = None) -> None:
    global _problem
    context = pyramid_word.Context(specs, ["areas"])
    if voicing is not None:
        context.fixed_voicing_mask = voicing
    if stage == "full":
        _problem = FullProblem(context, used, hint_tubes(context))
        _problem.event_knots = list(knots) if knots else []
    elif stage in ("body", "staged"):
        _problem = BodyProblem(context, used, context.word.hints)
        _problem.event_knots = list(knots) if knots else []
    else:
        _problem = TimeProblem(context, used)
        _problem.knots = list(knots) if knots else []
    _problem.window = window
    _problem.fixed_folds = fixed_folds


def graduated_worker_start(specs, knots, width, voicing, fixed_folds) -> None:
    SECANT_TICKS[0] = width
    worker_start(specs, knots, "time", 1, None, voicing, fixed_folds)


def worker_columns(task: Tuple[int, numpy.ndarray, List[int], numpy.ndarray]) -> numpy.ndarray:
    """Central-difference columns over the time quantum."""
    used, base, indices, base_residual = task
    if _problem.used != used:
        _problem.used = used
        _problem.steps = _problem.difference_steps()
    columns = numpy.empty((len(base_residual), len(indices)))
    for position, index in enumerate(indices):
        step = _problem.steps[index]
        ahead = base.copy()
        ahead[index] += step
        behind = base.copy()
        behind[index] -= step
        columns[:, position] = (_problem.residual(ahead) - _problem.residual(behind)) / (2.0 * step)
    return columns


def jacobian(pool, problem: TimeProblem, vector: numpy.ndarray, base_residual: numpy.ndarray) -> numpy.ndarray:
    """The Jacobian's columns (zero for a stage's frozen coefficients)."""
    frozen = getattr(problem, "frozen", None)
    active = numpy.flatnonzero(~frozen) if frozen is not None else numpy.arange(len(vector))
    chunks = [chunk for chunk in numpy.array_split(active, POOL * 2) if len(chunk) > 0]
    tasks = [(problem.used, vector, list(chunk), base_residual) for chunk in chunks]
    matrix = numpy.zeros((len(base_residual), len(vector)))
    matrix[:, active] = numpy.concatenate(pool.map(worker_columns, tasks), axis=1)
    return matrix


def report(problem: TimeProblem, vector: numpy.ndarray, cost: float, label: str, truth: numpy.ndarray) -> None:
    context = problem.context
    word = context.word
    time = problem.time(vector)
    field = pyramid_word.physical(context, problem.command(vector))
    true_field = pyramid_word.physical(context, context.true_command)
    events = []
    for node in sorted(set(node for k, node in word.stops.items() if k not in word.narrowings)):
        got = pyramid_word.crossings(field[:, node], word.seal)
        want = pyramid_word.crossings(true_field[:, node], word.seal)
        events.append("node %d %s" % (node, " ".join("%+.3f" % (g - w) for g, w in zip(got, want))
                                      if len(got) == len(want) else "%d of %d crossings" % (len(got), len(want))))
    # Clamped to the word's phones: past the last (or before the first) the
    # field reads the same posture and cannot be heard.
    last = context.word.phones - 1.0
    time = numpy.clip(time, 0.0, last)
    truth = numpy.clip(truth, 0.0, last)
    clocks = problem.times(vector)
    side = ""
    if len(clocks) > 1:
        folds = numpy.clip(clocks[1], 0.0, last)
        side = "; folds' clock off by %.4f max, %.4f rms" % (numpy.abs(folds - truth).max(),
                                                           numpy.sqrt(numpy.mean((folds - truth) ** 2)))
    if isinstance(problem, FullProblem):
        tubes = problem.tubes_of(vector)
        error = tubes[problem.vowels] - word.true_tubes[problem.vowels]
        displacement = problem.space_of(vector) - fields.identity_space(context.ticks, context.count)
        side += "; vowel areas off %.4f rms cm2; space field |d| max %.3f nodes" % (
            numpy.sqrt(numpy.mean(error ** 2)), numpy.abs(displacement).max())
    print("%s: sound %.1f dB; time field off the truth by %.4f max, %.4f rms (phones)%s; %s" % (
        label, cost, numpy.abs(time - truth).max(), numpy.sqrt(numpy.mean((time - truth) ** 2)), side,
        "; ".join(events)), flush=True)


def nearest(candidates: List[float], target: float) -> float:
    """The candidate closest to `target` (the target itself with none)."""
    best = target
    distance = float("inf")
    for candidate in candidates:
        if abs(candidate - target) < distance:
            best = candidate
            distance = abs(candidate - target)
    return best


def seal_nodes(context: pyramid_word.Context) -> List[int]:
    word = context.word
    return sorted(set(node for k, node in word.stops.items() if k not in word.narrowings))


def event_rows(problem: TimeProblem, vector: numpy.ndarray) -> numpy.ndarray:
    """The constraints that hold every seal crossing where it is, to first
    order: at each crossing tick t (fractional) of a seal node's command,
    the mouth clock's derivative there, d s(t) / d coefficients (s the
    running sum of the clipped rate: only ticks with a positive rate
    count). Zero for the folds' clock."""
    context = problem.context
    reads = problem.reads()
    width = reads.shape[1]
    rate = reads @ vector[:width]
    moving = (rate > 0.0).astype(numpy.float64)
    partial = numpy.vstack([numpy.zeros((1, width)), numpy.cumsum(reads * moving[:, None], axis=0)[:-1]])
    field = pyramid_word.physical(context, problem.command(vector))
    rows = []
    for node in seal_nodes(context):
        for crossing in pyramid_word.crossings(field[:, node], context.word.seal):
            tick = int(numpy.floor(crossing))
            share = crossing - tick
            row = (1.0 - share) * partial[tick] + share * partial[min(tick + 1, len(partial) - 1)]
            full = numpy.zeros(len(vector))
            full[:width] = row
            rows.append(full)
    return numpy.array(rows) if rows else numpy.zeros((0, len(vector)))


def event_steps(problem: TimeProblem, vector: numpy.ndarray, read: dict) -> dict:
    """The engine's sealed steps of each stop under this field, in the run
    holding the read events (pyramid_word.sealed_span)."""
    context = problem.context
    command = problem.command(vector)
    found = {}
    for k, (closure, release) in read.items():
        node = context.word.stops[k]
        middle = (closure + release) // 2 if closure >= 0 else release - 10 * 96  # a release alone: just inside its closure
        first, last = pyramid_word.sealed_span(context, command, node, middle)
        found[k] = (first if closure >= 0 else -1, last if release >= 0 else -1)
    return found


ACTION = os.environ.get("ACTION", "0") == "1"
HALVINGS = 12   # halvings of the action step before it is given up
AGREEMENT = 0.25   # a step is kept only when its actual gain is at least this share of the gain its linear
                   # model predicted (Nocedal and Wright's eta: outside that, the model does not hold and the
                   # step shrinks). (Keeping any gain at all kept a step that slid along the sound-alike
                   # valley for 0.1 dB: Papa's timing 0.020 -> 0.130 phones rms off, the folds' 0 -> 0.26.)


def action_matrix(problem: TimeProblem) -> numpy.ndarray:
    """The clocks' action as a matrix on the coefficients: each clock's rate
    change from tick to tick (its jerk of progress), rows for every clock;
    ||action_matrix @ vector||^2 is the action (the clipping ignored: the
    rate where it moves)."""
    reads = problem.reads()
    width = reads.shape[1]
    change = numpy.diff(reads, axis=0)
    blocks = []
    for clock in range(problem.clocks()):
        block = numpy.zeros((change.shape[0], problem.size()))
        block[:, clock * width:(clock + 1) * width] = change
        blocks.append(block)
    return numpy.vstack(blocks)


def levenberg(pool, problem: TimeProblem, vector: numpy.ndarray, truth: numpy.ndarray,
              read: Optional[dict] = None) -> numpy.ndarray:
    """LM steps; with `read` (the events read off the sound), each step in
    the null space of event_rows (the seal crossings held to first order)
    and kept only if the engine's sealed steps stay the ones read."""
    current = problem.residual(vector)
    cost = pyramid_word.level_db(current)
    for iteration in range(ITERATIONS):
        matrix = jacobian(pool, problem, vector, current)
        # The step in each coefficient's own quantum (problem.steps: the move
        # that shifts its time field by one engine step's worth), damped
        # uniformly in those units; the damping starts at tau times the
        # largest diagonal of the scaled J'J (Madsen, Nielsen and
        # Tingleff). (Marquardt's scaling by column norms left clipped
        # coefficients, columns near 7e-7, undamped: steps of norm 1e6 that
        # raised the cost.)
        # A column a failed render touched carries no information.
        matrix[:, ~numpy.all(numpy.isfinite(matrix), axis=0)] = 0.0
        norms = numpy.sqrt(numpy.sum(matrix ** 2, axis=0))
        kinds = problem.kinds()
        scale = numpy.where(kinds == 0.0, problem.steps,
                            numpy.where(norms > 0.0, 1.0 / numpy.maximum(norms, 1e-300), 0.0))
        # Each kind brought to unit mean diagonal, so one damping means the
        # same for both. (The rates' diagonals, far larger, set the damping
        # for all: the postures were over-damped and barely moved, Papa's
        # vowel areas 1.01 -> 0.99 cm2 in four steps.)
        for kind in (0.0, 1.0):
            members = (kinds == kind) & (norms > 0.0)
            if numpy.any(members):
                mean = float(numpy.mean((norms[members] * scale[members]) ** 2))
                if mean > 0.0:
                    scale[kinds == kind] /= numpy.sqrt(mean)
        scaled = matrix * scale[None, :]
        diagonal = numpy.sum(scaled ** 2, axis=0)
        live = diagonal > 0.0
        if not numpy.any(live):
            break
        damping = pyramid_word.DAMPING_START
        first = damping
        growth = 2.0
        accepted = False
        basis = numpy.eye(int(numpy.count_nonzero(live)))
        if read is not None and not problem.knots:
            constraints = (event_rows(problem, vector) * problem.steps[None, :])[:, live]
            if len(constraints):
                basis = scipy.linalg.null_space(constraints)
        reduced = scaled[:, live] @ basis
        left, singular, rows = numpy.linalg.svd(reduced, full_matrices=False)
        projected = left.T @ (-current)
        effort = None
        if ACTION:
            effort = action_matrix(problem)
            effort_z = effort[:, live] @ (scale[live][:, None] * basis)
        while not accepted and damping <= pyramid_word.DAMPING_LIMIT * first:
            # The damped Gauss-Newton step in the scaled coordinates (the
            # SVD form of the augmented least squares).
            z = rows.T @ (singular / (singular ** 2 + damping) * projected)
            step = numpy.zeros(len(vector))
            step[live] = scale[live] * (basis @ z)
            if read is not None and CHECK_EVENTS and event_steps(problem, vector + step, read) != read:
                damping *= growth
                growth *= 2.0
                continue
            trial = problem.residual(vector + step)
            predicted = float(numpy.sum(current ** 2) - numpy.sum((current + matrix @ step) ** 2))
            actual = float(numpy.sum(current ** 2) - numpy.sum(trial ** 2))
            ratio = actual / predicted if predicted > 0.0 else 0.0
            if actual > 0.0 and ratio >= AGREEMENT:
                vector, current = vector + step, trial
                cost = pyramid_word.level_db(current)
                if effort is not None:
                    # Lexicographic: the sound first; then, where the sound
                    # does not see (sigma^2 below this damping), the step
                    # that lowers the clocks' action (least jerk of
                    # progress), halved until the sound is no worse.
                    # The directions the sound does not hear at all (the
                    # numerical null space: sigma below sqrt(eps) of the
                    # largest; the observability map's -6000 dB stretches).
                    # (Below the damping instead, the sound rose by a hair
                    # in every step and none was kept.)
                    unseen = rows[singular <= numpy.sqrt(EPSILON) * singular[0]].T
                    if unseen.shape[1] > 0:
                        target = -(effort @ vector)
                        weights = numpy.linalg.lstsq(effort_z @ unseen, target, rcond=None)[0]
                        quiet = numpy.zeros(len(vector))
                        quiet[live] = scale[live] * (basis @ (unseen @ weights))
                        share = 1.0
                        for _ in range(HALVINGS):
                            calm = problem.residual(vector + share * quiet)
                            if numpy.sum(calm ** 2) <= numpy.sum(current ** 2):
                                before = float(numpy.sum((effort @ vector) ** 2))
                                vector, current = vector + share * quiet, calm
                                cost = pyramid_word.level_db(current)
                                print("  action step (share %.3g, %d unseen directions): action %.3e -> %.3e" % (
                                    share, unseen.shape[1], before, float(numpy.sum((effort @ vector) ** 2))),
                                    flush=True)
                                break
                            share *= 0.5
                damping *= max(1.0 / 3.0, 1.0 - (2.0 * ratio - 1.0) ** 3)
                accepted = True
            else:
                damping *= growth
                growth *= 2.0
        report(problem, vector, cost, "level %d step %d%s" % (problem.used, iteration,
                                                               "" if accepted else ": no damping lowered it"), truth)
        if not accepted:
            break
    return vector


def main() -> None:
    specs = sys.argv[1].split(",")
    context = pyramid_word.Context(specs, ["areas"])
    word = context.word
    truth = fields.progress(word.true_arrivals, word.true_leavings, context.ticks, 0.001)
    start = fields.progress(word.hinted_arrivals, word.hinted_leavings, context.ticks, 0.001)
    # Every level at once, as the paper trains its pyramid (solved level by
    # level, the coarse fits walked away from the truth: a coarse field's
    # best fit is not the truth's coarse part).
    problem = TimeProblem(context, min(LEVELS_MAX, context.levels))
    read = None
    if os.environ.get("EVENTS", "1") == "1" and seal_nodes(context):
        # The stops' events read off the sound, and the start's seals put on
        # them (pyramid_word.read_steps, align_events).
        hinted = problem.command(problem.vector_of(start))
        read = pyramid_word.read_steps(context, hinted, word.hinted_arrivals, word.hinted_leavings)
        if os.environ.get("EVENT_STEPS"):
            # Events recovered elsewhere (outside_in.py), as phone:closure:
            # release, -1 for an event not given: only these are knotted.
            read = {}
            for item in os.environ["EVENT_STEPS"].split(";"):
                phone, closure_step, release_step = (int(x) for x in item.split(":"))
                read[phone] = (closure_step, release_step)
        arrivals, leavings = pyramid_word.align_events(context, word.true_tubes, word.hinted_arrivals,
                                                       word.hinted_leavings, word.true_velum, read)
        start = fields.progress(arrivals, leavings, context.ticks, 0.001)
        print("events read (steps): %s" % read, flush=True)
        knots = []
        aligned = fields.read_postures(problem.rows, context.count, problem.space,
                                       fields.same_time(start, context.count))
        # Each knot's value exactly where the postures cross the seal (the
        # node's area linear in the time field between targets), its tick
        # where the aligned start crosses it. (Interpolated from the start
        # between ticks, the truth through the knots was 1.3e-2 off.)
        given = [step for pair in read.values() for step in pair if step >= 0]
        for node in seal_nodes(context):
            for crossing in pyramid_word.crossings(aligned[:, node], word.seal):
                if os.environ.get("EVENT_STEPS") and not any(abs(crossing * TICK_STEPS - step) < 5 * TICK_STEPS
                                                             for step in given):
                    continue
                near = float(numpy.interp(crossing, numpy.arange(context.ticks), start))
                phone = int(round(near))
                candidates = [value for value in fields.level_crossing(problem.rows, node, phone, word.seal)
                              if numpy.isfinite(value)]
                value = nearest(candidates, near)
                knots.append((crossing, value))
        problem.knots = sorted(knots)
        print("knots (tick, time field): %s" % ["(%.3f, %.4f)" % knot for knot in problem.knots], flush=True)
    vector = problem.vector_of(start)
    if pyramid_word.ONSET_MASK and context.has_folds:
        start_field = pyramid_word.physical(context, problem.command(vector))
        context.fixed_voicing_mask = pyramid_word.voicing_mask(context, start_field, len(context.true_out))
        print("voicing mask (from the start): %d of %d steps masked" % (
            int(numpy.sum(context.fixed_voicing_mask == 0.0)), len(context.fixed_voicing_mask)), flush=True)
    saved = os.path.join(pyramid_word.sequence_recover.SCRATCH,
                         "fields_time_%s.npy" % "_".join(specs).replace("!", ""))
    if os.environ.get("RESUME") == "1" and os.path.exists(saved):
        loaded = numpy.load(saved)
        # A one-clock solution starts both clocks.
        vector = loaded if len(loaded) == len(vector) else numpy.tile(loaded, len(vector) // len(loaded))
    print("%s: time field from the hinted times (max %.3f phones off the truth); levels %s, %d unknowns" % (
        " ".join(specs), numpy.abs(problem.time(vector) - truth).max(), context.stack.sizes[:problem.used],
        problem.size()), flush=True)
    fixed_folds = None
    if os.environ.get("FOLDS_TRUE") == "1":
        # Diagnostic only: the folds' clock held at the truth, the mouth's
        # solved (does the mouth converge once the folds' clock is right?).
        fixed_folds = truth.copy()
        problem.fixed_folds = fixed_folds
        problem.frozen = numpy.zeros(problem.size(), dtype=bool)
        problem.frozen[problem.reads().shape[1]:] = True
    with multiprocessing.Pool(POOL, initializer=worker_start,
                              initargs=(specs, problem.knots, "time", 1, None,
                                        getattr(context, "fixed_voicing_mask", None), fixed_folds)) as pool:
        report(problem, vector, pyramid_word.level_db(problem.residual(vector)), "start", truth)
        if read is not None and event_steps(problem, vector, read) != read:
            print("the start's sealed steps %s are not the ones read" % event_steps(problem, vector, read))
        if GRADUATED:
            pool.close()
            pool.join()
            for width in (10.0, 1.0, 1.0 / pyramid_word.TICK):
                SECANT_TICKS[0] = width
                problem.steps = problem.difference_steps()
                print("secant over %.4g ticks" % width, flush=True)
                with multiprocessing.Pool(POOL, initializer=graduated_worker_start,
                                          initargs=(specs, problem.knots, width,
                                                    getattr(context, "fixed_voicing_mask", None),
                                                    fixed_folds)) as rung_pool:
                    vector = levenberg(rung_pool, problem, vector, truth, read)
        else:
            vector = levenberg(pool, problem, vector, truth, read)
    numpy.save(os.path.join(pyramid_word.sequence_recover.SCRATCH,
                            "fields_time_%s.npy" % "_".join(specs).replace("!", "")), vector)


def main_full() -> None:
    """Stage `full`: the clocks, the vowels' postures and the space field
    from the hints (FullProblem), the events read off the sound and held as
    knots."""
    specs = sys.argv[1].split(",")
    context = pyramid_word.Context(specs, ["areas"])
    word = context.word
    truth = fields.progress(word.true_arrivals, word.true_leavings, context.ticks, 0.001)
    start = fields.progress(word.hinted_arrivals, word.hinted_leavings, context.ticks, 0.001)
    stage = sys.argv[2]
    if stage == "body":
        problem = BodyProblem(context, min(LEVELS_MAX, context.levels), word.hints)
        tubes = problem.base_tubes
    else:
        tubes = hint_tubes(context)
        problem = FullProblem(context, min(LEVELS_MAX, context.levels), tubes)
    hinted = problem.command(problem.vector_of_start(start))
    read = pyramid_word.read_steps(context, hinted, word.hinted_arrivals, word.hinted_leavings)
    arrivals, leavings = pyramid_word.align_events(context, tubes, word.hinted_arrivals, word.hinted_leavings,
                                                   word.true_velum, read)
    start = fields.progress(arrivals, leavings, context.ticks, 0.001)
    vector = problem.vector_of_start(start)
    # The unclipped command: the clock is linear between targets, so its
    # crossing between ticks is where the unclipped area crosses (the
    # clipped area's interpolation put the knots off: the truth's hold sound
    # at -59.2 dB instead of the floor).
    aligned = problem.command(vector)
    for node in seal_nodes(context):
        for crossing in pyramid_word.crossings(aligned[:, node], word.seal):
            lower = int(numpy.floor(numpy.interp(crossing, numpy.arange(context.ticks), start)))
            problem.event_knots.append((crossing, node, lower))
    problem.event_knots.sort()
    print("%s: full stage; events read (steps) %s; knots (tick, node, phone below) %s; %d unknowns" % (
        " ".join(specs), read, ["(%.3f, %d, %d)" % knot for knot in problem.event_knots], problem.size()),
        flush=True)
    saved = os.path.join(pyramid_word.sequence_recover.SCRATCH,
                         "fields_%s_%s.npy" % (stage, "_".join(specs).replace("!", "")))
    if os.environ.get("RESUME") == "1" and os.path.exists(saved):
        vector = numpy.load(saved)
    with multiprocessing.Pool(POOL, initializer=worker_start,
                              initargs=(specs, problem.event_knots, stage, problem.used)) as pool:
        report(problem, vector, pyramid_word.level_db(problem.residual(vector)), "start", truth)
        vector = levenberg(pool, problem, vector, truth, read)
        numpy.save(saved, vector)


FREEZE_FIELDS = [field for field in os.environ.get("FREEZE_FIELDS", "").split(",") if field]
# Articulatory parameters held at the start and reported unrecovered: those
# the sound cannot separate (MEASURED at Papa's truth: lip aperture against
# lip width, 4e-8 of the strongest direction, and protrusion, 2.2e-7; the
# sound hears the lip opening's area and barely its length).


def unheard(problem: "BodyProblem") -> numpy.ndarray:
    """The coefficients of the FREEZE_FIELDS parameters, every vowel's."""
    mask = numpy.zeros(problem.size(), dtype=bool)
    offset = problem.rate_size()
    width = problem.posture_width()
    for index in range(len(problem.vowels)):
        for j, field in enumerate(posture_recover.FIELDS):
            if field in FREEZE_FIELDS:
                mask[offset + index * width + j] = True
    return mask


def hold_window(context: pyramid_word.Context, arrivals: numpy.ndarray, leavings: numpy.ndarray) -> numpy.ndarray:
    """The central half of each vowel's hold (by the aligned hinted times),
    over the compared stretch: where the tube is still and the clocks do not
    matter, so the sound there pins the postures."""
    word = context.word
    window = numpy.zeros(word.size)
    for k in range(word.phones):
        if k in word.stops:
            continue
        quarter = 0.25 * (leavings[k] - arrivals[k])
        begin = int((arrivals[k] + quarter) * sequence_recover.STEP_RATE) - context.span.start
        end = int((leavings[k] - quarter) * sequence_recover.STEP_RATE) - context.span.start
        window[max(begin, 0):max(min(end, word.size), 0)] = 1.0
    return window


SURVEY = int(os.environ.get("SURVEY", "256"))   # postures a vowel's survey scores (a Sobol sequence)


def vowel_window(context: pyramid_word.Context, arrivals: numpy.ndarray, leavings: numpy.ndarray,
                 k: int) -> numpy.ndarray:
    """The central half of vowel k's hold alone."""
    word = context.word
    window = numpy.zeros(word.size)
    quarter = 0.25 * (leavings[k] - arrivals[k])
    begin = int((arrivals[k] + quarter) * sequence_recover.STEP_RATE) - context.span.start
    end = int((leavings[k] - quarter) * sequence_recover.STEP_RATE) - context.span.start
    window[max(begin, 0):max(min(end, word.size), 0)] = 1.0
    return window


SPECTRUM_TOP = 5000.0   # Hz: the harmonics compared (tube_recover.SpectrumTarget's top)


def hold_levels(flow: numpy.ndarray) -> Optional[numpy.ndarray]:
    """A hold's harmonic levels up to SPECTRUM_TOP, read at its own period."""
    if len(flow) == 0 or not numpy.all(numpy.isfinite(flow)) or numpy.ptp(flow) <= 0.0:
        return None
    try:
        period = posture_recover.period_of(flow)
    except IndexError:
        return None
    return numpy.abs(tube_recover.harmonics(flow, period, SPECTRUM_TOP))


def worker_score(task: Tuple[numpy.ndarray, Tuple[int, int]]) -> float:
    """A candidate's hold scored by its harmonic levels against the true
    hold's (rms of the log ratio, nepers): the coarse guide that ranks tubes
    by their distance from the truth far from it, where the equation error
    is flat (tube_recover.SpectrumTarget; MEASURED there, Spearman
    0.40-0.82; here the hold equation error sat at -25 to -28 dB for 257
    surveyed postures)."""
    vector, (begin, end) = task
    out = pyramid_word.render_command(_problem.context, _problem.command(vector))
    ours = hold_levels(out[begin:end])
    theirs = hold_levels(_problem.context.true_out[begin:end])
    if ours is None or theirs is None:
        return float("inf")
    size = min(len(ours), len(theirs))
    return float(numpy.sqrt(numpy.mean(numpy.log(numpy.maximum(ours[:size], 1e-12) /
                                                  numpy.maximum(theirs[:size], 1e-12)) ** 2)))


def survey_postures(pool, problem: "BodyProblem", vector: numpy.ndarray, arrivals: numpy.ndarray,
                    leavings: numpy.ndarray) -> numpy.ndarray:
    """Each vowel's posture by a survey of the body's range: SURVEY postures
    (a scrambled Sobol sequence over each heard parameter's range across the
    hint table's vowels; the held ones at the start) scored on that vowel's
    own hold, the best kept. (From the hints, the hold sound in body
    coordinates is a plateau, about -27 dB over the first 30% of the way to
    the truth: no slope for a local solver there.)"""
    context = problem.context
    word = context.word
    table = word.hints
    ranges = []
    for field in posture_recover.FIELDS:
        values = [table[name][field] for name in table]
        ranges.append((min(values), max(values)))
    offset = problem.rate_size()
    width = problem.posture_width()
    heard = [j for j, field in enumerate(posture_recover.FIELDS) if field not in FREEZE_FIELDS]
    sampler = scipy.stats.qmc.Sobol(d=len(heard), scramble=True, seed=1)
    points = sampler.random(SURVEY)
    vector = vector.copy()
    for index, k in enumerate(problem.vowels):
        window = vowel_window(context, arrivals, leavings, k)
        inside = numpy.flatnonzero(window) + context.span.start
        stretch = (int(inside[0]), int(inside[-1]) + 1)
        candidates = []
        for point in points:
            trial = vector.copy()
            for column, j in enumerate(heard):
                low, high = ranges[j]
                trial[offset + index * width + j] = low + point[column] * (high - low)
            candidates.append(trial)
        candidates.append(vector.copy())
        scores = pool.map(worker_score, [(candidate, stretch) for candidate in candidates])
        best = int(numpy.argmin(scores))
        print("  survey of %s: %d postures, best %.3f nepers rms (the start %.3f)" % (
            word.names[k], len(candidates), scores[best], scores[-1]), flush=True)
        vector[offset + index * width:offset + (index + 1) * width] = \
            candidates[best][offset + index * width:offset + (index + 1) * width]
    return vector


ROUNDS = int(os.environ.get("ROUNDS", "4"))   # read-and-solve rounds (main_iterate)


def release_knots(problem: TimeProblem, command: numpy.ndarray, read: dict) -> List[Tuple[float, float]]:
    """Knots at the releases read (a closure fades and reads late: Echo's
    velar closure read 39 steps late even through the true timing): each
    release step's tick, its value where the postures cross the seal from
    the stop toward the next phone."""
    context = problem.context
    word = context.word
    knots = []
    for k, (closure, release) in read.items():
        if release < 0:
            continue
        node = word.stops[k]
        values = fields.level_crossing(problem.rows, node, k, word.seal)
        value = values[1]
        if not numpy.isfinite(value):
            continue
        # The command tick whose engine release falls on the step read: the
        # glided area crosses the seal a few steps after the command does,
        # so the knot sits at the command's own crossing nearest the read.
        crossings = pyramid_word.crossings(command[:, node], word.seal)
        if not crossings:
            continue
        nearest_crossing = nearest(crossings, release / float(pyramid_word.TICK))
        shift = release / float(pyramid_word.TICK) - sealed_release_tick(problem, command, node, nearest_crossing)
        knots.append((nearest_crossing + shift, value))
    return sorted(knots)


def sealed_release_tick(problem: TimeProblem, command: numpy.ndarray, node: int, crossing: float) -> float:
    """The engine's release (in ticks) for the command's crossing at
    `crossing`: the sealed span's end after it."""
    first, last = pyramid_word.sealed_span(problem.context, command, node, int((crossing - 5.0) * pyramid_word.TICK))
    return (last if last >= 0 else crossing * pyramid_word.TICK) / float(pyramid_word.TICK)


def main_iterate() -> None:
    """Stage `iterate` (time fields, true postures): read the releases off
    the sound through the current timing, knot them, solve; again through the
    solved timing, until the reading holds. (A velar release read through
    the hinted timing was 3674 steps off, through a knot-free solve 223, and
    through the true timing exact: the front cavity's motion must be right
    for its release to read right.)"""
    specs = sys.argv[1].split(",")
    context = pyramid_word.Context(specs, ["areas"])
    word = context.word
    truth = fields.progress(word.true_arrivals, word.true_leavings, context.ticks, 0.001)
    start = fields.progress(word.hinted_arrivals, word.hinted_leavings, context.ticks, 0.001)
    problem = TimeProblem(context, min(LEVELS_MAX, context.levels))
    vector = problem.vector_of(start)
    true_read = pyramid_word.read_steps(context, context.true_command, word.hinted_arrivals, word.hinted_leavings)
    previous = None
    for round_index in range(ROUNDS):
        problem.knots = []
        command = problem.command(vector)
        read = pyramid_word.read_steps(context, command, word.hinted_arrivals, word.hinted_leavings)
        print("round %d: releases read %s (through the truth %s)" % (round_index, read, true_read), flush=True)
        if read == previous:
            break
        previous = read
        problem.knots = release_knots(problem, command, read) if round_index > 0 else []
        print("  knots %s" % ["(%.3f, %.4f)" % knot for knot in problem.knots], flush=True)
        with multiprocessing.Pool(POOL, initializer=worker_start, initargs=(specs, problem.knots)) as pool:
            report(problem, vector, pyramid_word.level_db(problem.residual(vector)), "round %d start" % round_index,
                   truth)
            vector = levenberg(pool, problem, vector, truth, None)
        numpy.save(os.path.join(pyramid_word.sequence_recover.SCRATCH,
                                "fields_iterate_%s.npy" % "_".join(specs).replace("!", "")), vector)


def adduction_shifted(times_folds: numpy.ndarray, arrivals: numpy.ndarray, leavings: numpy.ndarray, k: int,
                      shift: float, ticks: int) -> numpy.ndarray:
    """The folds' clock with the transition into phone k (leaving k - 1,
    arriving k) moved by `shift` seconds."""
    moved_arrivals = numpy.array(arrivals, dtype=numpy.float64)
    moved_leavings = numpy.array(leavings, dtype=numpy.float64)
    moved_leavings[k - 1] += shift
    moved_arrivals[k] += shift
    return fields.progress(moved_arrivals, moved_leavings, ticks, 0.001)


def main_onset() -> None:
    """Stage `onset` (time fields, true postures): each voicing onset after a
    voiceless phone read off the sound (the first collision in the glottal
    flow the sound implies through the start's tube) and the folds' clock
    aligned so the candidate's own first collision lands there (bisection
    on the adduction's timing), then held; the mouth clock solved, no mask.
    (Near the voicing threshold the folds' response turns on microseconds
    of phase: masks either left the settling cycles or hid the timing.)"""
    specs = sys.argv[1].split(",")
    context = pyramid_word.Context(specs, ["areas"])
    word = context.word
    truth = fields.progress(word.true_arrivals, word.true_leavings, context.ticks, 0.001)
    start = fields.progress(word.hinted_arrivals, word.hinted_leavings, context.ticks, 0.001)
    problem = TimeProblem(context, min(LEVELS_MAX, context.levels))
    vector = problem.vector_of(start)
    command = problem.command(vector)
    implied = pyramid_word.implied_flow(context, command)
    period = posture_recover.period_of(implied[context.span])
    arrivals = numpy.array(word.hinted_arrivals, dtype=numpy.float64)
    leavings = numpy.array(word.hinted_leavings, dtype=numpy.float64)
    folds = start.copy()
    for k in range(1, word.phones):
        if (k - 1) not in word.voiceless or k in word.voiceless:
            continue
        begin = int(leavings[k - 1] * sequence_recover.STEP_RATE)
        read = pyramid_word.first_collision(context, begin, period, implied)
        true_read = pyramid_word.first_collision(context, begin, period, context.true_in)

        def onset_at(shift: float) -> int:
            trial = adduction_shifted(folds, arrivals, leavings, k, shift, context.ticks)
            problem.fixed_folds = trial
            flow, _ = pyramid_word.render_both(context, problem.command(vector))
            return pyramid_word.first_collision(context, begin, period, flow)

        low, high = -0.02, 0.02
        for _ in range(30):
            middle = 0.5 * (low + high)
            if onset_at(middle) < read:
                low = middle
            else:
                high = middle
        folds = adduction_shifted(folds, arrivals, leavings, k, 0.5 * (low + high), context.ticks)
        problem.fixed_folds = folds
        print("onset into %s: read %.3f ms (through the truth %.3f); the folds' adduction moved %.3f ms" % (
            word.names[k], read / 96.0, true_read / 96.0, 1000 * 0.5 * (low + high)), flush=True)
    problem.fixed_folds = folds
    width = problem.reads().shape[1]
    problem.frozen = numpy.zeros(problem.size(), dtype=bool)
    problem.frozen[width:] = True
    with multiprocessing.Pool(POOL, initializer=worker_start,
                              initargs=(specs, [], "time", problem.used, None, None, folds)) as pool:
        report(problem, vector, pyramid_word.level_db(problem.residual(vector)), "start", truth)
        vector = levenberg(pool, problem, vector, truth, None)
    numpy.save(os.path.join(pyramid_word.sequence_recover.SCRATCH,
                            "fields_onset_%s.npy" % "_".join(specs).replace("!", "")), vector)


def main_staged() -> None:
    """Stage `staged`: the postures from the holds' central halves with the
    clocks held; then the clocks with the postures held; then all together.
    Each stage solves what its stretch of the sound determines. (Solved all
    at once, posture errors leaked into the timing: the first step put the
    20% basin test's timing 0.020 -> 0.107 phones off.)"""
    specs = sys.argv[1].split(",")
    context = pyramid_word.Context(specs, ["areas"])
    word = context.word
    truth = fields.progress(word.true_arrivals, word.true_leavings, context.ticks, 0.001)
    start = fields.progress(word.hinted_arrivals, word.hinted_leavings, context.ticks, 0.001)
    problem = BodyProblem(context, min(LEVELS_MAX, context.levels), word.hints)
    tubes = problem.base_tubes
    hinted = problem.command(problem.vector_of_start(start))
    read = pyramid_word.read_steps(context, hinted, word.hinted_arrivals, word.hinted_leavings)
    arrivals, leavings = pyramid_word.align_events(context, tubes, word.hinted_arrivals, word.hinted_leavings,
                                                   word.true_velum, read)
    start = fields.progress(arrivals, leavings, context.ticks, 0.001)
    vector = problem.vector_of_start(start)
    # The unclipped command: the clock is linear between targets, so its
    # crossing between ticks is where the unclipped area crosses (the
    # clipped area's interpolation put the knots off: the truth's hold sound
    # at -59.2 dB instead of the floor).
    aligned = problem.command(vector)
    for node in seal_nodes(context):
        for crossing in pyramid_word.crossings(aligned[:, node], word.seal):
            lower = int(numpy.floor(numpy.interp(crossing, numpy.arange(context.ticks), start)))
            problem.event_knots.append((crossing, node, lower))
    problem.event_knots.sort()
    rates = numpy.zeros(problem.size(), dtype=bool)
    rates[:problem.rate_size()] = True
    holds = hold_window(context, arrivals, leavings)
    stages = (("postures from the holds", rates, holds),
              ("clocks with the postures held", ~rates, None),
              ("all together", numpy.zeros(problem.size(), dtype=bool), None))
    print("%s: staged; events read (steps) %s; hold window %d of %d steps" % (
        " ".join(specs), read, int(holds.sum()), len(holds)), flush=True)
    saved = os.path.join(pyramid_word.sequence_recover.SCRATCH,
                         "fields_staged_%s.npy" % "_".join(specs).replace("!", ""))
    held = unheard(problem)
    if FREEZE_FIELDS:
        print("held at the start (unheard): %s" % ",".join(FREEZE_FIELDS), flush=True)
    if SURVEY > 0:
        print("stage: survey of the postures", flush=True)
        with multiprocessing.Pool(POOL, initializer=worker_start,
                                  initargs=(specs, problem.event_knots, "staged", problem.used, None)) as pool:
            vector = survey_postures(pool, problem, vector, arrivals, leavings)
        problem.window = None
        report(problem, vector, pyramid_word.level_db(problem.residual(vector)), "after the survey", truth)
    for label, frozen, window in stages:
        problem.frozen = frozen | held
        problem.window = window
        print("stage: %s" % label, flush=True)
        with multiprocessing.Pool(POOL, initializer=worker_start,
                                  initargs=(specs, problem.event_knots, "staged", problem.used, window)) as pool:
            report(problem, vector, pyramid_word.level_db(problem.residual(vector)), "stage start", truth)
            vector = levenberg(pool, problem, vector, truth, read)
        numpy.save(saved, vector)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[2] == "staged":
        main_staged()
        sys.exit(0)
    if len(sys.argv) > 2 and sys.argv[2] == "iterate":
        main_iterate()
        sys.exit(0)
    if len(sys.argv) > 2 and sys.argv[2] == "onset":
        main_onset()
        sys.exit(0)
    if len(sys.argv) > 2 and sys.argv[2] in ("full", "body"):
        main_full()
    else:
        main()
