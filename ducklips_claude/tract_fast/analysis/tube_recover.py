"""The tube from sound: the area function the sound implies, before any
posture (posture_recover's search over the organs found posture impostors
at -30 to -35 dB mismatch while the true posture scores -139 dB; twenty
starts of a local search over the nine organs never reached it).

The tube is searched directly. Its residual is the glottal flow our folds
make driving a candidate tube (sustain, the reversible engine) less the
flow the sound implies through that tube (reverse_core.peel_tube); the
true tube makes them the same flow. Levenberg-Marquardt on that residual,
column by column in parallel, coarse to fine: first only the residual's
low band (the first formants' phase only wraps above it), then wider.

Known-answer test: a posture's tube is rendered forward; its sound alone is
carried back to a tube, which is rendered forward again.
usage: tube_recover.py TRUE_POSTURE HINT_POSTURE   (key=value lists)"""
import os
import sys
from multiprocessing import Pool

from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy
from scipy.optimize import least_squares, leastsq
from scipy.signal import butter, sosfiltfilt
from scipy.stats import qmc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import posture_recover
import reverse_core
import driver_ledger
import column_locality

# The driver's shortcuts (NOTES_driver.md). EXACT: changes that leave every
# number the same (memoized residuals and Jacobians, columns of parameters
# that do not enter the sound, stages that would repeat one that did not
# move, renders stopped at the compared stretch's end). (EXACT, below, is
# the node search's depth in dB.) ALGO: named changes
# of the search's path, judged by the word's final round trip.
EXACT_SHORTCUTS = os.environ.get("FAST_EXACT", "1") == "1"
ALGO = set(name for name in os.environ.get("FAST_ALGO", "").split(",") if name)

NODE = posture_recover.NODE
STEP_RATE = reverse_core.STEP_RATE
CONTROLS = 16      # log-area control points along the tract (first stage)
DEFAULT_SEAL = 0.0 # the seal area targets take (a word's driver sets it)
DEFAULT_JET = 0.0  # the reversible jet's gain targets take (sustain REVJET; a word's driver sets it)
DEFAULT_VOICELESS = ()  # the voiceless phones' indices targets take (a word's driver sets it)
DEFAULT_FRICATIVES = ()  # the fricatives' indices targets take (a word's driver sets it)
BEND = 0.1         # weight of the profile's node-to-node bend on the spectrum's rungs (per neper)
STEP = 1e-3        # finite-difference step (log area, and length in cm)


class Target:
    """The sound being carried back and the stretch compared."""

    def __init__(self, flow_out, start, length, band, scale=1.0, equation=False, bend=0.0):
        self.flow_out = flow_out
        self.equation = equation
        self.bend = bend
        self.seal = DEFAULT_SEAL   # a node this narrow seals (sustain SEAL); closures then run forward
        self.jet = DEFAULT_JET     # the reversible jet's gain (sustain REVJET)
        self.voiceless = DEFAULT_VOICELESS   # phones whose folds are held apart (sustain FOLDSEQ)
        self.hidden_only = False   # a word's residual inside the folds' windows only (a release's time)
        self.fricatives = DEFAULT_FRICATIVES   # phones whose jet is on (sustain JETSEQ)
        self.output_from = None    # a step: our whole sound against the sound, our folds free from there
        self.context = None   # a sequence around the tube (sequence_recover.PrefixContext)
        self.start = start
        self.length = length
        self.band = band
        self.scale = scale
        self.seconds = max(0.3, (start + length) / STEP_RATE + 0.01)
        self.filter = butter(4, band / (0.5 * STEP_RATE), output="sos") if band > 0.0 else None


def harmonics(flow, period, top):
    """A stretch's harmonics (complex, Hann-windowed over whole cycles) up to `top` Hz."""
    cycles = len(flow) // period
    x = flow[:cycles * period] - numpy.mean(flow[:cycles * period])
    spectrum = numpy.fft.rfft(x * numpy.hanning(len(x)))
    count = int(top / (STEP_RATE / period))
    return spectrum[numpy.arange(1, count + 1) * cycles]


class SpectrumTarget:
    """The sound's harmonics over the stretch compared, the smooth form of
    the same test. With the peel the forward pass's exact inverse, our folds'
    flow equals the flow the sound implies exactly when our sound equals the
    sound; per harmonic, the log of our sound's amplitude over the sound's
    is that test, its real part the levels' difference (nepers), its
    imaginary part the phases' (radians, wrapped). The levels alone (no
    phase) are the coarse guide: over 512 surveyed postures their rms
    difference ranked the tubes by their log-area shape's distance from the
    truth with Spearman 0.70 / 0.82 / 0.40 (bdl5 /ae/ /ah/ /aa/), where the
    waveform test, exact at the truth, ranked them by nothing."""

    def __init__(self, flow_out, start, length, top=5000.0, phase=False):
        self.flow_out = flow_out
        self.start = start
        self.length = length
        self.top = top
        self.phase = phase
        self.context = None
        self.seconds = max(0.3, (start + length) / STEP_RATE + 0.01)
        self.period = posture_recover.period_of(flow_out[start:start + length])
        self.harmonics = harmonics(flow_out[start:start + length], self.period, top)


def residual_of_spectrum(tube, target):
    """Each sound read at its own period: our folds settle at a pitch the
    tube's load moves (period 761 against 769 steps 0.05 spans from the true
    /ae/), and read at the sound's period our upper harmonics were misread.
    Phases relative to the fundamental's (theta_k - k theta_1), so a shift
    in time is no difference."""
    if target.context is not None:
        sound = target.context.sound(tube, target.seconds)
    else:
        _, sound = posture_recover.render(tube, target.seconds)
    clip = sound[target.start:target.start + target.length]
    if len(clip) < target.length or not numpy.all(numpy.isfinite(clip)) or numpy.ptp(clip) <= 0.0:
        # No voice to read there (a render stopped short or silent): as far
        # as can be. (Oscar's /er/ survey crashed finding no period.)
        return numpy.full(len(target.harmonics) * (2 if target.phase else 1), 1e3)
    try:
        period = posture_recover.period_of(clip)
    except IndexError:
        return numpy.full(len(target.harmonics) * (2 if target.phase else 1), 1e3)
    ours = harmonics(clip, period, target.top)
    size = min(len(ours), len(target.harmonics))
    ratio = ours[:size] / target.harmonics[:size]
    levels = numpy.log(numpy.maximum(numpy.abs(ratio), 1e-12))
    if not target.phase:
        return pad(levels, len(target.harmonics))
    order = numpy.arange(1, size + 1)
    turn = numpy.angle(ratio) - order * numpy.angle(ratio[0])
    return numpy.concatenate([pad(levels, len(target.harmonics)), pad(numpy.angle(numpy.exp(1j * turn)), len(target.harmonics))])


def pad(values, size):
    """A residual of fixed length (a candidate with fewer harmonics below the
    top leaves zeros)."""
    out = numpy.zeros(size)
    out[:min(size, len(values))] = values[:size]
    return out


def residual_for(tube, target):
    if isinstance(target, SpectrumTarget):
        residual = residual_of_spectrum(tube, target)
    elif target.equation and target.context is not None:
        residual = target.context.equation(tube, target)
    elif target.equation:
        residual = residual_of_equation(tube, target)
    else:
        residual = residual_of_tube(tube, target)
    # Some postures' tubes do not render (the engine stops on a non-finite
    # sample): a wall the search turns back from.
    return numpy.where(numpy.isfinite(residual), residual, 1e3)


def fixed_clip(difference, target):
    """The compared stretch, mean removed, over the target's scale, always
    of the stretch's length (a render the engine stopped on a non-finite
    sample comes back short: the missing steps count as large; so do
    non-finite steps, a peel through a wrong closing that overflowed: as
    NaN, the timing grid's argmin took them, /aa b aa/'s closure 23 ms
    early)."""
    clip = difference[target.start:target.start + target.length]
    finite = numpy.isfinite(clip)
    if not numpy.all(finite):
        kept = clip[finite]
        mean = kept.mean() if len(kept) > 0 else 0.0
        clip = numpy.where(finite, (clip - mean) / target.scale, 1e3)
    else:
        clip = (clip - clip.mean()) / target.scale
    if len(clip) < target.length:
        clip = numpy.concatenate([clip, numpy.full(target.length - len(clip), 1e3)])
    return clip


def residual_of_equation(tube, target):
    """The folds driven by the tract the sound implies through this tube:
    the flow they pass less the implied flow (the equation error). The
    peel gives the first node's pressure and flow, so the wave the tube
    returns toward the glottis, (p - U rho c / A0) / 2; our folds read it
    (sustain FORCEWAVE) in place of their own tube's. Each step's error is
    then local: run free, our folds settle at a pitch the tube's load moves,
    and a nearly right tube's flow drifts in phase over the 15 periods."""
    areas, tract_length = tube
    implied, pressure = reverse_core.peel_tube_glottis(areas, target.flow_out, reverse_core.lining_constants(areas),
                                                       viscous=posture_recover.VISCOUS, hop=posture_recover.HOP,
                                                       end_extra=tract_length - len(areas) * NODE)
    returning = 0.5 * (pressure - implied * reverse_core.RHO_C / areas[0])
    forced, _ = posture_recover.render(tube, target.seconds, forced=returning)
    size = min(len(forced), len(implied))
    difference = forced[:size] - implied[:size]
    if target.filter is not None:
        difference = sosfiltfilt(target.filter, difference)
    return fixed_clip(difference, target)


def equation_mismatch(tube, target):
    """The equation error's level against the implied flow's swing (dB)."""
    driver_ledger.count_main(1)
    if target.context is not None:
        implied = target.context.implied(tube, target)
        error = target.context.equation(tube, target) * target.scale
    else:
        areas, tract_length = tube
        implied = reverse_core.peel_tube(areas, target.flow_out, reverse_core.lining_constants(areas),
                                         viscous=posture_recover.VISCOUS, hop=posture_recover.HOP,
                                         end_extra=tract_length - len(areas) * NODE)
        error = residual_of_equation(tube, target) * target.scale
    clip = implied[target.start:target.start + target.length]
    return 10.0 * numpy.log10(numpy.sum(error ** 2) / max(numpy.sum((clip - clip.mean()) ** 2), 1e-30))


def residual_of_tube(tube, target):
    """Driven less implied glottal flow over the compared stretch, low-passed
    to the stage's band, over the target's scale."""
    areas, tract_length = tube
    implied = reverse_core.peel_tube(areas, target.flow_out, reverse_core.lining_constants(areas),
                                     viscous=posture_recover.VISCOUS, hop=posture_recover.HOP,
                                     end_extra=tract_length - len(areas) * NODE)
    driven, _ = posture_recover.render(tube, target.seconds)
    size = min(len(driven), len(implied))
    difference = driven[:size] - implied[:size]
    if target.filter is not None:
        difference = sosfiltfilt(target.filter, difference)
    return fixed_clip(difference, target)


# Stage one: a smooth profile (log areas at CONTROLS points along a tract of
# a fixed count of whole nodes, and the fraction of a node past them). Stage
# two: every node's log area and that fraction. The count is tried a few
# ways: across a change of count the residual jumps, and differences taken
# across the jump misled the first search (from the hint it stopped at
# -27.9 dB with 47 nodes where the truth has 46).

class Shape:
    """A parameter vector's tube: `count` whole nodes, their log areas at
    `controls` evenly spaced points (interpolated) or, with no controls, one
    per node; the last parameter the fraction of a node past them."""

    def __init__(self, count, controls):
        self.count = count
        self.controls = controls

    def tube(self, parameters):
        fraction = float(numpy.clip(parameters[-1], 0.0, 2.45))
        if self.controls:
            centres = (numpy.arange(self.count) + 0.5) / self.count
            grid = numpy.linspace(0.0, 1.0, self.controls)
            log_areas = numpy.interp(centres, grid, parameters[:-1])
        else:
            log_areas = parameters[:-1]
        return numpy.exp(numpy.clip(log_areas, numpy.log(0.02), numpy.log(30.0))), (self.count + fraction) * NODE

    def residual(self, parameters, target):
        residual = residual_for(self.tube(parameters), target)
        bend = getattr(target, "bend", 0.0)
        if bend > 0.0 and not self.controls:
            residual = numpy.concatenate([residual, bend * numpy.diff(parameters[:-1], 2)])
        if isinstance(target, SpectrumTarget) and not self.controls:
            # Fewer harmonics than nodes: the spectrum's rungs also ask the
            # profile to bend little from node to node (the exact rungs do not).
            residual = numpy.concatenate([residual, BEND * numpy.diff(parameters[:-1], 2)])
        return residual

    def start_from(self, tube):
        """These parameters for the nearest tube of this shape to `tube`
        (its log-area profile stretched over this count)."""
        areas, tract_length = tube
        old = (numpy.arange(len(areas)) + 0.5) / len(areas)
        if self.controls:
            grid = numpy.linspace(0.0, 1.0, self.controls)
        else:
            grid = (numpy.arange(self.count) + 0.5) / self.count
        # Strictly inside the fraction's bounds: started on one, its
        # difference column read zero and the length never moved (bdl5 /aa/,
        # every posture basin at 2.5 nodes of air mass; no progress at all).
        fraction = float(numpy.clip(tract_length / NODE - self.count, 0.01, 2.4))
        return numpy.append(numpy.interp(grid, old, numpy.log(areas)), fraction)


class Problem:
    """One stage's residual and its columns (in parallel).

    Exact shortcuts (EXACT_SHORTCUTS): residuals and Jacobians memoized by the
    parameters' bytes (scipy's least_squares evaluated x0 and its Jacobian
    before MINPACK asked for them again, and every Jacobian's base point is
    the point MINPACK has just evaluated); a shape's inert columns
    (parameters that do not enter the sound: a stop's own log areas, which
    SequenceShape.split replaces by its vowel's) built from the base
    residual without a render."""

    MEMO = 6   # residuals kept (MINPACK revisits only the current point)

    def __init__(self, shape: Any, pool: Any, target: Any) -> None:
        self.shape = shape
        self.pool = pool
        self.target = target
        self.evaluations = 0
        self.residuals: Dict[bytes, numpy.ndarray] = {}
        self.residual_order: List[bytes] = []
        self.jacobians: Dict[bytes, numpy.ndarray] = {}
        self.first_residual: Optional[numpy.ndarray] = None   # the start's (MINPACK's first call)
        self.locality = column_locality.ColumnLocality(shape, target)

    def remember(self, key: bytes, value: numpy.ndarray) -> None:
        if key not in self.residuals:
            self.residual_order.append(key)
        self.residuals[key] = value
        while len(self.residual_order) > self.MEMO:
            self.residuals.pop(self.residual_order.pop(0), None)

    def residual(self, parameters: numpy.ndarray) -> numpy.ndarray:
        point = numpy.array(parameters, float)
        key = point.tobytes()
        if EXACT_SHORTCUTS and key in self.residuals:
            driver_ledger.count_saved(1)
            return self.residuals[key].copy()
        self.evaluations += 1
        driver_ledger.count_main(1)
        value = numpy.asarray(self.shape.residual(point, self.target), float)
        self.remember(key, value)
        if self.first_residual is None:
            self.first_residual = value
        return value.copy()

    def resumed_columns(self, point: numpy.ndarray,
                        moved_points: List[numpy.ndarray]) -> Tuple[List[numpy.ndarray], numpy.ndarray]:
        """The columns in batches, one per worker (RESUME_BATCHES): each
        batch computes the base with its record once, then its columns
        resumed. Columns dealt round robin in order of their first input
        divergence, so every batch holds early and late ones alike. Each
        batch's base counts as an evaluation (one more than a whole
        Jacobian's per batch but the first)."""
        batches = RESUME_BATCHES if RESUME_BATCHES > 0 else driver_ledger.POOL
        locality = column_locality.ColumnLocality(self.shape, self.target)
        steps = locality.resume_steps(point, moved_points)
        order = sorted(range(len(moved_points)), key=steps.__getitem__)
        shares: List[List[int]] = [order[index::batches] for index in range(batches)]
        shares = [share for share in shares if share]
        tasks = [(self.shape, point.copy(), [moved_points[position] for position in share], self.target)
                 for share in shares]
        results = driver_ledger.pool_map(self.pool, column_batch, tasks, len(moved_points) + len(shares))
        self.evaluations += len(shares)
        values: List[Optional[numpy.ndarray]] = [None] * len(moved_points)
        for share, (base, columns) in zip(shares, results):
            for position, value in zip(share, columns):
                values[position] = value
        return values, numpy.asarray(results[0][0], float)

    def jacobian(self, parameters: numpy.ndarray) -> numpy.ndarray:
        point = numpy.array(parameters, float)
        key = point.tobytes()
        if EXACT_SHORTCUTS and key in self.jacobians:
            driver_ledger.count_saved(len(point) + 1)
            return self.jacobians[key].copy()
        # Each parameter's own step where the shape gives them (a word's
        # times: STEP, 1 ms, was wider than a release's whole basin, and the
        # joint search never moved).
        steps = self.shape.steps(len(point)) if hasattr(self.shape, "steps") else numpy.full(len(point), STEP)
        inert: Set[int] = set()
        if EXACT_SHORTCUTS and hasattr(self.shape, "inert_columns"):
            inert = set(self.shape.inert_columns())
        base: Optional[numpy.ndarray] = self.residuals.get(key) if EXACT_SHORTCUTS else None
        tasks: List[Tuple[Any, numpy.ndarray, Any]] = []
        moved_points: List[numpy.ndarray] = []
        for column in range(len(point)):
            moved = point.copy()
            moved[column] += steps[column]
            moved_points.append(moved)
            if column not in inert:
                tasks.append((self.shape, moved, self.target))
        if base is None:
            tasks.append((self.shape, point.copy(), self.target))
        self.locality.observe(point, moved_points)
        if (EXACT_SHORTCUTS and hasattr(self.shape, "resumable") and self.shape.resumable(self.target)
                and len(tasks) > driver_ledger.POOL):
            values, computed_base = self.resumed_columns(point, [task[1] for task in tasks[:len(moved_points) - len(inert)]])
            self.evaluations += len(moved_points) - len(inert)
            if base is None:
                base = computed_base
                self.remember(key, base)
            else:
                driver_ledger.count_saved(1)
        else:
            values = driver_ledger.pool_map(self.pool, column_residual, tasks)
            self.evaluations += len(tasks)
            if base is None:
                base = numpy.asarray(values.pop(), float)
                self.remember(key, base)
            else:
                driver_ledger.count_saved(1)
        driver_ledger.count_saved(len(inert))
        columns: List[numpy.ndarray] = []
        computed = iter(values)
        for column in range(len(point)):
            if column in inert:
                value = self.shape.inert_residual(moved_points[column], self.target, base)
            else:
                value = numpy.asarray(next(computed), float)
            columns.append((value - base) / steps[column])
        jacobian = numpy.stack(columns, axis=1)
        if EXACT_SHORTCUTS:
            self.jacobians = {key: jacobian}
        return jacobian.copy()


RESUME_BATCHES = int(os.environ.get("RESUME_BATCHES", "0"))   # 0: one batch per worker


def column_batch(task: Tuple[Any, numpy.ndarray, List[numpy.ndarray], Any]) -> Tuple[numpy.ndarray, List[numpy.ndarray]]:
    """A worker's share of a Jacobian: the base point once, with its record,
    then its columns resumed from it (SequenceShape.resumed_residuals)."""
    shape, base_point, moved_points, target = task
    return shape.resumed_residuals(base_point, moved_points, target)


def column_residual(task: Tuple[Any, numpy.ndarray, Any]) -> numpy.ndarray:
    shape, parameters, target = task
    return shape.residual(parameters, target)


def mismatch(tube, flow_out, start, length):
    """The whole-band consistency (dB) as posture_recover measures it."""
    return posture_recover.consistency(tube, flow_out, start, length)


def stage(pool: Any, shape: Any, parameters: numpy.ndarray, target: Any, budget: int) -> Tuple[numpy.ndarray, int]:
    """Levenberg-Marquardt (unbounded: the parameters are log areas and the
    shape clips them; the bounded trust-region search stopped where it
    started). With EXACT_SHORTCUTS, MINPACK's lmder is called through
    leastsq: the same routine with the same arguments as
    least_squares(method='lm') (factor 100, no diag), without the Jacobian
    least_squares computes at the result for its report (a whole Jacobian
    per stage, discarded here: the stalled joints' 314 evaluations were 145
    at x0, 24 trials and 145 for that report)."""
    found, used, gain = stage_gain(pool, shape, parameters, target, budget)
    return found, used


def stage_gain(pool: Any, shape: Any, parameters: numpy.ndarray, target: Any,
               budget: int) -> Tuple[numpy.ndarray, int, float]:
    """stage, and what it gained on its own residual (dB, from the start's
    sum of squares to the result's; no extra evaluation)."""
    problem = Problem(shape, pool, target)
    start = numpy.array(parameters, float)
    if EXACT_SHORTCUTS:
        found = leastsq(problem.residual, start, Dfun=problem.jacobian, full_output=True, col_deriv=False,
                        ftol=1e-12, xtol=1e-12, gtol=1e-14, maxfev=budget, factor=100.0, diag=None)
        result = numpy.asarray(found[0], float)
        final = numpy.asarray(found[2]["fvec"], float)
    else:
        solved = least_squares(problem.residual, start, jac=problem.jacobian, method="lm", max_nfev=budget,
                               xtol=1e-12, ftol=1e-12, gtol=1e-14)
        result = solved.x
        final = solved.fun
    initial = problem.first_residual
    gain = 0.0
    if initial is not None:
        gain = 10.0 * numpy.log10(max(float(numpy.sum(initial ** 2)), 1e-300) / max(float(numpy.sum(final ** 2)), 1e-300))
    return result, problem.evaluations, gain


def tube_rms(found, true):
    size = min(len(found[0]), len(true[0]))
    return float(numpy.sqrt(numpy.mean((found[0][:size] - true[0][:size]) ** 2)))


class PostureShape:
    """A point of the posture search (each organ -1..1 of its span about the
    hint) as its tube."""

    def __init__(self, hint, widen=1.0, count=None):
        self.hint = hint
        self.widen = widen
        self.count = count

    def tube(self, parameters):
        posture = dict(self.hint)
        for field, unit in zip(posture_recover.FIELDS, parameters):
            reach = self.widen * posture_recover.SPAN[field] * float(numpy.clip(unit, -1.0, 1.0))
            posture[field] = max(self.hint[field] + reach, posture_recover.FLOOR.get(field, -1e9))
        return posture_recover.nodes(posture, self.count)

    def residual(self, parameters, target):
        return residual_for(self.tube(parameters), target)


def survey_cost(task):
    """rms dB of a posture's harmonic levels against the sound's."""
    shape, parameters, target = task
    return 8.686 * float(numpy.sqrt(numpy.mean(shape.residual(parameters, target) ** 2)))


SURVEY = 1024       # hinted postures sampled (Sobol) for the global stage
WIDEN = 1.5         # the survey's reach about the hint, in posture_recover.SPAN
                    # (bdl5's vowels sit up to 1.71 spans from organs_vtl_tongue's)
COUNTS = 2          # node counts tried either side of each basin's
KEEP = 6            # of them, carried into the posture search
FINISH = 6          # of those, carried into the tube search (/er/'s right-count
                    # basin, 0.17 cm² off, ranked 4th by the posture stage)
BENDS = (1.0, 0.1, 0.01)   # the node search's first runs: weight of the profile's
                    # node-to-node bend (log area), then none. Its zigzags are the
                    # equation error's near-null directions, bumpy at -80 dB: /uw/
                    # from 2% noise stopped at -79.6 dB unbent, -112.6 bent so.
EXACT = -150.0      # dB: a node search this deep has found the tube; the other basins are skipped
NODE_RUNGS = 10     # at most this many Levenberg-Marquardt runs of 30 over the nodes (two without 1 dB of gain end it)
STALLS = 1 if "stall1" in ALGO else 2   # unbent rungs without 1 dB of gain that end a node search
HALVING_KEEP = 2    # FAST_ALGO halving: basins carried past the node search's first rung


def recover(pool, hint, flow_out, start, length, truth=None, count=None, context=None):
    """The tube the sound implies: a survey of postures about the hint, the
    posture search from the best of them, then every node free. (From one
    posture basin of five the node search went from -33.0 to -119.8 dB and
    from 0.52 to 0.057 cm² rms off the truth; two others stopped at -45.7
    and -34.8 dB. So several basins are carried and the deepest kept.)"""
    shape = PostureShape(hint, WIDEN, count)
    whole = Target(flow_out, start, length, 0.0, 1.0, True)
    guide = SpectrumTarget(flow_out, start, length)
    ladder = (guide, SpectrumTarget(flow_out, start, length, 5000.0, True),
              SpectrumTarget(flow_out, start, length, 10000.0, True), whole, whole)
    # In a sequence, every candidate is heard after the phones before it
    # (the peel and our folds need the history: searched alone from the
    # sound's start, every phone after the first failed, -13 dB).
    for rung in ladder:
        rung.context = context
    points = qmc.Sobol(len(posture_recover.FIELDS), seed=7).random(SURVEY) * 2.0 - 1.0
    points[0] = 0.0
    with driver_ledger.stage("survey"):
        costs = numpy.array(driver_ledger.pool_map(pool, survey_cost, [(shape, point, guide) for point in points]))
    order = numpy.argsort(costs)[:KEEP]
    print("survey: %d postures, levels rms best %.2f dB, hint %.2f dB, kept %s" % (
        SURVEY, costs[order[0]], costs[0], " ".join("%.1f" % costs[i] for i in order)))
    sys.stdout.flush()
    carried = []
    for basin, index in enumerate(order):
        with driver_ledger.stage("posture %d" % basin):
            carried.append(posture_search(pool, hint, shape, points[index], costs[index], guide, ladder, whole, truth,
                                          count))
    carried.sort(key=first_of)
    return finish_basins(pool, carried, flow_out, start, length, truth, count, context)


def posture_search(pool: Any, hint: Dict[str, float], shape: Any, point: numpy.ndarray, cost: float, guide: Any,
                   ladder: Sequence[Any], whole: Any, truth: Any, count: Optional[int]) -> Tuple[float, Any]:
    """One surveyed posture's basin: its count chosen by the levels, then
    the spectrum's rungs and the equation error."""
    # The posture searched at a fixed count of nodes (the tube then moves
    # continuously with it), the count chosen by the levels: each nearby
    # count gets the levels' rung, the best the whole ladder.
    natural = len(shape.tube(point)[0])
    counts = range(natural - COUNTS, natural + COUNTS + 1) if count is None else [count]
    trials = []
    for trial_count in counts:
        fixed = PostureShape(hint, WIDEN, trial_count)
        parameters, used = stage(pool, fixed, point, guide, 15)
        if EXACT_SHORTCUTS and len(counts) == 1:
            # One count (a word's): its levels decide nothing; not evaluated.
            driver_ledger.count_saved(1)
            trials.append((0.0, trial_count, parameters))
            continue
        driver_ledger.count_main(1)
        trials.append((float(numpy.sum(fixed.residual(parameters, guide) ** 2)), trial_count, parameters))
    trials.sort(key=first_of)
    levels, chosen, parameters = trials[0]
    fixed = PostureShape(hint, WIDEN, chosen)
    for rung in ladder[1:3]:
        parameters, used = stage(pool, fixed, parameters, rung, 25)
    coarse = equation_mismatch(fixed.tube(parameters), whole)
    parameters, used = stage(pool, fixed, parameters, whole, 25)
    found = fixed.tube(parameters)
    value = equation_mismatch(found, whole)
    note = "" if truth is None else "; tube rms %.3f (true count %d)" % (tube_rms(found, truth), len(truth[0]))
    print("  posture search from %.2f dB levels: %d nodes (natural %d), by the spectrum %.1f dB, then %.1f dB%s" % (
        cost, chosen, natural, coarse, value, note))
    sys.stdout.flush()
    return value, found


def finish_basins(pool: Any, carried: List[Tuple[float, Any]], flow_out: numpy.ndarray, start: int, length: int,
                  truth: Any, count: Optional[int], context: Any) -> Tuple[float, Any]:
    """The node search from each carried basin (deepest first), the best
    kept."""
    best = None
    memo: Dict[Tuple[bytes, float], Tuple[numpy.ndarray, int]] = {}
    carried = carried[:FINISH]
    if "halving" in ALGO and len(carried) > HALVING_KEEP:
        # Every basin's first rung, then only the deepest HALVING_KEEP on
        # (their first rungs reused from the memo). In the reference's
        # /aa m aa/ second /aa/, the basin that won (-79.0 dB) was second
        # by its posture (-41.3 against -60.3) but first after one rung
        # (-78.7); the basins that ended at -39.6 to -50.7 dB were already
        # 12 to 44 dB behind after it.
        firsts: List[Tuple[float, int]] = []
        for basin, (value, found) in enumerate(carried):
            with driver_ledger.stage("node search %d first rung" % basin):
                result, finished = node_search(pool, found, flow_out, start, length, count, context, memo, 1)
            firsts.append((result, basin))
        firsts.sort(key=first_of)
        print("  halving: first rungs %s; carried on: basins %s" % (
            " ".join("%.1f" % value for value, basin in sorted(firsts, key=second_of)),
            " ".join("%d" % basin for value, basin in firsts[:HALVING_KEEP])))
        sys.stdout.flush()
        carried = [carried[basin] for value, basin in firsts[:HALVING_KEEP]]
    for basin, (value, found) in enumerate(carried):
        with driver_ledger.stage("node search %d" % basin):
            result, finished = node_search(pool, found, flow_out, start, length, count, context, memo)
        note = "" if truth is None else "; tube rms %.3f" % tube_rms(finished, truth)
        print("  node search from %.1f dB: %.1f dB%s" % (value, result, note))
        sys.stdout.flush()
        if best is None or result < best[0]:
            best = (result, finished)
        if best[0] < EXACT:
            break   # one basin has reached the truth; the rest cannot beat it by anything audible
    return best


def node_search(pool, found, flow_out, start, length, count=None, context=None,
                memo: Optional[Dict[Tuple[bytes, float], Tuple[numpy.ndarray, int]]] = None,
                rungs: int = NODE_RUNGS):
    """Every node free from `found`, on the equation error; returns its
    level (dB) and the best tube."""
    whole = Target(flow_out, start, length, 0.0, 1.0, True)
    whole.context = context
    # Every node free, on the exact test only, the best tube kept. (With
    # the spectrum's rungs first, a posture's tube at -63.3 dB came out
    # at -44.7: near the truth they pull it away.)
    # Free to choose its length (no count given), the search turns air
    # mass into nodes past 1.2 nodes of it (the last area continued: what
    # the air mass stands for) and gives a node back when the air mass
    # is pinned at nothing. (Rounded vowels' small lip openings trade
    # against the open end's air mass; /ow/'s basins had 45-50 nodes, the
    # truth 48, and the search could not change its count.)
    current = (Shape(len(found[0]), 0).tube(Shape(len(found[0]), 0).start_from(found)))
    current_value = equation_mismatch(current, whole)
    finished = current
    result = current_value
    stalled = 0
    rung = 0
    # A rung from the same start with the same bend repeats an earlier
    # one exactly (the search is deterministic): its result is reused.
    # (A rung that did not improve leaves `current` as it was, and the next
    # unbent rung starts where it did: /aa b aa/'s second /aa/ ran rungs
    # 3 and 4 from one start in several basins.)
    rungs_run: Dict[Tuple[bytes, float], Tuple[numpy.ndarray, int]] = memo if memo is not None else {}
    while rung < min(rungs, NODE_RUNGS) and stalled < STALLS:
        nodes = Shape(len(current[0]), 0)
        bend = BENDS[rung] if rung < len(BENDS) else 0.0
        target = Target(flow_out, start, length, 0.0, 1.0, True, bend)
        target.context = context
        begin = nodes.start_from(current)
        key = (begin.tobytes(), bend)
        with driver_ledger.stage("rung %d" % rung):
            if EXACT_SHORTCUTS and key in rungs_run:
                parameters, used = rungs_run[key]
                driver_ledger.count_saved(used)
            else:
                parameters, used = stage(pool, nodes, begin, target, 30)
                rungs_run[key] = (parameters, used)
        attempt = nodes.tube(parameters)
        attempt_value = equation_mismatch(attempt, whole)
        extra = attempt[1] / NODE - len(attempt[0])
        print("    node rung %d: %.1f dB, %d nodes and %.2f of air mass" % (rung, attempt_value, len(attempt[0]), extra))
        stalled = stalled + 1 if attempt_value > current_value - 1.0 and bend == 0.0 else 0
        if attempt_value < current_value:
            current = attempt
            current_value = attempt_value
        if current_value < result:
            finished = current
            result = current_value
        if count is None:
            extra = current[1] / NODE - len(current[0])
            if extra > 1.2:
                current = recount(current, len(current[0]) + int(extra - 0.2))
                current_value = equation_mismatch(current, whole)
                stalled = 0
            elif extra < 0.05:
                current = recount(current, len(current[0]) - 1)
                current_value = equation_mismatch(current, whole)
                stalled = 0
        rung += 1
    return result, finished


def recount(tube, count):
    """A tube over `count` whole nodes, its length kept: nodes past the
    old count take the last area; air mass makes up the rest."""
    areas, tract_length = tube
    if count > len(areas):
        areas = numpy.concatenate([areas, numpy.full(count - len(areas), areas[-1])])
    else:
        areas = areas[:count]
    return areas, max(tract_length, (count + 0.01) * NODE)


def first_of(pair):
    return pair[0]


def second_of(pair):
    return pair[1]


def main():
    true = posture_recover.parse(sys.argv[1])
    hint = posture_recover.parse(sys.argv[2])
    tube = posture_recover.nodes(true)
    _, flow_out = posture_recover.render(tube)
    start = int(0.1 * STEP_RATE)
    period = posture_recover.period_of(flow_out[start:])
    length = 15 * period
    hint_tube = posture_recover.nodes(hint)
    print("true tube %.2f cm, %d nodes; hint %.2f cm, %d nodes; mismatch true %.1f hint %.1f dB; tube rms %.3f" % (
        tube[1], len(tube[0]), hint_tube[1], len(hint_tube[0]), mismatch(tube, flow_out, start, length),
        mismatch(hint_tube, flow_out, start, length), tube_rms(hint_tube, tube)))
    sys.stdout.flush()
    with Pool(10) as pool:
        value, found = recover(pool, hint, flow_out, start, length, tube)
    _, again = posture_recover.render(found)
    print("recovered: %.1f dB; tube rms %.3f; round trip: the recovered tube's sound against the true one's: %.1f dB" % (
        value, tube_rms(found, tube), reverse_core.snr(flow_out[start:start + length], again[start:start + length])))
    print("true  " + " ".join("%.2f" % a for a in tube[0]))
    print("found " + " ".join("%.2f" % a for a in found[0]))


if __name__ == "__main__":
    main()
