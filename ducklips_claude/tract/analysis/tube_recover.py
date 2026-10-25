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

import numpy
from scipy.optimize import least_squares
from scipy.signal import butter, sosfiltfilt
from scipy.stats import qmc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import posture_recover
import reverse_core

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
    ours = harmonics(clip, posture_recover.period_of(clip), target.top)
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
    """One stage's residual and its columns (in parallel)."""

    def __init__(self, shape, pool, target):
        self.shape = shape
        self.pool = pool
        self.target = target
        self.evaluations = 0

    def residual(self, parameters):
        self.evaluations += 1
        return self.shape.residual(parameters, self.target)

    def jacobian(self, parameters):
        # Each parameter's own step where the shape gives them (a word's
        # times: STEP, 1 ms, was wider than a release's whole basin, and the
        # joint search never moved).
        steps = self.shape.steps(len(parameters)) if hasattr(self.shape, "steps") else numpy.full(len(parameters), STEP)
        tasks = []
        for column in range(len(parameters)):
            moved = numpy.array(parameters, float)
            moved[column] += steps[column]
            tasks.append((self.shape, moved, self.target))
        tasks.append((self.shape, numpy.array(parameters, float), self.target))
        columns = self.pool.map(column_residual, tasks)
        self.evaluations += len(tasks)
        base = columns[-1]
        return numpy.stack([(c - base) / step for c, step in zip(columns[:-1], steps)], axis=1)


def column_residual(task):
    shape, parameters, target = task
    return shape.residual(parameters, target)


def mismatch(tube, flow_out, start, length):
    """The whole-band consistency (dB) as posture_recover measures it."""
    return posture_recover.consistency(tube, flow_out, start, length)


def stage(pool, shape, parameters, target, budget):
    """Levenberg-Marquardt (unbounded: the parameters are log areas and the
    shape clips them; the bounded trust-region search stopped where it
    started)."""
    problem = Problem(shape, pool, target)
    result = least_squares(problem.residual, parameters, jac=problem.jacobian, method="lm", max_nfev=budget,
                           xtol=1e-12, ftol=1e-12, gtol=1e-14)
    return result.x, problem.evaluations


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
    costs = numpy.array(pool.map(survey_cost, [(shape, point, guide) for point in points]))
    order = numpy.argsort(costs)[:KEEP]
    print("survey: %d postures, levels rms best %.2f dB, hint %.2f dB, kept %s" % (
        SURVEY, costs[order[0]], costs[0], " ".join("%.1f" % costs[i] for i in order)))
    sys.stdout.flush()
    carried = []
    for index in order:
        # The posture searched at a fixed count of nodes (the tube then moves
        # continuously with it), the count chosen by the levels: each nearby
        # count gets the levels' rung, the best the whole ladder.
        natural = len(shape.tube(points[index])[0])
        counts = range(natural - COUNTS, natural + COUNTS + 1) if count is None else [count]
        trials = []
        for trial_count in counts:
            fixed = PostureShape(hint, WIDEN, trial_count)
            parameters, used = stage(pool, fixed, points[index], guide, 15)
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
            costs[index], chosen, natural, coarse, value, note))
        sys.stdout.flush()
        carried.append((value, found))
    carried.sort(key=first_of)
    best = None
    for value, found in carried[:FINISH]:
        result, finished = node_search(pool, found, flow_out, start, length, count, context)
        note = "" if truth is None else "; tube rms %.3f" % tube_rms(finished, truth)
        print("  node search from %.1f dB: %.1f dB%s" % (value, result, note))
        sys.stdout.flush()
        if best is None or result < best[0]:
            best = (result, finished)
        if best[0] < EXACT:
            break   # one basin has reached the truth; the rest cannot beat it by anything audible
    return best


def node_search(pool, found, flow_out, start, length, count=None, context=None):
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
    while rung < NODE_RUNGS and stalled < 2:
        nodes = Shape(len(current[0]), 0)
        bend = BENDS[rung] if rung < len(BENDS) else 0.0
        target = Target(flow_out, start, length, 0.0, 1.0, True, bend)
        target.context = context
        parameters, used = stage(pool, nodes, nodes.start_from(current), target, 30)
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
