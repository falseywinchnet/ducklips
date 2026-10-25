"""A moving mouth from sound: a sequence of held tubes and the transitions
between them (sustain AREASEQ keyframes: each phone's tube held from its
arrival to its leaving, linear between).

The engine's frame motion depends only on the keyframes, so for a
hypothesis reverse_core.moving_geometry gives every step's geometry and
reverse_core.peel_moving_glottis carries the sound back through it. The
test is the equation error (tube_recover.residual_of_equation): our folds
driven by the tract the sound implies, against the flow it implies; the
true hypothesis makes them one flow.

Stages: each phone's tube from its hold (tube_recover.recover, hinted by the
phone's posture, at the first phone's node count: the moving tube keeps one
count and one length); each transition's leaving and arrival times on a
coarse-to-fine grid over a window around it; then Levenberg-Marquardt over
every tube and time.

Known-answer test: bdl5 profile postures for the phones, held and moved on
a plan; the sound alone carried back; the recovered sequence rendered again.
usage: sequence_recover.py PHONE,PHONE,...   (names in the bdl5 profile)"""
import hashlib
import os
import sys
from multiprocessing import Pool
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import posture_recover
import reverse_core
import tube_recover
import observer
import marching
import driver_ledger
import column_locality

STEP_RATE = reverse_core.STEP_RATE
DATA = os.path.join(posture_recover.ROOT, "data")
OUT = os.path.join(posture_recover.ROOT, "out", "reversible")
HOLD = 0.07             # s each phone is held in the known answer
TRANSITIONS = (0.05, 0.07, 0.06, 0.08)   # s each transition takes (cycled)
HINT_OFFSET = 0.0073    # s the hinted times are off the truth (alternating sign; off the timing
                        # grids' 0.2 ms lattice: at 8 ms the grids held the truth exactly)
SETTLE = 0.02           # s after the hinted arrival before a hold's window opens
CLEAR = 0.01            # s before the hinted leaving that it closes
GAP = 0.001             # s: keyframes never closer than this
TIME_STEP = 2e-5        # s: a time's finite-difference step (two engine steps; a release's basin is ~0.1-0.25 ms)
POLISH = 2              # rounds of a later phone's times, then its tube again
STOP_AREA = 1e-4        # a stop's sealed node (cm^2)
RELEASE_AREA = 0.05     # cm^2 its seal opens to, RELEASE_STEP after it leaves
RELEASE_STEP = 0.001    # s
# Then open to RELEASE_WIDE at RELEASE_OPEN after the leaving (or to the
# next phone's area there, if narrower), and on to the next phone: the
# opening's own motion, so the stretch the folds decide (narrower than
# HIDE_AREA) depends on the release's time alone. (Opened straight toward
# the next phone, its burst was shaped by that phone: with the next vowel at
# its hint the release came out 0.5 ms late, by the window's sound alone as
# by the whole.) 0: straight on.
RELEASE_OPEN = 0.004    # s
RELEASE_WIDE = 0.3      # cm^2
NASAL_HINT = 0.4        # a nasal's hinted velum port (cm^2)
JET_GAIN = 0.4          # the reversible jet's gain in a word with a fricative, gated on for the
                        # fricatives (jet_keyframes). (Its source is the flow of the step before, an
                        # explicit feedback: 0.8 ran away through 0.08 cm^2 with the folds apart;
                        # ungated, 0.3 ran away at /k!/'s release into /er/, and next to an open
                        # velum the jet in /ao/ and /er/ ran the march away at the truth.)
# A voiceless phone's folds (cord.cpp's postures: lower rest, upper rest
# (cm), chink (cm^2)): held apart from its arrival less ABDUCT_LEAD, brought
# together again from its leaving plus VOICE_ONSET, each move MOVE long.
MODAL_FOLDS = (0.00495, 0.0004, 0.0002)
PLOSIVE_FOLDS = (0.0501, 0.0501, 0.0002)
FRICATIVE_FOLDS = (0.10225, 0.10225, 0.0002)
ABDUCT_LEAD = 0.01      # s
VOICE_ONSET = 0.0       # s
FOLD_MOVE = 0.02        # s
NASAL_INVERSE = "march"  # the split march with folds-decided release windows ('window': the whole velum-open stretch forward)
SEAL = 0.01             # a node this narrow seals
# The forced folds' observer gain (sustain OBSERVE): pulled toward the flow
# the sound implies wherever it is revealed (0: none). Read from the
# environment so the pool's workers take it.
OBSERVE_GAIN = float(os.environ.get("OBSERVE_GAIN", "0"))
HIDE_AREA = 0.1         # cm^2: a node this narrow hides the inside as a seal does (0: the seal's).
                        # Hidden only once sealed, /aa d aa/'s closure 0.05 ms late cost -23.6 dB
                        # (truth -209.5) and its release was rugged; at 0.1 both rise monotonically
JOINT_BENDS = (1.0, 0.1, 0.0, 0.0)   # the prefix's joint runs (bend weights, then none)


def table(path):
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        name, rest = line.split(None, 1)
        out[name] = posture_recover.parse(rest.split("#")[0])
    return out


def plan(count):
    """Arrival and leaving times (s) of `count` phones."""
    arrivals = [0.0]
    leavings = [0.03 + HOLD]
    for index in range(1, count):
        arrivals.append(leavings[-1] + TRANSITIONS[(index - 1) % len(TRANSITIONS)])
        leavings.append(arrivals[-1] + HOLD)
    return numpy.array(arrivals), numpy.array(leavings)


def vowel_before(index, stops):
    """The vowel a consonant's tube is made from (its node then changed):
    the last phone before `index` that is not a consonant (a stop, nasal or
    fricative; the phone before's carried its own constriction on: /s k/),
    or for a word's leading consonants the first vowel after them."""
    base = index - 1
    while base in stops:
        base -= 1
    if base >= 0:
        return base
    base = index + 1
    while base in stops:
        base += 1
    return base


def jet_keyframes(arrivals, leavings, fricatives):
    """The reversible jet's gate keyframes (posture_recover.render's
    `jet_gate`): frication as a gesture, on through each fricative and
    ramped over the transitions into and out of it, off elsewhere (the
    vowels' own narrowings, /ao/'s 0.15 cm^2 pharynx, would otherwise
    hiss); None with no fricative."""
    if not [k for k in fricatives if k < len(arrivals)]:
        return None
    count = len(arrivals)
    shapes = []
    # (Only the phones there are: a prefix joint fits the first few.)
    for index in [k for k in sorted(fricatives) if k < count]:
        rise = leavings[index - 1] if index > 0 else arrivals[index]
        fall = arrivals[index + 1] if index + 1 < count else leavings[index] + 0.05
        shapes.append((rise, arrivals[index], leavings[index], fall))
    times = sorted(set([0.0] + [moment for shape in shapes for moment in shape]))
    gates = []
    for moment in times:
        value = 0.0
        for rise, full, held, fall in shapes:
            if full <= moment <= held:
                value = max(value, 1.0)
            elif rise < moment < full:
                value = max(value, (moment - rise) / (full - rise))
            elif held < moment < fall:
                value = max(value, (fall - moment) / (fall - held))
        gates.append(value)
    return times, gates


def jet_gate_for(arrivals, leavings, fricatives, steps):
    """jet_keyframes per step as the engine reads them (None with none)."""
    frames = jet_keyframes(arrivals, leavings, fricatives)
    if frames is None:
        return None
    return reverse_core.jet_gate_steps(frames[0], frames[1], steps)


def fold_keyframes(arrivals, leavings, voiceless):
    """The folds' posture keyframes (posture_recover.render's `folds`,
    sustain FOLDKEYS: linear between, so smooth in the phone times) for the
    voiceless phones ({index: posture}): from modal ABDUCT_LEAD before the
    arrival to the phone's posture FOLD_MOVE later, held, back to modal over
    FOLD_MOVE from the leaving plus VOICE_ONSET; where two overlap, the one
    further along. None with none. (As triggers at a step, a time moved by
    less than a step did nothing, then jumped: Echo's joint searches never
    moved, its time derivatives being spikes.)"""
    if not [k for k in voiceless if k < len(arrivals)]:
        return None
    shapes = []
    # (Only the phones there are: a prefix joint fits the first few; Papa's
    # second /p!/ crashed the first prefix.)
    for index in [k for k in sorted(voiceless) if k < len(arrivals)]:
        rise = arrivals[index] - ABDUCT_LEAD
        fall = leavings[index] + VOICE_ONSET
        shapes.append((rise, rise + FOLD_MOVE, fall, fall + FOLD_MOVE, numpy.asarray(voiceless[index], float)))
    times = sorted(set([0.0] + [max(moment, 0.0) for shape in shapes for moment in shape[:4]]))
    modal = numpy.asarray(MODAL_FOLDS, float)
    postures = []
    for moment in times:
        weight = 0.0
        target = modal
        for rise, full, held, fall, posture in shapes:
            if full <= moment <= held:
                share = 1.0
            elif rise < moment < full:
                share = (moment - rise) / (full - rise)
            elif held < moment < fall:
                share = (fall - moment) / (fall - held)
            else:
                share = 0.0
            if share > weight:
                weight = share
                target = posture
        postures.append(tuple(modal + weight * (target - modal)))
    return times, postures


def keyframes(tubes, arrivals, leavings):
    """The AREASEQ keyframes: each tube at its arrival and its leaving; a
    sealed phone (a node at STOP_AREA) opens its seal to RELEASE_AREA
    RELEASE_STEP after its leaving, then moves on to the next phone. (Opened
    only toward the next phone's area, the moment the seal broke depended on
    that area as much as on the leaving time: with the next vowel still at
    its hint, the release came out 7 ms late and the vowel after it was
    lost.)"""
    times = []
    frames = []
    count = len(tubes)
    for index in range(count):
        tube = numpy.asarray(tubes[index], float)
        times.extend([arrivals[index], leavings[index]])
        frames.extend([tube, tube])
        sealed = tube <= 2.0 * STOP_AREA
        if index + 1 < count and numpy.any(sealed):
            opened = tube.copy()
            opened[sealed] = RELEASE_AREA
            moment = min(leavings[index] + RELEASE_STEP, 0.5 * (leavings[index] + arrivals[index + 1]))
            times.append(moment)
            frames.append(opened)
            if RELEASE_OPEN > RELEASE_STEP:
                wider = tube.copy()
                following = numpy.asarray(tubes[index + 1], float)
                wider[sealed] = numpy.minimum(RELEASE_WIDE, following[sealed])
                later = min(leavings[index] + RELEASE_OPEN, 0.5 * (moment + arrivals[index + 1]))
                times.append(later)
                frames.append(wider)
    return numpy.array(times), numpy.array(frames)


class SequenceShape:
    """A hypothesis's parameters: every phone's node log areas (one count),
    then each transition's leaving and arrival times; the first arrival and
    the length fixed. Times kept in order, GAP apart."""

    def __init__(self, phones, count, length, stops=None, nasals=None, narrowings=None):
        self.phones = phones
        self.count = count
        self.length = length
        self.stops = stops if stops is not None else {}   # phone index -> seal node
        # Phone index -> its node's area where it is not a seal (a fricative).
        self.narrowings = narrowings if narrowings is not None else {}
        # Nasal phones (index order): each one's velum port area, a parameter
        # after the times.
        self.nasals = sorted(nasals) if nasals is not None else []
        self.time_count = 2 * phones - 1

    def split(self, parameters):
        tubes = numpy.exp(numpy.clip(parameters[:self.phones * self.count].reshape(self.phones, self.count),
                                     numpy.log(0.02), numpy.log(30.0)))
        for index, node in self.stops.items():
            tubes[index] = tubes[vowel_before(index, self.stops)]
            tubes[index, node] = self.narrowings.get(index, STOP_AREA)
        times = parameters[self.phones * self.count:self.phones * self.count + self.time_count]
        arrivals = numpy.zeros(self.phones)
        leavings = numpy.zeros(self.phones)
        for index in range(self.phones):
            if index > 0:
                arrivals[index] = max(times[2 * index - 1], leavings[index - 1] + GAP)
            leavings[index] = max(times[2 * index], arrivals[index] + GAP)
        return tubes, arrivals, leavings

    def join(self, tubes, arrivals, leavings, velum_areas=None):
        times = [leavings[0]]
        for index in range(1, self.phones):
            times.extend([arrivals[index], leavings[index]])
        ports = [] if velum_areas is None else [velum_areas[k] for k in self.nasals]
        return numpy.concatenate([numpy.log(numpy.asarray(tubes)).ravel(), times, ports])

    def steps(self, size):
        """The finite-difference step of each parameter: log areas and
        velum areas tube_recover.STEP, times TIME_STEP."""
        steps = numpy.full(size, tube_recover.STEP)
        steps[self.phones * self.count:self.phones * self.count + self.time_count] = TIME_STEP
        return steps

    def inert_columns(self) -> List[int]:
        """The parameters that do not enter the sound: a stop's (or a word's
        leading consonant's) own log areas, which split replaces by its
        vowel's. Their columns are the base residual (with a bend, its bend
        rows moved): /aa b aa/'s joint, 46 of 143 columns."""
        columns: List[int] = []
        for index in sorted(self.stops):
            columns.extend(range(index * self.count, (index + 1) * self.count))
        return columns

    def inert_residual(self, moved: numpy.ndarray, target: Any, base: numpy.ndarray) -> numpy.ndarray:
        """An inert column's residual from the base point's: the sound's
        part the same bytes (split gives the same tubes and times), the bend
        rows (as residual's) from the moved parameters."""
        bend = getattr(target, "bend", 0.0)
        if bend <= 0.0:
            return base.copy()
        logs = moved[:self.phones * self.count].reshape(self.phones, self.count)
        rows = self.phones * (self.count - 2)
        return numpy.concatenate([base[:len(base) - rows], bend * numpy.diff(logs, 2, axis=1).ravel()])

    def input_tracks(self, parameters: numpy.ndarray, target: Any) -> Dict[str, Any]:
        """The keyframe tracks the engine and the peels read for these
        parameters (column_locality.first_divergence)."""
        tubes, arrivals, leavings = self.split(parameters)
        times, frames = keyframes(tubes, arrivals, leavings)
        tracks: Dict[str, Any] = {"areas": (times, frames)}
        if self.nasals:
            velum = velum_keyframes(arrivals, leavings, self.velum(parameters))
            tracks["velum"] = (numpy.asarray(velum[0], float), numpy.asarray(velum[1], float))
        folds = fold_keyframes(arrivals, leavings, getattr(target, "voiceless", {}))
        tracks["folds"] = None if folds is None else (numpy.asarray(folds[0], float), numpy.asarray(folds[1], float))
        jet = jet_keyframes(arrivals, leavings, getattr(target, "fricatives", ()))
        tracks["jet"] = None if jet is None else (numpy.asarray(jet[0], float), numpy.asarray(jet[1], float))
        return tracks

    def resumable(self, target: Any) -> bool:
        """Whether this target's residual is sealed_core's (the one with a
        record and resume)."""
        return (RESUME and not self.nasals and getattr(target, "seal", 0.0) > 0.0 and target.equation
                and getattr(target, "output_from", None) is None)

    def with_bend(self, residual: numpy.ndarray, parameters: numpy.ndarray, target: Any) -> numpy.ndarray:
        bend = getattr(target, "bend", 0.0)
        if bend > 0.0:
            logs = parameters[:self.phones * self.count].reshape(self.phones, self.count)
            residual = numpy.concatenate([residual, bend * numpy.diff(logs, 2, axis=1).ravel()])
        return residual

    def resumed_residuals(self, base_point: numpy.ndarray, moved_points: Sequence[numpy.ndarray],
                          target: Any) -> Tuple[numpy.ndarray, List[numpy.ndarray]]:
        """The residual at base_point and at each moved point, the moved ones
        resumed from the base's record where their inputs first differ
        (column_locality's contract; observer.implied_resume skips the
        free render of every closure wholly before it). In the order given;
        computed in decreasing order of the first difference (the geometry
        workspace then copies nothing back)."""
        if not self.resumable(target):
            return (self.residual(base_point, target), [self.residual(point, target) for point in moved_points])
        tubes, arrivals, leavings = self.split(base_point)
        base_tracks = self.input_tracks(base_point, target)
        base_residual, record = sealed_core(tubes, arrivals, leavings, self.length, target, True, None, base_tracks)
        tracks = [self.input_tracks(point, target) for point in moved_points]
        order = sorted(range(len(moved_points)), key=DivergenceOrder(base_tracks, tracks), reverse=True)
        results: List[Optional[numpy.ndarray]] = [None] * len(moved_points)
        for position in order:
            point = moved_points[position]
            moved_tubes, moved_arrivals, moved_leavings = self.split(point)
            value, _ = sealed_core(moved_tubes, moved_arrivals, moved_leavings, self.length, target, False, record,
                                   tracks[position])
            results[position] = self.with_bend(value, point, target)
        return self.with_bend(base_residual, base_point, target), results

    def velum(self, parameters):
        """Each phone's velum port area (0 for an oral phone)."""
        areas = numpy.zeros(self.phones)
        start = self.phones * self.count + self.time_count
        for offset, index in enumerate(self.nasals):
            areas[index] = max(float(parameters[start + offset]), 0.0)
        return areas

    def residual(self, parameters, target):
        tubes, arrivals, leavings = self.split(parameters)
        if self.nasals and target.equation:
            residual = branched_residual(tubes, arrivals, leavings, self.velum(parameters), self.length, target)
        else:
            residual = sequence_residual(tubes, arrivals, leavings, self.length, target)
        bend = getattr(target, "bend", 0.0)
        if bend > 0.0:
            # Each tube's node-to-node bend (as tube_recover.Shape's).
            logs = parameters[:self.phones * self.count].reshape(self.phones, self.count)
            residual = numpy.concatenate([residual, bend * numpy.diff(logs, 2, axis=1).ravel()])
        return residual


# Finite-difference columns resumed from the base point's record (exact:
# driver_bench resume); env RESUME=0 computes every column whole.
RESUME = tube_recover.EXACT_SHORTCUTS and os.environ.get("RESUME", "1") == "1"


class DivergenceOrder:
    """Sort key: a moved point's first input divergence (s) from the base."""

    def __init__(self, base_tracks: Dict[str, Any], tracks: Sequence[Dict[str, Any]]) -> None:
        self.base_tracks = base_tracks
        self.tracks = tracks

    def __call__(self, position: int) -> float:
        moment = column_locality.first_divergence(self.base_tracks, self.tracks[position])
        return 1e9 if moment == numpy.inf else moment


class FreeFlow:
    """Our folds' flow for a sequence hypothesis, forced by given returning
    waves (NaN: free)."""

    def __init__(self, first, length, seconds, times, frames, seal, velum=None, folds=None, jet=0.0, jet_gate=None):
        self.folds = folds
        self.jet = jet
        self.jet_gate = jet_gate
        self.first = first
        self.length = length
        self.seconds = seconds
        self.times = times
        self.frames = frames
        self.seal = seal
        self.velum = velum

    def __call__(self, forced, observed=None, seconds=None):
        seconds = self.seconds if seconds is None else seconds
        flow, _ = posture_recover.render((self.first, self.length), seconds, keyframes=(self.times, self.frames),
                                         forced=forced, seal=self.seal, velum=self.velum, folds=self.folds,
                                         jet=self.jet, jet_gate=self.jet_gate,
                                         observed=None if observed is None else (observed, OBSERVE_GAIN))
        if len(flow) < len(forced):
            flow = numpy.concatenate([flow, numpy.zeros(len(forced) - len(flow))])
        return flow[:len(forced)]


SCRATCH = posture_recover.SCRATCH   # this tree's scratch (the reference's is /tmp/es)
# Renders stop at the compared stretch's end (exact: the engine is causal
# and its output does not depend on the run's length; a 0.41 s /aa b aa/
# render and a 0.2 s one agreed bitwise over the 0.2 s; driver_bench.py
# checks each residual). A low-passed target (filtfilt reads both ways)
# keeps the whole sound.
TRUNCATE = tube_recover.EXACT_SHORTCUTS


def horizon(target: Any) -> int:
    """The steps a residual's renders must reach for the target's compared
    stretch."""
    if TRUNCATE and target.filter is None:
        return min(target.start + target.length, len(target.flow_out))
    return len(target.flow_out)


def reach(target: Any, nodes: int) -> int:
    """The steps of the sound a residual's peels must read for the target's
    compared stretch: its horizon and the peels' reach past it
    (column_locality.LOOKAHEAD_NODES node counts; a 46-node tube's peel was
    bitwise the same on the sound cut 50 steps past the stretch, not 25).
    The whole sound when not truncating."""
    if TRUNCATE and target.filter is None:
        return min(horizon(target) + column_locality.LOOKAHEAD_NODES * nodes, len(target.flow_out))
    return len(target.flow_out)


NOSE = None   # the reversible nose's constants (reverse_core.read_nose), read once


def nose_constants():
    global NOSE
    if NOSE is None:
        path = os.path.join(SCRATCH, "nose_constants_%d.txt" % os.getpid())
        posture_recover.render((numpy.full(46, 2.0), 46.5 * posture_recover.NODE), 0.01, velum=([0.0, 1.0], [0.0, 0.0]),
                               nose_dump=path)
        NOSE = reverse_core.read_nose(path)
        os.remove(path)
    return NOSE


def velum_keyframes(arrivals, leavings, velum_areas):
    """The velum's keyframes: each phone's port area held from its arrival
    to its leaving (0 for an oral phone), linear between."""
    times = []
    areas = []
    for arrive, leave, area in zip(arrivals, leavings, velum_areas):
        times.extend([arrive, leave])
        areas.extend([area, area])
    return times, areas


RELEASE_MARGIN = 10   # steps either side of a nasal release's window
# With the velum open, a mouth narrower than this in front of the port is
# the weak path: the march (oral peeled, nose forward) inverts it inside
# the loop and ran away (/uw n ih ao/: 341 dB at the truth, /uw/'s lips
# 0.333 cm^2), so the folds decide there too.
WEAK_ORAL = 0.7       # cm^2 (0.5 left /iy uw n/ running away: /uw/ at /iy/'s count narrows to 0.518 just in front of the port)


def march_residual(tubes, arrivals, leavings, velum_areas, length, target):
    """The equation error of a nasal word: the split march (marching) as the
    inverse wherever the sound decides it, and at each release of a seal
    with the velum open a short window decided by our folds (a render free
    there, its PORTWAVES driving the march through; the sound cannot decide
    those ~nose-length steps). Inside the windows, our sound against the
    sound."""
    times, frames = keyframes(tubes, arrivals, leavings)
    count = reach(target, frames.shape[1])
    sound = target.flow_out[:count]
    schedule = reverse_core.keyframe_schedule(times, frames, count)
    geometry = numpy.ascontiguousarray(reverse_core.moving_geometry(schedule, posture_recover.VISCOUS))
    velum = velum_keyframes(arrivals, leavings, velum_areas)
    port = reverse_core.velum_port(velum[0], velum[1], count)
    nose = nose_constants()
    port_node = nose[0]
    m = len(nose[1])
    end_extra = length - frames.shape[1] * posture_recover.NODE
    windows = []
    if HIDE_AREA > target.seal:
        # The folds also decide where a node in front of the port is narrow
        # but open (closing and opening): the oral peel through it amplified
        # a wrong hypothesis's mismatch until the march ran away (/aa m aa/,
        # the closure 0.2 ms late: +27.9 dB, 0.5 ms: +1580 dB). Closing:
        # from when the narrowing enters the oral peel's reach (port node at
        # t - 1, open end at t - 1 + oral_lag) to when the seal does (the
        # sealed mode then peels the nose); opening: from the release's
        # window to past the narrowing.
        oral_lag = frames.shape[1] - 1 - port_node
        front = numpy.zeros_like(geometry[:, :, 0])
        front[:, port_node + 1:] = geometry[:, port_node + 1:, 0]
        front[:, :port_node + 1] = numpy.inf
        narrow = (front.min(axis=1) <= HIDE_AREA) | ((front.min(axis=1) <= WEAK_ORAL) & (port >= 1e-6))
        shut = front.min(axis=1) <= target.seal
        t = 0
        while t < count:
            if not narrow[t]:
                t += 1
                continue
            start = t
            while t < count and narrow[t]:
                t += 1
            stretch = numpy.nonzero(shut[start:t])[0]
            begin = max(start - oral_lag + 1 - RELEASE_MARGIN, 1)
            if len(stretch) == 0 or numpy.any(port[start + stretch] < 1e-6):
                # Narrow and never shut, or shut with the velum shut for any
                # of it (the sealed mode peels the nose: /n d/ held one seal
                # while the velum closed, and the march fell back to the
                # oral peel through it): all of it the folds'.
                windows.append((begin, min(t + 2 + RELEASE_MARGIN, count)))
                continue
            first_shut = start + int(stretch[0])
            last_shut = start + int(stretch[-1])
            closing_end = min(first_shut - oral_lag + 1 + RELEASE_MARGIN, count)
            if closing_end > begin:
                windows.append((begin, closing_end))
            windows.append((max(last_shut - m + 2 - RELEASE_MARGIN, 1), min(t + 2 + RELEASE_MARGIN, count)))
        merged = []
        for begin, end in sorted(windows):
            if merged and begin <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((begin, end))
        windows = merged
    else:
        for first, last, node in observer.closures(geometry, target.seal):
            if last + 1 < count and port[last + 1] >= 1e-6:
                windows.append((max(last - m + 2 - RELEASE_MARGIN, 1), min(last + 2 + RELEASE_MARGIN, count)))
    if os.environ.get("SHOW_WINDOWS"):
        print("march windows (s):", " ".join("%.4f-%.4f" % (b / STEP_RATE, e / STEP_RATE) for b, e in windows))
    seconds = (count + 1) / STEP_RATE
    area0 = geometry[:, 0, 0]
    port_y = port / reverse_core.RHO_C
    gate = jet_gate_for(arrivals, leavings, target.fricatives, count)
    jet = (target.jet, reverse_core.JET_SEED, posture_recover.JET_CLOSING if target.seal > 0.0 else 0.0, gate)
    jet_steps = marching.jet_table(numpy.ascontiguousarray(geometry[:, :, 0]), jet, gate)

    def returning(arrays, until=None):
        flow, pressure = marching._peel_below(geometry, port_node, arrays["port_pressure"], arrays["from_front"],
                                              arrays["returned"], port_y, posture_recover.HOP, True, jet_steps,
                                              until=until)
        return flow, pressure, 0.5 * (pressure - flow * reverse_core.RHO_C / area0)

    def render_waves(forced, steps=None):
        path = os.path.join(SCRATCH, "waves_%d_%d.raw" % (os.getpid(), numpy.random.randint(1 << 30)))
        span = seconds if steps is None else (steps + 1) / STEP_RATE
        flow, _ = posture_recover.render((frames[0], length), span, keyframes=(times, frames), forced=forced,
                                         seal=target.seal, velum=velum, port_waves=path, jet=target.jet,
                                         jet_gate=jet_keyframes(arrivals, leavings, target.fricatives),
                                         folds=fold_keyframes(arrivals, leavings, target.voiceless))
        waves = numpy.fromfile(path, dtype=numpy.float64).reshape(-1, 6)
        os.remove(path)
        return flow, waves

    # The peel below the port reads up to port_node steps ahead: the waves
    # it implies are good only that far short of what the march has
    # reached. Our folds run free from there (forced with the last steps'
    # waves, computed from zeros past the window, the window was lost).
    free = numpy.zeros(count, bool)
    for begin, end in windows:
        free[max(begin - port_node - 2, 0):end] = True

    # A window's render is read up to its end plus the nose's lag
    # (marching's windowed step reads share at k = t + m - 1).
    ends = {begin: end for begin, end in windows}

    def render_window(arrays, begin):
        # (The peel below the port only as far as the forced waves are read:
        # NaN from there; NOTES_kernels.md's recipe.)
        flow, pressure, wave = returning(arrays, max(begin - port_node - 2, 0) if TRUNCATE else None)
        forced = wave.copy()
        forced[max(begin - port_node - 2, 0):] = numpy.nan
        if TRUNCATE:
            return render_waves(forced, min(ends[begin] + m + 2, count))[1]
        return render_waves(forced)[1]

    share, returned, port_pressure, from_front, mode, waves, in_window = marching.march_windows(
        sound, geometry, port, nose, posture_recover.HOP, end_extra, target.seal, windows, render_window,
        posture_recover.NOSE_POLE, posture_recover.NOSE_ZERO, jet)
    arrays = {"port_pressure": port_pressure, "from_front": from_front, "returned": returned}
    implied, pressure, wave = returning(arrays)
    forced = wave.copy()
    forced[free] = numpy.nan
    driven, waves_b = render_waves(forced)
    size = min(len(driven), len(implied), len(waves_b))
    difference = driven[:size] - implied[:size]
    hidden = free[:size]
    difference[hidden] = waves_b[:size, 4][hidden] + waves_b[:size, 5][hidden] - target.flow_out[:size][hidden]
    if target.hidden_only:
        difference[~hidden] = 0.0
    if target.filter is not None:
        difference = tube_recover.sosfiltfilt(target.filter, difference)
    return tube_recover.fixed_clip(difference, target)


def branched_residual(tubes, arrivals, leavings, velum_areas, length, target):
    if NASAL_INVERSE == "march":
        return march_residual(tubes, arrivals, leavings, velum_areas, length, target)
    return window_residual(tubes, arrivals, leavings, velum_areas, length, target)


def window_residual(tubes, arrivals, leavings, velum_areas, length, target):
    """The equation error where the sound reveals the inside; inside the
    windows (a seal or the velum open) the forward sound against the sound."""
    times, frames = keyframes(tubes, arrivals, leavings)
    count = len(target.flow_out)
    schedule = reverse_core.keyframe_schedule(times, frames, count)
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    velum = velum_keyframes(arrivals, leavings, velum_areas)
    port = reverse_core.velum_port(velum[0], velum[1], count)
    free = FreeFlow(frames[0], length, (len(target.flow_out) + 1) / STEP_RATE, times, frames, target.seal, velum,
                    fold_keyframes(arrivals, leavings, target.voiceless))
    implied, pressure, revealed, modelled = observer.implied_branched(
        geometry, target.flow_out, port, nose_constants(), target.seal, posture_recover.HOP,
        length - frames.shape[1] * posture_recover.NODE, free, posture_recover.NOSE_POLE, posture_recover.NOSE_ZERO)
    forced = 0.5 * (pressure - implied * reverse_core.RHO_C / geometry[:, 0, 0])
    forced[~revealed] = numpy.nan
    driven = free(forced)
    size = min(len(driven), len(implied))
    difference = driven[:size] - implied[:size]
    hidden = ~revealed[:size]
    difference[hidden] = modelled[:size][hidden] - target.flow_out[:size][hidden]
    if target.filter is not None:
        difference = tube_recover.sosfiltfilt(target.filter, difference)
    return tube_recover.fixed_clip(difference, target)


class SealedRecord:
    """A sealed residual's base run kept for its finite-difference columns
    (observer.implied_record, reverse_core.moving_geometry_record): the
    schedule, the input tracks' times, the geometry workspace and the
    observer's record."""

    def __init__(self, schedule: numpy.ndarray, tracks: Dict[str, Any], workspace: Any, implied: Any) -> None:
        self.schedule = schedule
        self.tracks = tracks
        self.workspace = workspace
        self.implied = implied


def sealed_residual(tubes, arrivals, leavings, length, target):
    """The equation error with closures (observer.implied): zero where the
    sound cannot reveal the inside and the folds run free."""
    return sealed_core(tubes, arrivals, leavings, length, target)[0]


def last_hidden_end(geometry: numpy.ndarray, seal: float, count: int) -> int:
    """The last closure's hidden stretch's end (observer.implied's
    hidden_to = last sealed step + 2), 0 with none."""
    end = 0
    for first, last, node in observer.closures(geometry, seal, HIDE_AREA):
        end = max(end, min(last + 2, count))
    return end


def sealed_core(tubes, arrivals, leavings, length, target, keep: bool = False, base: Optional[SealedRecord] = None,
                tracks: Optional[Dict[str, Any]] = None) -> Tuple[numpy.ndarray, Optional[SealedRecord]]:
    """sealed_residual, with its record kept (keep) or resumed from a base
    record (base; `tracks`, this hypothesis's input tracks, for where it
    first differs). The same numbers either way (driver_bench resume)."""
    times, frames = keyframes(tubes, arrivals, leavings)
    steps = reach(target, frames.shape[1])
    sound = target.flow_out[:steps]
    schedule = reverse_core.keyframe_schedule(times, frames, steps)
    first_changed = 0
    if base is not None:
        first_changed = reverse_core.first_difference(schedule, base.schedule)
        # The folds' and the jet's keyframes reach the engine and the peels
        # too (the areas' are in the schedule).
        for name in ("folds", "jet", "velum"):
            moment = column_locality.track_divergence(base.tracks.get(name), tracks.get(name))
            if moment < numpy.inf:
                first_changed = min(first_changed, int(numpy.floor(moment * STEP_RATE)))
        if first_changed < base.workspace.record.every:
            # Changed before the geometry record's first boundary: computed
            # whole. (GeometryWorkspace.resume from boundary 0 starts from
            # the base's own start state, made from its first rest row; a
            # first phone's column changes that row: MEASURED, columns of
            # /aa b aa/'s first /aa/ 0.0016-0.024 off.)
            base = None
            geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
        else:
            geometry = base.workspace.resume(schedule, first_changed)
    elif keep:
        geometry, geometry_record = reverse_core.moving_geometry_record(schedule, posture_recover.VISCOUS)
    else:
        geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    # (The free renders span the whole sound: a closure past the window's end
    # still needs our folds' flow. Truncated, they span the compared stretch,
    # and a closure's free render only to the last closure's hidden end
    # (observer.implied reads a closure's folds over [hidden_from,
    # hidden_to); which closure a call is for is not passed, so the last
    # one's end bounds them all): /aa b aa/'s closure ends near 0.22 s of
    # the 0.40 s stretch.)
    closure_seconds = (horizon(target) + 1) / STEP_RATE
    if TRUNCATE:
        closure_seconds = (min(max(last_hidden_end(geometry, target.seal, steps), 1), horizon(target)) + 1) / STEP_RATE
    free = FreeFlow(frames[0], length, closure_seconds, times, frames, target.seal,
                    folds=fold_keyframes(arrivals, leavings, target.voiceless), jet=target.jet,
                    jet_gate=jet_keyframes(arrivals, leavings, target.fricatives))
    gate = jet_gate_for(arrivals, leavings, target.fricatives, steps)
    end_extra = length - frames.shape[1] * posture_recover.NODE
    record = None
    if base is not None:
        implied, pressure, revealed, modelled, implied_record = observer.implied_resume(
            base.implied, geometry, sound, free, first_changed, gate)
    elif keep:
        implied, pressure, revealed, modelled, implied_record = observer.implied_record(
            geometry, sound, target.seal, posture_recover.HOP, end_extra, free, hide=HIDE_AREA, jet=target.jet,
            jet_closing=posture_recover.JET_CLOSING, jet_gate=gate)
        record = SealedRecord(schedule, tracks if tracks is not None else {},
                              reverse_core.GeometryWorkspace(geometry, geometry_record), implied_record)
    else:
        implied, pressure, revealed, modelled = observer.implied(geometry, sound, target.seal, posture_recover.HOP,
                                                                 end_extra, free, hide=HIDE_AREA, jet=target.jet,
                                                                 jet_closing=posture_recover.JET_CLOSING,
                                                                 jet_gate=gate)
    forced = 0.5 * (pressure - implied * reverse_core.RHO_C / geometry[:, 0, 0])
    forced[~revealed] = numpy.nan
    observed = None
    if OBSERVE_GAIN > 0.0:
        observed = implied.copy()
        observed[~revealed] = numpy.nan
    driven = free(forced, observed, (horizon(target) + 1) / STEP_RATE)
    size = min(len(driven), len(implied))
    difference = driven[:size] - implied[:size]
    # Inside a closure's window the forward sound against the sound: the
    # release is the loudest thing a stop does. (Left out, the release's
    # time was 7 ms off and the vowel after it was lost.)
    hidden = ~revealed[:size]
    difference[hidden] = modelled[:size][hidden] - target.flow_out[:size][hidden]
    if target.hidden_only:
        difference[~hidden] = 0.0
    if target.filter is not None:
        difference = tube_recover.sosfiltfilt(target.filter, difference)
    return tube_recover.fixed_clip(difference, target), record


def sequence_residual(tubes, arrivals, leavings, length, target):
    if getattr(target, "seal", 0.0) > 0.0 and target.equation:
        return sealed_residual(tubes, arrivals, leavings, length, target)
    """The equation error of a sequence hypothesis over the target's window
    (or, for a target without `equation`, our folds running free)."""
    times, frames = keyframes(tubes, arrivals, leavings)
    schedule = reverse_core.keyframe_schedule(times, frames, len(target.flow_out))
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    implied, pressure = reverse_core.peel_moving_glottis(geometry, target.flow_out, hop=posture_recover.HOP,
                                                         end_extra=length - frames.shape[1] * posture_recover.NODE,
                                                         jet=target.jet,
                                                         gate=jet_gate_for(arrivals, leavings, target.fricatives,
                                                                           len(target.flow_out)))
    forced = None
    if target.equation:
        forced = 0.5 * (pressure - implied * reverse_core.RHO_C / geometry[:, 0, 0])
    driven, _ = posture_recover.render((frames[0], length), target.seconds, keyframes=(times, frames), forced=forced,
                                       jet=target.jet, folds=fold_keyframes(arrivals, leavings, target.voiceless),
                                      jet_gate=jet_keyframes(arrivals, leavings, target.fricatives))
    size = min(len(driven), len(implied))
    difference = driven[:size] - implied[:size]
    if target.filter is not None:
        difference = tube_recover.sosfiltfilt(target.filter, difference)
    return tube_recover.fixed_clip(difference, target)


class PrefixContext:
    """A phone's candidate tube heard in its sequence: the phones before it
    as recovered, it held from its arrival on (to the end), the times as
    they stand, the length fixed."""

    def __init__(self, before, arrivals, leavings, length, seal=0.0, velum_areas=None, jet=0.0, voiceless=None,
                 leading=None, fricatives=()):
        self.fricatives = fricatives
        # A word's leading consonants ({index: (node, area)}), made from the
        # candidate itself (the first vowel's): `before` is then empty.
        self.leading = leading if leading is not None else {}
        self.jet = jet
        self.voiceless = voiceless if voiceless is not None else {}
        self.before = [numpy.asarray(tube, float) for tube in before]
        self.arrivals = arrivals
        self.leavings = leavings
        self.length = length
        self.seal = seal
        self.velum_areas = velum_areas   # per phone; None: no nose

    def tubes(self, tube):
        rest = len(self.arrivals) - len(self.before)
        tubes = numpy.array(self.before + [numpy.asarray(tube[0], float)] * rest)
        for index, (node, area) in self.leading.items():
            tubes[index, node] = area
        return tubes

    def sound(self, tube, seconds):
        times, frames = keyframes(self.tubes(tube), self.arrivals, self.leavings)
        velum = None
        if self.velum_areas is not None:
            velum = velum_keyframes(self.arrivals, self.leavings, self.velum_areas)
        _, sound = posture_recover.render((frames[0], self.length), seconds, keyframes=(times, frames), seal=self.seal,
                                          velum=velum, jet=self.jet,
                                          jet_gate=jet_keyframes(self.arrivals, self.leavings, self.fricatives),
                                          folds=fold_keyframes(self.arrivals, self.leavings, self.voiceless))
        return sound

    def equation(self, tube, target):
        if self.velum_areas is not None:
            return branched_residual(self.tubes(tube), self.arrivals, self.leavings, self.velum_areas, self.length,
                                     target)
        return sequence_residual(self.tubes(tube), self.arrivals, self.leavings, self.length, target)

    def implied(self, tube, target):
        """The flow entering, for the error's scale (with a seal or the nose:
        our folds' own, run free; the peel cannot pass those)."""
        times, frames = keyframes(self.tubes(tube), self.arrivals, self.leavings)
        if self.seal > 0.0 or self.velum_areas is not None:
            velum = None if self.velum_areas is None else velum_keyframes(self.arrivals, self.leavings,
                                                                           self.velum_areas)
            flow, _ = posture_recover.render((frames[0], self.length), target.seconds, keyframes=(times, frames),
                                             seal=self.seal, velum=velum, jet=self.jet,
                                             jet_gate=jet_keyframes(self.arrivals, self.leavings, self.fricatives),
                                             folds=fold_keyframes(self.arrivals, self.leavings, self.voiceless))
            return flow
        schedule = reverse_core.keyframe_schedule(times, frames, len(target.flow_out))
        geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
        return reverse_core.peel_moving(geometry, target.flow_out, hop=posture_recover.HOP,
                                        end_extra=self.length - frames.shape[1] * posture_recover.NODE, jet=self.jet,
                                        gate=jet_gate_for(self.arrivals, self.leavings, self.fricatives,
                                                          len(target.flow_out)))


def level(residual, reference):
    """A residual's level against the reference flow's swing (dB)."""
    return 10.0 * numpy.log10(numpy.sum(residual ** 2) / numpy.sum((reference - reference.mean()) ** 2))


def window(begin, end):
    return int(begin * STEP_RATE), int((end - begin) * STEP_RATE)


def output_residual(tubes, arrivals, leavings, velum_areas, length, target):
    """Our whole engine's sound against the sound over the target's stretch:
    our folds forced by the waves the sound implies up to target.output_from
    (a step before anything there hides the inside: the plain peel holds),
    free from there. A pure forward run from a state the sound gave, the
    same stretch for every hypothesis (a release's time: the folds' window
    alone grew with a later release, and its minimum slid 11 ms early)."""
    times, frames = keyframes(tubes, arrivals, leavings)
    begin = target.output_from
    # The peel only as far as the forced waves are wanted (it reads ~a node
    # count of steps ahead).
    peel_reach = min(max(begin, 1) + 4 * frames.shape[1], len(target.flow_out))
    schedule = reverse_core.keyframe_schedule(times, frames, peel_reach)
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    implied, pressure = reverse_core.peel_moving_glottis(geometry, target.flow_out[:peel_reach], hop=posture_recover.HOP,
                                                         end_extra=length - frames.shape[1] * posture_recover.NODE,
                                                         jet=target.jet,
                                                         gate=jet_gate_for(arrivals, leavings, target.fricatives,
                                                                           len(target.flow_out)))
    forced = numpy.full(len(target.flow_out), numpy.nan)
    forced[:begin] = (0.5 * (pressure - implied * reverse_core.RHO_C / geometry[:, 0, 0]))[:begin]
    velum = None
    if velum_areas is not None and numpy.any(numpy.asarray(velum_areas) > 0.0):
        velum = velum_keyframes(arrivals, leavings, velum_areas)
    _, sound = posture_recover.render((frames[0], length), (horizon(target) + 1) / STEP_RATE,
                                      keyframes=(times, frames), forced=forced, seal=target.seal, velum=velum,
                                      jet=target.jet, folds=fold_keyframes(arrivals, leavings, target.voiceless),
                                      jet_gate=jet_keyframes(arrivals, leavings, target.fricatives))
    size = min(len(sound), len(target.flow_out))
    difference = sound[:size] - target.flow_out[:size]
    return tube_recover.fixed_clip(difference, target)


def word_residual(tubes, arrivals, leavings, velum_areas, length, target):
    """The residual of a word's hypothesis, with the nose when any phone
    opens the velum."""
    if target.output_from is not None:
        return output_residual(tubes, arrivals, leavings, velum_areas, length, target)
    if velum_areas is not None and numpy.any(numpy.asarray(velum_areas) > 0.0) and target.equation:
        return branched_residual(tubes, arrivals, leavings, velum_areas, length, target)
    return sequence_residual(tubes, arrivals, leavings, length, target)


def timing_cost(task):
    tubes, arrivals, leavings, velum_areas, length, target = task
    return float(numpy.sum(word_residual(tubes, arrivals, leavings, velum_areas, length, target) ** 2))


class GridMemo:
    """The timing grids' costs within one transition's refinement (the tubes
    fixed): each rung's centre is the rung before's best, evaluated again on
    the same target (exact: the same bytes in, the same cost out)."""

    def __init__(self) -> None:
        self.costs: Dict[Tuple[bytes, bytes, bytes, int], float] = {}
        self.targets: List[Any] = []   # kept alive: the keys hold their ids

    def key(self, arrivals: numpy.ndarray, leavings: numpy.ndarray, velum_areas: Any, target: Any) -> Tuple[bytes, bytes, bytes, int]:
        if all(kept is not target for kept in self.targets):
            self.targets.append(target)
        velum = b"" if velum_areas is None else numpy.asarray(velum_areas, float).tobytes()
        return (numpy.asarray(arrivals, float).tobytes(), numpy.asarray(leavings, float).tobytes(), velum, id(target))


def timing_grid(pool, tubes, arrivals, leavings, length, index, leaving_offsets, arrival_offsets, target,
                velum_areas=None, memo: Optional[GridMemo] = None):
    """Transition `index` (into phone `index`): the best leaving and arrival
    over a grid of offsets from the current ones."""
    tasks = []
    for leaving_offset in leaving_offsets:
        for arrival_offset in arrival_offsets:
            trial_leavings = leavings.copy()
            trial_arrivals = arrivals.copy()
            trial_leavings[index - 1] += leaving_offset
            trial_arrivals[index] += arrival_offset
            if trial_arrivals[index] > trial_leavings[index - 1] + GAP:
                tasks.append((tubes, trial_arrivals, trial_leavings, velum_areas, length, target))
    costs = numpy.zeros(len(tasks))
    pending: List[int] = []
    for position, task in enumerate(tasks):
        if tube_recover.EXACT_SHORTCUTS and memo is not None and memo.key(task[1], task[2], velum_areas, target) in memo.costs:
            costs[position] = memo.costs[memo.key(task[1], task[2], velum_areas, target)]
            driver_ledger.count_saved(1)
        else:
            pending.append(position)
    values = driver_ledger.pool_map(pool, timing_cost, [tasks[position] for position in pending])
    for position, value in zip(pending, values):
        costs[position] = value
        if memo is not None:
            memo.costs[memo.key(tasks[position][1], tasks[position][2], velum_areas, target)] = value
    costs[~numpy.isfinite(costs)] = numpy.inf
    best = tasks[int(numpy.argmin(costs))]
    return best[1], best[2]


JOINT_STOP_DB = 0.5   # FAST_ALGO jointstop: an unbent round gaining less than this ends its rounds
# FAST_ALGO jointmore: the final joint goes on past its 4 rounds while a round
# gains at least JOINT_MORE_DB, to JOINT_ROUNDS_MAX. (MEASURED, /aa b aa/ on
# the corrected model, joint_continue.py from the 4-round result: rounds 5-8
# gained 2.9, 2.3, 2.9 and 0.9 dB; round trip 130.0 -> 133.5, 136.3, 139.6,
# 140.9 dB; tubes rms max 0.0077 -> 0.0029 cm^2.)
JOINT_MORE_DB = 1.0
JOINT_ROUNDS_MAX = 8
# FAST_ALGO blocks: after the final joint, rounds of the times alone (their
# finite-difference step BLOCK_TIME_STEP) then the tubes and velum alone,
# while a round gains at least BLOCK_STOP_DB, to BLOCK_ROUNDS_MAX. (MEASURED,
# Echo on the corrected model from its stalled joint: the times' columns
# were ~2000x the tubes' (1.2e7 against ~5e3) and nonlinear inside the 2e-5
# s step (the voiceless closure's abduction ends the voicing: a bifurcation),
# so every whole step failed; alternating, with the times stepped by 1e-6 s:
# round trip 48.4 -> 53.2, 56.2, 59.3 dB in three rounds, still climbing.)
BLOCK_TIME_STEP = 1e-6
# The times' step then shrinks with their moves (BLOCK_STEP_SHARE of the
# last round's largest time move, within BLOCK_STEP_FLOOR..BLOCK_TIME_STEP).
# (MEASURED, Golf polished at -107 dB, truth -245.9: toward the truth the
# cost fell 0.0137 -> 0.0096 over 75% of the way and five orders over the
# last 10%, ~6 ns of time: its word-leading /g/'s release was 60 ns off and
# the voicing's start from rest made that an even 83 dB floor from 0.10 s
# on; the release's column was nonlinear inside the 1e-6 s step.)
BLOCK_STEP_SHARE = 0.1
BLOCK_STEP_FLOOR = 1e-9
BLOCK_STOP_DB = 0.3
BLOCK_ROUNDS_MAX = 12


class BlockShape:
    """A shape's parameters with all but one block held: the free block's
    indices, the rest from `held`."""

    def __init__(self, shape: Any, held: numpy.ndarray, free: numpy.ndarray, steps: numpy.ndarray) -> None:
        self.shape = shape
        self.held = numpy.array(held, float)
        self.free = numpy.asarray(free, int)
        self.free_steps = numpy.asarray(steps, float)

    def full(self, parameters: numpy.ndarray) -> numpy.ndarray:
        out = self.held.copy()
        out[self.free] = parameters
        return out

    def residual(self, parameters: numpy.ndarray, target: Any) -> numpy.ndarray:
        return self.shape.residual(self.full(parameters), target)

    def steps(self, size: int) -> numpy.ndarray:
        return self.free_steps

    def inert_columns(self) -> List[int]:
        """The shape's inert columns that are free here, as block indices
        (a stop's own areas never reach the sound)."""
        inert = set(self.shape.inert_columns()) if hasattr(self.shape, "inert_columns") else set()
        return [position for position, index in enumerate(self.free) if int(index) in inert]

    def inert_residual(self, parameters: numpy.ndarray, target: Any, base: numpy.ndarray) -> numpy.ndarray:
        return self.shape.inert_residual(self.full(parameters), target, base)


SWEEP_SCALES = (1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-9)   # s
SWEEP_REPEATS = 3


def sweep_cost(task: Tuple[Any, numpy.ndarray, Any]) -> float:
    shape, parameters, target = task
    return float(numpy.sum(shape.residual(parameters, target) ** 2))


def time_sweep(pool: Any, shape: Any, parameters: numpy.ndarray, whole: Any) -> Tuple[numpy.ndarray, float]:
    """The times one at a time by +-h at shrinking scales (SWEEP_SCALES),
    every candidate of a scale evaluated together: each time moves to its
    better side when that lowers the cost; the moves are kept together when
    they lower it together, else only the best one. Repeated at a scale while
    it gains (to SWEEP_REPEATS). Returns the parameters and the largest time
    move. (LM's model cannot see a canyon narrower than its step: Golf's
    word-leading release was 60 ns off, and at a 1 us step LM moved the
    times 16 ps.)"""
    first = shape.phones * shape.count
    indices = list(range(first, first + shape.time_count))
    start = numpy.array(parameters, float)
    current = start.copy()
    cost = sweep_cost((shape, current, whole))
    for scale in SWEEP_SCALES:
        for repeat in range(SWEEP_REPEATS):
            tasks = []
            for index in indices:
                for sign in (-1.0, 1.0):
                    moved = current.copy()
                    moved[index] += sign * scale
                    tasks.append((shape, moved, whole))
            costs = driver_ledger.pool_map(pool, sweep_cost, tasks)
            moves = numpy.zeros(len(current))
            best_single = None
            for position, index in enumerate(indices):
                lower, upper = costs[2 * position], costs[2 * position + 1]
                if min(lower, upper) < cost:
                    moves[index] = -scale if lower < upper else scale
                    single = min(lower, upper)
                    if best_single is None or single < best_single[0]:
                        best_single = (single, index, moves[index])
            if best_single is None:
                break
            combined = current + moves
            combined_cost = sweep_cost((shape, combined, whole))
            if combined_cost < best_single[0]:
                current, cost = combined, combined_cost
            else:
                current = current.copy()
                current[best_single[1]] += best_single[2]
                cost = best_single[0]
    return current, float(numpy.max(numpy.abs(current - start)))


def block_rounds(pool: Any, shape: Any, parameters: numpy.ndarray, whole: Any, reference: numpy.ndarray,
                 report: Any) -> numpy.ndarray:
    """FAST_ALGO blocks: the times alone, then the tubes and velum alone,
    round after round while a round gains BLOCK_STOP_DB."""
    size = len(parameters)
    steps = shape.steps(size)
    times = numpy.arange(shape.phones * shape.count, shape.phones * shape.count + shape.time_count)
    rest = numpy.array([k for k in range(size) if k < times[0] or k > times[-1]])
    before = level(shape.residual(parameters, whole), reference)
    time_step = BLOCK_TIME_STEP
    for round_index in range(BLOCK_ROUNDS_MAX):
        moved_times = 0.0
        for label, free, free_steps in (("times", times, numpy.full(len(times), time_step)),
                                        ("tubes", rest, steps[rest])):
            if label == "times":
                with driver_ledger.stage("block %d times" % round_index):
                    parameters, moved_times = time_sweep(pool, shape, parameters, whole)
                continue
            block = BlockShape(shape, parameters, free, free_steps)
            with driver_ledger.stage("block %d %s" % (round_index, label)):
                found, used = tube_recover.stage(pool, block, parameters[free], whole, 30)
            parameters = block.full(found)
        moved = moved_times >= 1e-12   # (LM's last-place moves are not moves)
        if moved:
            time_step = min(max(BLOCK_STEP_SHARE * moved_times, BLOCK_STEP_FLOOR), BLOCK_TIME_STEP)
        else:
            # The times did not move at this step (a canyon narrower than
            # it): ten times finer.
            time_step = max(0.1 * time_step, BLOCK_STEP_FLOOR)
        after = level(shape.residual(parameters, whole), reference)
        report("block round %d" % round_index, parameters, after)
        print("  (block round %d: times moved at most %.3g s; next time step %.3g s)" % (
            round_index, moved_times, time_step))
        if before - after < BLOCK_STOP_DB:
            break
        before = min(before, after)
    return parameters


class StillRounds:
    """Rounds of Levenberg-Marquardt where a round that did not move makes
    the next round on the same target an exact repeat (deterministic
    residuals; the stage starts afresh from the same bytes): skipped, its
    evaluations counted as saved. (/aa m aa/'s and Echo's four final joint
    rounds did not move: 4 x 314 evaluations, three of them repeats.)
    With FAST_ALGO jointstop, an unbent round that gained less than
    JOINT_STOP_DB on its own residual also ends the rounds on that target
    (an algorithmic stop: the reference's /aa b aa/ final rounds gained 6.2,
    5.2, 1.9 and 2.1 dB, so it does not stop there)."""

    def __init__(self) -> None:
        self.still: Optional[Tuple[Any, bytes, int]] = None   # (target signature, parameters' bytes, evaluations)
        self.stopped: Optional[Any] = None   # a signature whose rounds jointstop ended
        self.reason = ""
        self.gain = 0.0   # the last round's (dB)

    def run(self, pool: Any, shape: Any, parameters: numpy.ndarray, target: Any, signature: Any,
            budget: int) -> Tuple[numpy.ndarray, int, bool]:
        start = numpy.array(parameters, float)
        if (tube_recover.EXACT_SHORTCUTS and self.still is not None and self.still[0] == signature
                and self.still[1] == start.tobytes()):
            driver_ledger.count_saved(self.still[2])
            self.reason = "an exact repeat of a round that did not move"
            return start, 0, True
        if self.stopped is not None and self.stopped == signature:
            self.reason = "jointstop: the round before gained less than %.1f dB" % JOINT_STOP_DB
            return start, 0, True
        found, used, gain = tube_recover.stage_gain(pool, shape, start, target, budget)
        found = numpy.array(found, float)
        self.still = (signature, start.tobytes(), used) if found.tobytes() == start.tobytes() else None
        if "jointstop" in tube_recover.ALGO and getattr(target, "bend", 0.0) == 0.0 and gain < JOINT_STOP_DB:
            self.stopped = signature
        self.gain = gain
        return found, used, False


def prefix_joint(pool, tubes, arrivals, leavings, length, flow_out, end, true_tubes, true_arrivals, true_leavings,
                 stops=None, velum_areas=None, narrowings=None):
    """Every phone so far and the times between them, together, over the
    sound up to `end` (s): the next phone is then heard after a history as
    exact as can be made. (Phone by phone alone, the errors piled up through
    our folds' phase: /aa/ exact, /iy/ -96.4 dB, /ae/ -57.4, /uh/ -23.5.)"""
    phones = len(tubes)
    nasal = [k for k in range(phones) if velum_areas is not None and velum_areas[k] > 0.0]
    shape = SequenceShape(phones, len(tubes[0]), length,
                          {k: node for k, node in (stops or {}).items() if k < phones}, nasal, narrowings)
    start, size = window(0.02, end)
    parameters = shape.join(tubes, arrivals[:phones], leavings[:phones],
                            None if velum_areas is None else velum_areas[:phones])
    exact = tube_recover.Target(flow_out, start, size, 0.0, 1.0, True)
    initial = parameters.copy()
    initial_cost = float(numpy.sum(shape.residual(initial, exact) ** 2))
    rounds = StillRounds()
    for round_index, bend in enumerate(JOINT_BENDS):
        target = tube_recover.Target(flow_out, start, size, 0.0, 1.0, True, bend)
        with driver_ledger.stage("round %d" % round_index):
            parameters, used, skipped = rounds.run(pool, shape, parameters, target, bend, 30)
    reference = shape.residual(parameters, exact)
    if float(numpy.sum(reference ** 2)) > initial_cost or not joint_sane(shape, initial, parameters):
        # The bent rounds traded the sound for smoothness and the unbent
        # ones did not win it back: kept as it was. (Papa: from a last vowel
        # 2.98 cm^2 off, the four-phone joint left every tube 4.2 off and a
        # time 76 ms off, -5.9 dB.)
        print("  joint over %d phones: worse than before it (%.4g against %.4g), not taken" % (
            phones, float(numpy.sum(reference ** 2)), initial_cost))
        parameters = initial
        reference = shape.residual(parameters, exact)
    found_tubes, found_arrivals, found_leavings = shape.split(parameters)
    clip = flow_out[start:start + size]
    print("  joint over %d phones to %.3f s: %.1f dB (against the sound's swing); tubes rms %s; times off %s ms" % (
        phones, end, 10.0 * numpy.log10(numpy.sum(reference ** 2) / numpy.sum((clip - clip.mean()) ** 2)),
        " ".join("%.3f" % tube_recover.tube_rms((found_tubes[k], length), (true_tubes[k], length)) for k in range(phones)),
        " ".join("%.2f" % (1e3 * value) for value in numpy.concatenate([
            found_leavings[:-1] - true_leavings[:phones - 1], found_arrivals[1:] - true_arrivals[1:phones]]))))
    sys.stdout.flush()
    arrivals = arrivals.copy()
    leavings = leavings.copy()
    arrivals[:phones] = found_arrivals
    leavings[:phones - 1] = found_leavings[:-1]
    if velum_areas is not None:
        found_velum = shape.velum(parameters)
        for k in nasal:
            velum_areas[k] = found_velum[k]
    return [found_tubes[k] for k in range(phones)], arrivals, leavings


JOINT_TIME_LIMIT = 0.003   # s: a joint that moves any time further is not taken
JOINT_TUBE_LIMIT = 1.0     # cm^2 rms: nor one that moves any phone's tube further
# (India's third run: from a vowel searched to -12 dB, the four-phone joint
# scored a little better than its start with every tube 23 cm^2 off and
# times 99 to 139 ms off. A joint refines; it does not relocate.)


def joint_sane(shape: Any, before: numpy.ndarray, after: numpy.ndarray) -> bool:
    """Whether a joint's result stays a refinement of its start."""
    first_tubes, first_arrivals, first_leavings = shape.split(before)
    last_tubes, last_arrivals, last_leavings = shape.split(after)
    moved = numpy.concatenate([last_arrivals - first_arrivals, last_leavings[:-1] - first_leavings[:-1]])
    if numpy.max(numpy.abs(moved)) > JOINT_TIME_LIMIT:
        return False
    for k in range(len(first_tubes)):
        if numpy.sqrt(numpy.mean((last_tubes[k] - first_tubes[k]) ** 2)) > JOINT_TUBE_LIMIT:
            return False
    return True


REUSE_DEPTH = -60.0     # dB: FAST_ALGO reuse accepts the earlier same vowel's node search this deep
                        # (the full search reached -59.8 and -79.0 dB for the second /aa/ of
                        # /aa b aa/ and /aa m aa/ in the reference)
EARLY_START = 0.005     # s: after a release, the next phone's first arrival this much before its hint
ARRIVAL_OFFSETS = numpy.arange(-0.012, 0.01201, 0.002)   # s about the hint (arrival_by_posture)
ARRIVAL_BUDGET = 20     # the posture fit's evaluations at each


def arrival_by_posture(pool: Any, hint: Dict[str, float], count: int, tubes: List[numpy.ndarray],
                       arrivals: numpy.ndarray, leavings: numpy.ndarray, length: float, index: int,
                       flow_out: numpy.ndarray, seal: float, port_areas: Any, jet: float, voiceless: Any,
                       fricatives: Any, leading: Dict[int, Any], true_arrivals: numpy.ndarray) -> numpy.ndarray:
    """A vowel's arrival after a release: at each candidate (the arrival as
    it stands less EARLY_START, i.e. the hint's, plus ARRIVAL_OFFSETS), a
    quick fit of its posture (the organs, from the hint) on its window's
    equation error in its sequence; the arrival whose fit goes deepest.
    (MEASURED, Papa's second /ah/ with the earlier phones true: the fit
    reached -12.7, -22.9, -24.3, -34.9, -45.0, -33.4, -43.7, -28.6 dB at
    -12, -8, -4, -2, 0, +2, +4, +8 ms from the true arrival; the hint's tube
    alone carries nothing of it, and the early start, 12.3 ms before the
    truth there, left even the true tube at -11.5 dB in its window.)"""
    shape = tube_recover.PostureShape(hint, tube_recover.WIDEN, count)
    center = arrivals[index] + EARLY_START
    best = (numpy.inf, arrivals[index])
    rows = []
    for offset in ARRIVAL_OFFSETS:
        trial = arrivals.copy()
        trial[index] = max(center + offset, leavings[index - 1] + RELEASE_OPEN + GAP)
        begin, span = window(trial[index] + SETTLE, leavings[index] - CLEAR)
        target = tube_recover.Target(flow_out, begin, span, 0.0, 1.0, True)
        if leading:
            target.context = PrefixContext([], trial, leavings, length, seal, port_areas, jet, voiceless, leading,
                                           fricatives)
        else:
            target.context = PrefixContext(tubes, trial, leavings, length, seal, port_areas, jet, voiceless,
                                           fricatives=fricatives)
        found, used = tube_recover.stage(pool, shape, numpy.zeros(len(posture_recover.FIELDS)), target,
                                         ARRIVAL_BUDGET)
        residual = shape.residual(found, target)
        clip = flow_out[begin:begin + span]
        decibels = 10.0 * numpy.log10(max(float(numpy.sum((residual * target.scale) ** 2)), 1e-300)
                                      / max(float(numpy.sum((clip - clip.mean()) ** 2)), 1e-300))
        rows.append("%+.0f:%.1f" % (1e3 * offset, decibels))
        if decibels < best[0]:
            best = (decibels, trial[index])
    out = arrivals.copy()
    out[index] = best[1]
    print("arrival by posture (ms about the hint: dB): %s; arrival %.5f (true %.5f)" % (
        " ".join(rows), out[index], true_arrivals[index]))
    sys.stdout.flush()
    return out


def refine_transition(pool, before, tube, arrivals, leavings, length, index, flow_out, true_arrivals,
                      true_leavings, velum_areas=None, polish=False):
    """Transition `index`'s leaving and arrival, coarse to fine, over a
    window around it, with the phones before as recovered and `tube` for
    the phone it reaches (and after); `polish`: `tube` is a found one (else
    a hint's)."""
    tubes = numpy.array([numpy.asarray(t, float) for t in before] + [numpy.asarray(tube, float)] * (len(arrivals) - len(before)))
    memo = GridMemo()
    begin, span = window(leavings[index - 1] - 0.02, arrivals[index] + 0.03)
    low = tube_recover.Target(flow_out, begin, span, 1500.0, 1.0, True)
    exact = tube_recover.Target(flow_out, begin, span, 0.0, 1.0, True)
    if numpy.any(tubes[index - 1] <= 2.0 * STOP_AREA) and numpy.any(tubes[index] <= 2.0 * STOP_AREA):
        # Seal to seal (India's /n/ -> /d/: one alveolar seal held while the
        # velum closes): no release; the times show only in the nasal
        # murmur's fading, inside the folds' window, which runs from the
        # first seal's closing to the second's release whatever these times
        # are: that window's sound against the sound, on a 2-D grid. (Scanned
        # as a release, the /n/'s leaving came out 17.6 ms early and the word
        # was lost: 5.0 dB.)
        # (MEASURED, true tubes: -16.4 at the truth, the arrival a funnel
        # 5.2 / 5.3 at +-1 ms to 7.2 at +-12, the leaving nearly flat beyond
        # +-1 ms (5.9 -> 6.4): the arrival alone first, then the leaving
        # alone, then 2-D to 10 us.)
        inside = tube_recover.Target(flow_out, begin, span, 0.0, 1.0, True)
        inside.hidden_only = True
        sweep = numpy.arange(-0.016, 0.01601, 0.001)
        with driver_ledger.stage("grid seal to seal"):
            arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, [0.0], sweep, inside,
                                             velum_areas, memo)
            arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, sweep, [0.0], inside,
                                             velum_areas, memo)
            for offsets in (numpy.arange(-0.001, 0.00101, 0.00025), numpy.arange(-0.00025, 0.000251, 0.00005),
                            numpy.arange(-0.00005, 0.0000501, 0.00001)):
                arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, offsets, offsets,
                                                 inside, velum_areas, memo)
        print("transition %d: leaving %.5f (true %.5f), arrival %.5f (true %.5f)" % (
            index, leavings[index - 1], true_leavings[index - 1], arrivals[index], true_arrivals[index]))
        sys.stdout.flush()
        return arrivals, leavings
    if numpy.any(tubes[index - 1] <= 2.0 * STOP_AREA):
        # Out of a seal: the release's time first, by the folds' window
        # alone (its opening's own motion makes it the release's alone: with
        # the next vowel at its hint, that minimum fell on the true release
        # to 0.05 ms where the whole residual's fell 0.5 ms late); then the
        # next phone's arrival by the whole.
        # The stretch: from 5 ms before the closure arrived (the plain peel
        # holds there) to 15 ms past the release as hinted, the same for
        # every candidate.
        # (A word's leading closure: from rest at 0, the forward run exact.)
        # Before the whole chain of sealed phones this release ends (India's
        # /d/ follows /n/ under one seal), and before the velum starts to
        # open if any of them is nasal (it opens over the transition into
        # the nasal, from the leaving before it): the plain peel holds only
        # there. (From 5 ms before the /d/'s own arrival, the velum still
        # open, the scan put the release 3 ms off and the word was lost.)
        chain = index - 1
        while chain - 1 >= 0 and numpy.any(tubes[chain - 1] <= 2.0 * STOP_AREA):
            chain -= 1
        moment = arrivals[chain] - 0.005
        for k in range(chain, index):
            if velum_areas is not None and velum_areas[k] > 0.0 and k > 0:
                moment = min(moment, leavings[k - 1] - 0.005)
        start = max(int(moment * STEP_RATE), 0)
        inside = tube_recover.Target(flow_out, start, int((leavings[index - 1] + 0.015) * STEP_RATE) - start, 0.0, 1.0,
                                     False)
        inside.output_from = start
        if velum_areas is not None and velum_areas[index - 1] > 0.0:
            # A nasal: the velum opens before the seal, so the plain peel
            # does not hold before it; its release's window (the march's)
            # keeps one length whatever the release's time.
            inside = tube_recover.Target(flow_out, begin, span, 0.0, 1.0, True)
            inside.hidden_only = True
        for rung, offsets in enumerate((numpy.arange(-0.01, 0.01001, 0.00025), numpy.arange(-0.0003, 0.000301, 0.000025))):
            with driver_ledger.stage("release scan %d" % rung):
                arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, offsets, [0.0], inside,
                                                 velum_areas, memo)
        if not polish:
            # A hint's tube says nothing of the next phone's arrival (/aa b
            # aa/: the cost flat at ~1e7 from -16 to +24 ms, against a clean
            # minimum with the true tube), and a late arrival spoils the
            # hold's window (+5 ms: -28 dB there, -5 ms: -66 dB): the hint's,
            # early, and the tube searched there.
            arrivals = arrivals.copy()
            arrivals[index] = max(arrivals[index] - EARLY_START, leavings[index - 1] + RELEASE_OPEN + GAP)
            rungs = ()
        else:
            # With a found tube: near the arrival as it stands (the coarse
            # rung jumped 18 ms on a tube 0.43 cm^2 off).
            rungs = (numpy.arange(-0.003, 0.0031, 0.001), numpy.arange(-0.0008, 0.00081, 0.0002),
                     numpy.arange(-0.00016, 0.000161, 0.00004))
        for rung, offsets in enumerate(rungs):
            with driver_ledger.stage("arrival grid %d" % rung):
                arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, [0.0], offsets, exact,
                                                 velum_areas, memo)
        print("transition %d: leaving %.5f (true %.5f), arrival %.5f (true %.5f)" % (
            index, leavings[index - 1], true_leavings[index - 1], arrivals[index], true_arrivals[index]))
        sys.stdout.flush()
        return arrivals, leavings
    if numpy.any(tubes[index] <= 2.0 * STOP_AREA) and not polish:
        # Into a seal: the arrival is a narrow canyon (Echo's /eh/ -> /k!/,
        # corrected model: ~1e-17 at the true pair, a plateau ~10^6.4-7 at
        # +-2 ms of arrival whatever the leaving, a clean V along the
        # leaving at the true arrival; with the leaving 8 ms off, a funnel
        # +-2.5 ms wide around the true arrival). The 2-D coarse grid (4 ms,
        # the hint 7.3 ms off) never sat in it and chose 0.0943 for 0.1000.
        # So the arrival alone first, then the leaving alone, then fine 2-D.
        sweep = numpy.arange(-0.012, 0.01201, 0.0005)
        with driver_ledger.stage("grid closing arrival"):
            arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, [0.0], sweep, low,
                                             velum_areas, memo)
        with driver_ledger.stage("grid closing leaving"):
            arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, sweep, [0.0], low,
                                             velum_areas, memo)
    else:
        coarse = numpy.arange(-0.016, 0.0161, 0.004)
        with driver_ledger.stage("grid coarse"):
            arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, coarse, coarse, low,
                                             velum_areas, memo)
    fine = numpy.arange(-0.003, 0.0031, 0.001)
    with driver_ledger.stage("grid fine"):
        arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, fine, fine, exact, velum_areas,
                                         memo)
    finest = numpy.arange(-0.0008, 0.00081, 0.0002)
    with driver_ledger.stage("grid finest"):
        arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, finest, finest, exact,
                                         velum_areas, memo)
    # A release's basin is ~0.1-0.25 ms wide (the voicing onset after it
    # amplifies the folds' state difference): one rung finer.
    finer = numpy.arange(-0.00016, 0.000161, 0.00004)
    with driver_ledger.stage("grid finer"):
        arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, finer, finer, exact,
                                         velum_areas, memo)
    print("transition %d: leaving %.4f (true %.4f), arrival %.4f (true %.4f)" % (
        index, leavings[index - 1], true_leavings[index - 1], arrivals[index], true_arrivals[index]))
    sys.stdout.flush()
    return arrivals, leavings


class Word:
    """A known-answer word: its phones' kinds, the true tubes and times, the
    true sound, and the hints the search starts from."""

    def __init__(self, specs: Sequence[str]) -> None:
        self.specs = list(specs)


def build_word(specs: Sequence[str]) -> Word:
    """The word's known answer from its phone specs (see main)."""
    word = Word(specs)
    profile = table(os.path.join(DATA, "profile_bdl5.organs"))
    hints = table(os.path.join(DATA, "organs_vtl_tongue.txt"))
    # A stop is written NAME@place: its tube the phone before's, one node
    # sealed (the lips: the node before the open end). A nasal adds ^area:
    # the velum's port open to that area while it is held. A fricative is
    # NAME@place~area: that node narrowed to the area (not sealed), the
    # reversible jet gated on through it; NAME@place=area a narrowing without
    # it (an /l/). A trailing ! makes it voiceless (the folds held
    # apart through it: fold_keyframes).
    specs = list(specs)
    names = [spec.split("@")[0] for spec in specs]
    phones = len(names)
    # The count and length from the first vowel (a word may lead with a consonant).
    first, length = posture_recover.nodes(profile[[spec for spec in specs if "@" not in spec][0]])
    count = len(first)
    # Seal nodes by place, from the open end (46 nodes over 17.0 cm, 0.37 cm
    # a node): the lips the node before it; the alveolar ridge ~1.8 cm
    # behind them; the velum ~6.5 cm behind them.
    places = {"lips": count - 2, "alveolar": count - 6, "velar": count - 18}
    stops = {}
    narrowings = {}
    fricatives = []
    voiceless = {}
    true_velum = numpy.zeros(phones)
    for index, spec in enumerate(specs):
        if "@" in spec:
            place = spec.split("@")[1]
            unvoiced = place.endswith("!")
            place = place.rstrip("!")
            if "^" in place:
                place, area = place.split("^")
                true_velum[index] = float(area)
            if "~" in place:
                place, area = place.split("~")
                narrowings[index] = float(area)
                fricatives.append(index)
            if "=" in place:
                place, area = place.split("=")
                narrowings[index] = float(area)
            stops[index] = places[place]
            if unvoiced:
                voiceless[index] = FRICATIVE_FOLDS if index in narrowings else PLOSIVE_FOLDS
    nasals = [k for k in range(phones) if true_velum[k] > 0.0]
    seal = SEAL if len(stops) > len(narrowings) else 0.0
    jet = JET_GAIN if fricatives else 0.0
    fricatives = tuple(fricatives)
    # The vowels first: a leading consonant is made from the vowel after it.
    true_tubes = numpy.zeros((phones, count))
    for index, name in enumerate(names):
        if index not in stops:
            true_tubes[index] = posture_recover.nodes(profile[name], count)[0]
    for index in sorted(stops):
        true_tubes[index] = true_tubes[vowel_before(index, stops)]
        true_tubes[index, stops[index]] = narrowings.get(index, STOP_AREA)
    true_arrivals, true_leavings = plan(phones)
    seconds = true_leavings[-1] + 0.05
    times, frames = keyframes(true_tubes, true_arrivals, true_leavings)
    true_port = velum_keyframes(true_arrivals, true_leavings, true_velum) if nasals else None
    flow_in, flow_out = posture_recover.render((true_tubes[0], length), seconds, keyframes=(times, frames), seal=seal,
                                               velum=true_port, jet=jet,
                                               jet_gate=jet_keyframes(true_arrivals, true_leavings, fricatives),
                                               folds=fold_keyframes(true_arrivals, true_leavings, voiceless))
    word.seal, word.jet, word.voiceless, word.fricatives = seal, jet, voiceless, fricatives
    configure_targets(word)
    start, size = window(0.02, seconds - 0.01)
    whole = tube_recover.Target(flow_out, start, size, 0.0, 1.0, True)
    reference = flow_in[start:start + size]
    shape = SequenceShape(phones, count, length, stops, nasals, narrowings)
    truth = shape.join(true_tubes, true_arrivals, true_leavings, true_velum)
    # The hinted times: the truth's, each transition's leaving and arrival
    # pushed HINT_OFFSET apart or together; a nasal's velum hinted at NASAL_HINT.
    arrivals = true_arrivals.copy()
    leavings = true_leavings.copy()
    for index in range(1, phones):
        sign = 1.0 if index % 2 else -1.0
        leavings[index - 1] += sign * HINT_OFFSET
        arrivals[index] -= sign * HINT_OFFSET
    velum = numpy.zeros(phones)
    for k in nasals:
        velum[k] = NASAL_HINT
    port_areas = velum if nasals else None
    word.truth = truth
    word.profile = profile
    word.hints = hints
    word.names = names
    word.phones = phones
    word.count = count
    word.length = length
    word.stops = stops
    word.narrowings = narrowings
    word.fricatives = fricatives
    word.voiceless = voiceless
    word.true_velum = true_velum
    word.nasals = nasals
    word.seal = seal
    word.jet = jet
    word.true_tubes = true_tubes
    word.true_arrivals = true_arrivals
    word.true_leavings = true_leavings
    word.true_port = true_port
    word.seconds = seconds
    word.flow_in = flow_in
    word.flow_out = flow_out
    word.start = start
    word.size = size
    word.whole = whole
    word.reference = reference
    word.shape = shape
    word.hinted_arrivals = arrivals
    word.hinted_leavings = leavings
    word.hinted_velum = velum
    return word


def configure_targets(word: Word) -> None:
    """The word's kinds as tube_recover's target defaults (every Target made
    after this carries them)."""
    tube_recover.DEFAULT_SEAL = word.seal
    tube_recover.DEFAULT_JET = word.jet
    tube_recover.DEFAULT_VOICELESS = word.voiceless
    tube_recover.DEFAULT_FRICATIVES = word.fricatives


# The first vowel's search (no context: the sound's start, the tube alone)
# depends only on its hint, its window, the node count, the sound up to the
# window's end plus the peel's reach, and the engine: the same in every word
# that starts with that vowel and that leaving (/aa b aa/ and /aa m aa/ ran
# it identically, line for line). Kept on disk under the engine's hash.
FIRST_VOWEL_CACHE = os.path.join(posture_recover.ROOT, "out", "cache")
FIRST_VOWEL_REACH = 4   # node counts of steps past the window the key holds (column_locality.LOOKAHEAD_NODES)


def engine_digest() -> str:
    """The engine's bytes (the tool, the probe, and with RENDERER=engine the
    in-process library)."""
    digest = hashlib.sha256()
    paths = [posture_recover.SUSTAIN, posture_recover.PROBE]
    digest.update(posture_recover.RENDERER.encode())
    if posture_recover.RENDERER == "engine":
        import engine
        paths.append(engine.LIBRARY)
    for path in paths:
        with open(path, "rb") as handle:
            digest.update(handle.read())
    return digest.hexdigest()[:16]


def first_vowel_key(hint: Dict[str, float], flow_out: numpy.ndarray, begin: int, span: int, count: int) -> str:
    digest = hashlib.sha256()
    digest.update(repr(sorted(hint.items())).encode())
    digest.update(repr((begin, span, count, tube_recover.SURVEY, tube_recover.WIDEN, tube_recover.KEEP,
                        tube_recover.FINISH, tube_recover.BENDS, tube_recover.EXACT, tube_recover.NODE_RUNGS,
                        sorted(tube_recover.ALGO), posture_recover.VISCOUS, posture_recover.HOP,
                        posture_recover.PITCH_EAR, getattr(reverse_core, "ENGINE_GLIDE_PHASE", None))).encode())
    reach = min(begin + span + FIRST_VOWEL_REACH * count, len(flow_out))
    digest.update(numpy.ascontiguousarray(flow_out[:reach], float).tobytes())
    digest.update(engine_digest().encode())
    return digest.hexdigest()[:24]


def recover_first_vowel(pool: Any, hint: Dict[str, float], flow_out: numpy.ndarray, begin: int, span: int,
                        truth: Any, count: int) -> Tuple[float, Any]:
    """tube_recover.recover for a word's first vowel without context, from
    the cache when this exact search has been run."""
    if not tube_recover.EXACT_SHORTCUTS or os.environ.get("FIRST_VOWEL_CACHE", "1") != "1":
        return tube_recover.recover(pool, hint, flow_out, begin, span, truth, count, None)
    path = os.path.join(FIRST_VOWEL_CACHE, "first_vowel_%s.npz" % first_vowel_key(hint, flow_out, begin, span, count))
    if os.path.exists(path):
        kept = numpy.load(path)
        print("(first vowel from the cache: %s)" % os.path.basename(path))
        return float(kept["value"]), (numpy.array(kept["areas"]), float(kept["length"]))
    value, found = tube_recover.recover(pool, hint, flow_out, begin, span, truth, count, None)
    os.makedirs(FIRST_VOWEL_CACHE, exist_ok=True)
    numpy.savez(path, value=value, areas=found[0], length=found[1])
    return value, found


def recover_word(pool: Any, word: Word) -> Tuple[numpy.ndarray, Any]:
    """Phone by phone, then the joint: the recovered parameters (word.shape)."""
    names = word.names
    hints = word.hints
    count = word.count
    length = word.length
    stops = word.stops
    narrowings = word.narrowings
    nasals = word.nasals
    seal = word.seal
    jet = word.jet
    voiceless = word.voiceless
    fricatives = word.fricatives
    flow_out = word.flow_out
    true_tubes = word.true_tubes
    true_arrivals = word.true_arrivals
    true_leavings = word.true_leavings
    true_velum = word.true_velum
    phones = word.phones
    shape = word.shape
    whole = word.whole
    reference = word.reference
    arrivals = word.hinted_arrivals.copy()
    leavings = word.hinted_leavings.copy()
    velum = word.hinted_velum.copy()
    port_areas = velum if nasals else None
    tubes = []
    # The word's leading consonants: made from the first vowel's tube,
    # searched with it (each candidate makes them).
    leading = {}
    for index in range(phones):
        if index not in stops:
            break
        leading[index] = (stops[index], narrowings.get(index, STOP_AREA))
    first_vowel = len(leading)
    for index, name in enumerate(names):
        with driver_ledger.stage("phone %d %s" % (index, name)):
            if index in leading:
                # Its tube for now from the first vowel's hint.
                tube = posture_recover.nodes(hints[names[first_vowel]], count)[0]
                tube[stops[index]] = narrowings.get(index, STOP_AREA)
                tubes.append(tube)
                continue
            if index in stops:
                # A stop (or a nasal): its tube is known from the phone before;
                # the closure's times are searched (its hold is unheard by the
                # peel), and a nasal's velum by the window's forward sound.
                guess = numpy.array(tubes[vowel_before(index, stops)])
                guess[stops[index]] = narrowings.get(index, STOP_AREA)
                with driver_ledger.stage("transition"):
                    arrivals, leavings = refine_transition(pool, tubes, guess, arrivals, leavings, length, index,
                                                           flow_out, true_arrivals, true_leavings, port_areas)
                if index in nasals:
                    with driver_ledger.stage("velum grid"):
                        velum[index] = velum_grid(pool, tubes, guess, arrivals, leavings, velum, length, index,
                                                  flow_out)
                    print("%s: velum %.3f (true %.3f)" % (name, velum[index], true_velum[index]))
                tubes.append(guess)
                print("%s: %s at node %d" % (name, "narrowed" if index in narrowings else "closed, sealed",
                                             stops[index]))
                sys.stdout.flush()
                continue
            if index > 0:
                guess = posture_recover.nodes(hints[name], count)[0]
                with driver_ledger.stage("transition"):
                    arrivals, leavings = refine_transition(pool, tubes, guess, arrivals, leavings, length, index,
                                                           flow_out, true_arrivals, true_leavings, port_areas)
                if numpy.any(numpy.asarray(tubes[index - 1]) <= 2.0 * STOP_AREA):
                    with driver_ledger.stage("arrival by posture"):
                        arrivals = arrival_by_posture(pool, hints[name], count, tubes, arrivals, leavings, length,
                                                      index, flow_out, seal, port_areas, jet, voiceless, fricatives,
                                                      leading if index == first_vowel else {}, true_arrivals)
            begin, span = window(arrivals[index] + SETTLE, leavings[index] - CLEAR)
            context = None
            if index == first_vowel and leading:
                context = PrefixContext([], arrivals, leavings, length, seal, port_areas, jet, voiceless, leading,
                                        fricatives)
            elif index > 0 or jet > 0.0:
                context = PrefixContext(tubes, arrivals, leavings, length, seal, port_areas, jet, voiceless,
                                        fricatives=fricatives)
            reused = None
            earlier = [k for k in range(index) if names[k] == name and k not in stops]
            if "reuse" in tube_recover.ALGO and earlier and context is not None:
                # The same vowel earlier in the word: its tube as a node
                # search's start; deep enough, the survey is skipped.
                with driver_ledger.stage("reuse"):
                    reused = tube_recover.node_search(pool, (numpy.array(tubes[earlier[-1]]), length), flow_out,
                                                      begin, span, count, context)
                print("%s: from the earlier %s (phone %d): %.1f dB, tube rms %.3f" % (
                    name, name, earlier[-1], reused[0], tube_recover.tube_rms(reused[1], (true_tubes[index], length))))
                sys.stdout.flush()
            if reused is not None and reused[0] < REUSE_DEPTH:
                value, found = reused
            else:
                with driver_ledger.stage("recover"):
                    if context is None:
                        value, found = recover_first_vowel(pool, hints[name], flow_out, begin, span,
                                                           (true_tubes[index], length), count)
                    else:
                        value, found = tube_recover.recover(pool, hints[name], flow_out, begin, span,
                                                            (true_tubes[index], length), count, context)
                if reused is not None and reused[0] < value:
                    value, found = reused
            print("%s: %.1f dB, tube rms %.3f" % (name, value, tube_recover.tube_rms(found, (true_tubes[index], length))))
            sys.stdout.flush()
            for repeat in range(POLISH if index > 0 else 0):
                # The times again with this tube, then the tube again at them
                # (/iy/: searched with the transition 3 ms off, -42.9 dB and
                # 0.275 cm² off; the times then came within 0.2 ms).
                if index == first_vowel:
                    for k, (node, area) in leading.items():
                        tubes[k] = numpy.array(found[0])
                        tubes[k][node] = area
                with driver_ledger.stage("polish %d transition" % repeat):
                    arrivals, leavings = refine_transition(pool, tubes, found[0], arrivals, leavings, length,
                                                           index, flow_out, true_arrivals, true_leavings,
                                                           port_areas, True)
                begin, span = window(arrivals[index] + SETTLE, leavings[index] - CLEAR)
                if index == first_vowel and leading:
                    polish_context = PrefixContext([], arrivals, leavings, length, seal, port_areas, jet, voiceless,
                                                   leading, fricatives)
                else:
                    polish_context = PrefixContext(tubes, arrivals, leavings, length, seal, port_areas, jet, voiceless,
                                                   fricatives=fricatives)
                with driver_ledger.stage("polish %d node search" % repeat):
                    value, found = tube_recover.node_search(pool, found, flow_out, begin, span, count,
                                                            polish_context)
                print("%s polished: %.1f dB, tube rms %.3f" % (
                    name, value, tube_recover.tube_rms(found, (true_tubes[index], length))))
                sys.stdout.flush()
            if index == first_vowel:
                for k, (node, area) in leading.items():
                    tubes[k] = numpy.array(found[0])
                    tubes[k][node] = area
            tubes.append(found[0])
            if index > 0:
                with driver_ledger.stage("prefix joint"):
                    tubes, arrivals, leavings = prefix_joint(pool, tubes, arrivals, leavings, length, flow_out,
                                                             leavings[index] - CLEAR, true_tubes, true_arrivals,
                                                             true_leavings, stops, port_areas, narrowings)
    tubes = numpy.array(tubes)
    parameters = shape.join(tubes, arrivals, leavings, velum)
    print("before the joint search: %.1f dB" % level(shape.residual(parameters, whole), reference))
    sys.stdout.flush()
    rounds = StillRounds()
    more = "jointmore" in tube_recover.ALGO
    for rung in range(JOINT_ROUNDS_MAX if more else 4):
        if more and rung >= 4 and (rounds.gain < JOINT_MORE_DB or rounds.stopped is not None):
            print("joint: ended after %d rounds (the last gained %.2f dB)" % (rung, rounds.gain))
            break
        start_point = parameters.copy()
        with driver_ledger.stage("final joint %d" % rung):
            parameters, used, skipped = rounds.run(pool, shape, parameters, whole, "whole", 30)
        if not joint_sane(shape, start_point, parameters):
            print("joint %d: relocated rather than refined (a time over %.0f ms or a tube over %.1f cm^2), not taken; "
                  "rounds ended" % (rung, 1e3 * JOINT_TIME_LIMIT, JOINT_TUBE_LIMIT))
            parameters = start_point
            break
        found_tubes, found_arrivals, found_leavings = shape.split(parameters)
        if skipped:
            print("joint %d: skipped (%s)" % (rung, rounds.reason))
            sys.stdout.flush()
            continue
        print("joint %d: %d evaluations; %.1f dB; tubes rms %s; times off %s ms; velum %s" % (
            rung, used, level(shape.residual(parameters, whole), reference),
            " ".join("%.3f" % tube_recover.tube_rms((found_tubes[k], length), (true_tubes[k], length))
                     for k in range(phones)),
            " ".join("%.2f" % (1e3 * value) for value in numpy.concatenate([
                found_leavings - true_leavings, found_arrivals[1:] - true_arrivals[1:]])),
            " ".join("%.3f" % v for v in shape.velum(parameters)[nasals])))
        sys.stdout.flush()
    if "blocks" in tube_recover.ALGO:
        parameters = block_rounds(pool, shape, parameters, whole, reference,
                                  BlockReport(shape, length, true_tubes, true_arrivals, true_leavings, nasals))
    return parameters, velum


class BlockReport:
    """block_rounds' line per round: the error, the tubes and times off."""

    def __init__(self, shape: Any, length: float, true_tubes: numpy.ndarray, true_arrivals: numpy.ndarray,
                 true_leavings: numpy.ndarray, nasals: List[int]) -> None:
        self.shape = shape
        self.length = length
        self.true_tubes = true_tubes
        self.true_arrivals = true_arrivals
        self.true_leavings = true_leavings
        self.nasals = nasals

    def __call__(self, label: str, parameters: numpy.ndarray, decibels: float) -> None:
        found_tubes, found_arrivals, found_leavings = self.shape.split(parameters)
        print("%s: %.1f dB; tubes rms %s; times off %s us; velum %s" % (
            label, decibels,
            " ".join("%.4f" % tube_recover.tube_rms((found_tubes[k], self.length), (self.true_tubes[k], self.length))
                     for k in range(len(found_tubes))),
            " ".join("%.2f" % (1e6 * value) for value in numpy.concatenate([
                found_leavings[:-1] - self.true_leavings[:-1], found_arrivals[1:] - self.true_arrivals[1:]])),
            " ".join("%.4f" % v for v in self.shape.velum(parameters)[self.nasals])))
        sys.stdout.flush()


def main():
    word = build_word(sys.argv[1].split(","))
    configure_targets(word)
    specs = word.specs
    length = word.length
    seal = word.seal
    jet = word.jet
    voiceless = word.voiceless
    fricatives = word.fricatives
    nasals = word.nasals
    true_tubes = word.true_tubes
    true_arrivals = word.true_arrivals
    true_leavings = word.true_leavings
    true_velum = word.true_velum
    true_port = word.true_port
    seconds = word.seconds
    flow_out = word.flow_out
    start = word.start
    size = word.size
    shape = word.shape
    print("pool %d workers (POOL), chunk %d; exact shortcuts %s (resume %s); algorithmic changes %s; renderer %s; "
          "model: pitch ear %s, glide phase %s" % (
              driver_ledger.POOL, driver_ledger.CHUNK, "on" if tube_recover.EXACT_SHORTCUTS else "off",
              "on" if RESUME else "off", ",".join(sorted(tube_recover.ALGO)) or "none", posture_recover.RENDERER,
              posture_recover.PITCH_EAR_SETTING, getattr(reverse_core, "ENGINE_GLIDE_PHASE", "-")))
    print("%s: %d nodes, %.2f cm, %.2f s; the truth's error %.1f dB" % (
        " ".join(specs), word.count, length, seconds, level(shape.residual(word.truth, word.whole), word.reference)))
    sys.stdout.flush()
    with Pool(driver_ledger.POOL) as pool:
        with driver_ledger.stage("word"):
            parameters, velum = recover_word(pool, word)
    found_tubes, found_arrivals, found_leavings = shape.split(parameters)
    found_velum = shape.velum(parameters)
    times, frames = keyframes(found_tubes, found_arrivals, found_leavings)
    found_port = velum_keyframes(found_arrivals, found_leavings, found_velum) if nasals else None
    _, again = posture_recover.render((found_tubes[0], length), seconds, keyframes=(times, frames), seal=seal,
                                      velum=found_port, jet=jet,
                                      jet_gate=jet_keyframes(found_arrivals, found_leavings, fricatives),
                                      folds=fold_keyframes(found_arrivals, found_leavings, voiceless))
    print("round trip: the recovered sequence's sound against the true one's: %.1f dB" % reverse_core.snr(
        flow_out[start:start + size], again[start:start + size]))
    report_errors(word, found_tubes, found_arrivals, found_leavings, found_velum)
    report_ledger(word)
    # Kept for listening and for later stages: the recovered keyframes, and
    # both sounds as the engine's own wave files (sustain's 48 kHz output).
    label = "_".join(spec.replace("@", "-").replace("^", "v").replace("~", "n").replace("!", "x")
                     for spec in specs).lower()
    os.makedirs(OUT, exist_ok=True)
    numpy.savez(os.path.join(OUT, "sequence_%s.npz" % label), tubes=found_tubes, arrivals=found_arrivals,
                leavings=found_leavings, velum=found_velum, length=length, true_tubes=true_tubes,
                true_arrivals=true_arrivals, true_leavings=true_leavings, true_velum=true_velum)
    true_times, true_frames = keyframes(true_tubes, true_arrivals, true_leavings)
    posture_recover.render((true_tubes[0], length), seconds, keyframes=(true_times, true_frames),
                           wave=os.path.join(OUT, "sequence_%s_original.wav" % label), seal=seal, velum=true_port,
                           jet=jet, folds=fold_keyframes(true_arrivals, true_leavings, voiceless),
                           jet_gate=jet_keyframes(true_arrivals, true_leavings, fricatives))
    posture_recover.render((found_tubes[0], length), seconds, keyframes=(times, frames),
                           wave=os.path.join(OUT, "sequence_%s_recovered.wav" % label), seal=seal, velum=found_port,
                           jet=jet, folds=fold_keyframes(found_arrivals, found_leavings, voiceless),
                           jet_gate=jet_keyframes(found_arrivals, found_leavings, fricatives))


def report_errors(word: Word, found_tubes: numpy.ndarray, found_arrivals: numpy.ndarray,
                  found_leavings: numpy.ndarray, found_velum: numpy.ndarray) -> None:
    """The recovered answer against the truth, in one line."""
    rms = [tube_recover.tube_rms((found_tubes[k], word.length), (word.true_tubes[k], word.length))
           for k in range(word.phones)]
    times = numpy.concatenate([found_leavings[:-1] - word.true_leavings[:-1], found_arrivals[1:] - word.true_arrivals[1:]])
    print("errors: tubes rms max %.4f cm^2 (%s); times max %.4f ms; velum %s" % (
        max(rms), " ".join("%.4f" % value for value in rms), 1e3 * float(numpy.max(numpy.abs(times))),
        " ".join("%.4f (true %.3f)" % (found_velum[k], word.true_velum[k]) for k in word.nasals) or "-"))
    sys.stdout.flush()


def report_ledger(word: Word) -> None:
    """Where the time went (driver_ledger), the columns' projected
    locality, and the ledger kept beside the word's output."""
    entry = driver_ledger.LEDGER.records.get("word")
    if entry is not None:
        print("total: %d evaluations (%d saved by exact shortcuts), wall %.1f min, CPU %.1f min (pool %d)" % (
            entry.evaluations, entry.saved, entry.wall / 60.0, (entry.worker_cpu + entry.main_cpu) / 60.0,
            driver_ledger.POOL))
    driver_ledger.LEDGER.report(depth=int(os.environ.get("LEDGER_DEPTH", "4")))
    driver_ledger.report_kinds()
    column_locality.report(sys.stdout)
    label = "_".join(spec.replace("@", "-").replace("^", "v").replace("~", "n").replace("!", "x")
                     for spec in word.specs).lower()
    os.makedirs(OUT, exist_ok=True)
    driver_ledger.LEDGER.dump(os.path.join(OUT, "ledger_%s%s.json" % (label, os.environ.get("LEDGER_TAG", ""))))


def velum_cost(task):
    tubes, arrivals, leavings, velum_areas, length, target = task
    return float(numpy.sum(word_residual(tubes, arrivals, leavings, velum_areas, length, target) ** 2))


def velum_grid(pool, before, tube, arrivals, leavings, velum_areas, length, index, flow_out):
    """A nasal phone's velum area: the best of a grid by the equation error
    over its hold (sealed, the march peels the nose: /aa m aa/ with the
    vowel after at its hint, an exact minimum at the true 0.50; through the
    release and on, where the phones after are still guesses, the cost was
    flat from 0.45 to 0.55 and the search's grid chose 0.84), coarse to
    fine."""
    tubes = numpy.array([numpy.asarray(t, float) for t in before] + [numpy.asarray(tube, float)] * (len(arrivals) - len(before)))
    begin, span = window(arrivals[index] + 0.005, leavings[index] - 0.002)
    target = tube_recover.Target(flow_out, begin, span, 0.0, 1.0, True)
    best = velum_areas[index]
    for rung in range(3):
        if rung == 0:
            grid = numpy.arange(0.05, 1.01, 0.05)
        elif rung == 1:
            grid = best + numpy.arange(-0.04, 0.041, 0.005)
        else:
            grid = best + numpy.arange(-0.005, 0.0051, 0.0005)
        tasks = []
        for area in grid:
            trial = numpy.array(velum_areas, float)
            trial[index] = max(float(area), 0.005)
            tasks.append((tubes, arrivals, leavings, trial, length, target))
        with driver_ledger.stage("velum rung %d" % rung):
            costs = numpy.array(driver_ledger.pool_map(pool, velum_cost, tasks))
        costs[~numpy.isfinite(costs)] = numpy.inf
        best = tasks[int(numpy.argmin(costs))][3][index]
    return best


if __name__ == "__main__":
    main()
