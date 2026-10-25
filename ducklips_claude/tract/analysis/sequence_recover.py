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
import os
import sys
from multiprocessing import Pool

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import posture_recover
import reverse_core
import tube_recover
import observer
import marching

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

    def __call__(self, forced, observed=None):
        flow, _ = posture_recover.render((self.first, self.length), self.seconds, keyframes=(self.times, self.frames),
                                         forced=forced, seal=self.seal, velum=self.velum, folds=self.folds,
                                         jet=self.jet, jet_gate=self.jet_gate,
                                         observed=None if observed is None else (observed, OBSERVE_GAIN))
        if len(flow) < len(forced):
            flow = numpy.concatenate([flow, numpy.zeros(len(forced) - len(flow))])
        return flow[:len(forced)]


NOSE = None   # the reversible nose's constants (reverse_core.read_nose), read once


def nose_constants():
    global NOSE
    if NOSE is None:
        path = os.path.join("/tmp/es", "nose_constants_%d.txt" % os.getpid())
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
    count = len(target.flow_out)
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

    def returning(arrays):
        flow, pressure = marching._peel_below(geometry, port_node, arrays["port_pressure"], arrays["from_front"],
                                              arrays["returned"], port_y, posture_recover.HOP, True, jet_steps)
        return flow, pressure, 0.5 * (pressure - flow * reverse_core.RHO_C / area0)

    def render_waves(forced):
        path = os.path.join("/tmp/es", "waves_%d_%d.raw" % (os.getpid(), numpy.random.randint(1 << 30)))
        flow, _ = posture_recover.render((frames[0], length), seconds, keyframes=(times, frames), forced=forced,
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

    def render_window(arrays, begin):
        flow, pressure, wave = returning(arrays)
        forced = wave.copy()
        forced[max(begin - port_node - 2, 0):] = numpy.nan
        return render_waves(forced)[1]

    share, returned, port_pressure, from_front, mode, waves, in_window = marching.march_windows(
        target.flow_out, geometry, port, nose, posture_recover.HOP, end_extra, target.seal, windows, render_window,
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


def sealed_residual(tubes, arrivals, leavings, length, target):
    """The equation error with closures (observer.implied): zero where the
    sound cannot reveal the inside and the folds run free."""
    times, frames = keyframes(tubes, arrivals, leavings)
    schedule = reverse_core.keyframe_schedule(times, frames, len(target.flow_out))
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    # (The free renders span the whole sound: a closure past the window's end
    # still needs our folds' flow.)
    free = FreeFlow(frames[0], length, (len(target.flow_out) + 1) / STEP_RATE, times, frames, target.seal,
                    folds=fold_keyframes(arrivals, leavings, target.voiceless), jet=target.jet,
                    jet_gate=jet_keyframes(arrivals, leavings, target.fricatives))
    implied, pressure, revealed, modelled = observer.implied(geometry, target.flow_out, target.seal,
                                                             posture_recover.HOP,
                                                             length - frames.shape[1] * posture_recover.NODE, free,
                                                             hide=HIDE_AREA, jet=target.jet,
                                                             jet_closing=posture_recover.JET_CLOSING,
                                                             jet_gate=jet_gate_for(arrivals, leavings,
                                                                                   target.fricatives,
                                                                                   len(target.flow_out)))
    forced = 0.5 * (pressure - implied * reverse_core.RHO_C / geometry[:, 0, 0])
    forced[~revealed] = numpy.nan
    observed = None
    if OBSERVE_GAIN > 0.0:
        observed = implied.copy()
        observed[~revealed] = numpy.nan
    driven = free(forced, observed)
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
    return tube_recover.fixed_clip(difference, target)


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
    reach = min(max(begin, 1) + 4 * frames.shape[1], len(target.flow_out))
    schedule = reverse_core.keyframe_schedule(times, frames, reach)
    geometry = reverse_core.moving_geometry(schedule, posture_recover.VISCOUS)
    implied, pressure = reverse_core.peel_moving_glottis(geometry, target.flow_out[:reach], hop=posture_recover.HOP,
                                                         end_extra=length - frames.shape[1] * posture_recover.NODE,
                                                         jet=target.jet,
                                                         gate=jet_gate_for(arrivals, leavings, target.fricatives,
                                                                           len(target.flow_out)))
    forced = numpy.full(len(target.flow_out), numpy.nan)
    forced[:begin] = (0.5 * (pressure - implied * reverse_core.RHO_C / geometry[:, 0, 0]))[:begin]
    velum = None
    if velum_areas is not None and numpy.any(numpy.asarray(velum_areas) > 0.0):
        velum = velum_keyframes(arrivals, leavings, velum_areas)
    _, sound = posture_recover.render((frames[0], length), (len(target.flow_out) + 1) / STEP_RATE,
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


def timing_grid(pool, tubes, arrivals, leavings, length, index, leaving_offsets, arrival_offsets, target,
                velum_areas=None):
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
    costs = numpy.array(pool.map(timing_cost, tasks))
    costs[~numpy.isfinite(costs)] = numpy.inf
    best = tasks[int(numpy.argmin(costs))]
    return best[1], best[2]


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
    for bend in JOINT_BENDS:
        target = tube_recover.Target(flow_out, start, size, 0.0, 1.0, True, bend)
        parameters, used = tube_recover.stage(pool, shape, parameters, target, 30)
    found_tubes, found_arrivals, found_leavings = shape.split(parameters)
    exact = tube_recover.Target(flow_out, start, size, 0.0, 1.0, True)
    reference = shape.residual(parameters, exact)
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


EARLY_START = 0.005     # s: after a release, the next phone's first arrival this much before its hint


def refine_transition(pool, before, tube, arrivals, leavings, length, index, flow_out, true_arrivals,
                      true_leavings, velum_areas=None, polish=False):
    """Transition `index`'s leaving and arrival, coarse to fine, over a
    window around it, with the phones before as recovered and `tube` for
    the phone it reaches (and after); `polish`: `tube` is a found one (else
    a hint's)."""
    tubes = numpy.array([numpy.asarray(t, float) for t in before] + [numpy.asarray(tube, float)] * (len(arrivals) - len(before)))
    begin, span = window(leavings[index - 1] - 0.02, arrivals[index] + 0.03)
    low = tube_recover.Target(flow_out, begin, span, 1500.0, 1.0, True)
    exact = tube_recover.Target(flow_out, begin, span, 0.0, 1.0, True)
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
        start = max(int((arrivals[index - 1] - 0.005) * STEP_RATE), 0)
        inside = tube_recover.Target(flow_out, start, int((leavings[index - 1] + 0.015) * STEP_RATE) - start, 0.0, 1.0,
                                     False)
        inside.output_from = start
        if velum_areas is not None and velum_areas[index - 1] > 0.0:
            # A nasal: the velum opens before the seal, so the plain peel
            # does not hold before it; its release's window (the march's)
            # keeps one length whatever the release's time.
            inside = tube_recover.Target(flow_out, begin, span, 0.0, 1.0, True)
            inside.hidden_only = True
        for offsets in (numpy.arange(-0.01, 0.01001, 0.00025), numpy.arange(-0.0003, 0.000301, 0.000025)):
            arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, offsets, [0.0], inside,
                                             velum_areas)
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
        for offsets in rungs:
            arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, [0.0], offsets, exact,
                                             velum_areas)
        print("transition %d: leaving %.5f (true %.5f), arrival %.5f (true %.5f)" % (
            index, leavings[index - 1], true_leavings[index - 1], arrivals[index], true_arrivals[index]))
        sys.stdout.flush()
        return arrivals, leavings
    coarse = numpy.arange(-0.016, 0.0161, 0.004)
    arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, coarse, coarse, low, velum_areas)
    fine = numpy.arange(-0.003, 0.0031, 0.001)
    arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, fine, fine, exact, velum_areas)
    finest = numpy.arange(-0.0008, 0.00081, 0.0002)
    arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, finest, finest, exact, velum_areas)
    # A release's basin is ~0.1-0.25 ms wide (the voicing onset after it
    # amplifies the folds' state difference): one rung finer.
    finer = numpy.arange(-0.00016, 0.000161, 0.00004)
    arrivals, leavings = timing_grid(pool, tubes, arrivals, leavings, length, index, finer, finer, exact, velum_areas)
    print("transition %d: leaving %.4f (true %.4f), arrival %.4f (true %.4f)" % (
        index, leavings[index - 1], true_leavings[index - 1], arrivals[index], true_arrivals[index]))
    sys.stdout.flush()
    return arrivals, leavings


def main():
    profile = table(os.path.join(DATA, "profile_bdl5.organs"))
    hints = table(os.path.join(DATA, "organs_vtl_tongue.txt"))
    # A stop is written NAME@place: its tube the phone before's, one node
    # sealed (the lips: the node before the open end). A nasal adds ^area:
    # the velum's port open to that area while it is held. A fricative is
    # NAME@place~area: that node narrowed to the area (not sealed), the
    # reversible jet gated on through it; NAME@place=area a narrowing without
    # it (an /l/). A trailing ! makes it voiceless (the folds held
    # apart through it: fold_keyframes).
    specs = sys.argv[1].split(",")
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
    tube_recover.DEFAULT_SEAL = seal
    tube_recover.DEFAULT_JET = jet
    tube_recover.DEFAULT_VOICELESS = voiceless
    tube_recover.DEFAULT_FRICATIVES = fricatives
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
    start, size = window(0.02, seconds - 0.01)
    whole = tube_recover.Target(flow_out, start, size, 0.0, 1.0, True)
    reference = flow_in[start:start + size]
    shape = SequenceShape(phones, count, length, stops, nasals, narrowings)
    truth = shape.join(true_tubes, true_arrivals, true_leavings, true_velum)
    print("%s: %d nodes, %.2f cm, %.2f s; the truth's error %.1f dB" % (
        " ".join(specs), count, length, seconds, level(shape.residual(truth, whole), reference)))
    sys.stdout.flush()
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
    with Pool(10) as pool:
        # Phone by phone. Our folds are an oscillator whose phase is neutral:
        # an earlier error (the hinted transition 8 ms off) leaves them out
        # of phase for good, and a later phone searched after it failed
        # (-12.9 dB) though heard in its sequence. So each transition's
        # times come first (they were found within 1 ms with the next
        # phone's tube still wrong), then the phone's tube in its sequence,
        # then the times again with it.
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
                arrivals, leavings = refine_transition(pool, tubes, guess, arrivals, leavings, length, index,
                                                       flow_out, true_arrivals, true_leavings, port_areas)
                if index in nasals:
                    velum[index] = velum_grid(pool, tubes, guess, arrivals, leavings, velum, length, index, flow_out)
                    print("%s: velum %.3f (true %.3f)" % (name, velum[index], true_velum[index]))
                tubes.append(guess)
                print("%s: %s at node %d" % (name, "narrowed" if index in narrowings else "closed, sealed",
                                             stops[index]))
                sys.stdout.flush()
                continue
            if index > 0:
                guess = posture_recover.nodes(hints[name], count)[0]
                arrivals, leavings = refine_transition(pool, tubes, guess, arrivals, leavings, length, index,
                                                       flow_out, true_arrivals, true_leavings, port_areas)
            begin, span = window(arrivals[index] + SETTLE, leavings[index] - CLEAR)
            context = None
            if index == first_vowel and leading:
                context = PrefixContext([], arrivals, leavings, length, seal, port_areas, jet, voiceless, leading,
                                        fricatives)
            elif index > 0 or jet > 0.0:
                context = PrefixContext(tubes, arrivals, leavings, length, seal, port_areas, jet, voiceless,
                                        fricatives=fricatives)
            value, found = tube_recover.recover(pool, hints[name], flow_out, begin, span,
                                                (true_tubes[index], length), count, context)
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
                arrivals, leavings = refine_transition(pool, tubes, found[0], arrivals, leavings, length, index,
                                                       flow_out, true_arrivals, true_leavings, port_areas, True)
                begin, span = window(arrivals[index] + SETTLE, leavings[index] - CLEAR)
                if index == first_vowel and leading:
                    polish_context = PrefixContext([], arrivals, leavings, length, seal, port_areas, jet, voiceless,
                                                   leading, fricatives)
                else:
                    polish_context = PrefixContext(tubes, arrivals, leavings, length, seal, port_areas, jet, voiceless,
                                                   fricatives=fricatives)
                value, found = tube_recover.node_search(pool, found, flow_out, begin, span, count, polish_context)
                print("%s polished: %.1f dB, tube rms %.3f" % (
                    name, value, tube_recover.tube_rms(found, (true_tubes[index], length))))
                sys.stdout.flush()
            if index == first_vowel:
                for k, (node, area) in leading.items():
                    tubes[k] = numpy.array(found[0])
                    tubes[k][node] = area
            tubes.append(found[0])
            if index > 0:
                tubes, arrivals, leavings = prefix_joint(pool, tubes, arrivals, leavings, length, flow_out,
                                                         leavings[index] - CLEAR, true_tubes, true_arrivals,
                                                         true_leavings, stops, port_areas, narrowings)
        tubes = numpy.array(tubes)
        parameters = shape.join(tubes, arrivals, leavings, velum)
        print("before the joint search: %.1f dB" % level(shape.residual(parameters, whole), reference))
        sys.stdout.flush()
        for rung in range(4):
            parameters, used = tube_recover.stage(pool, shape, parameters, whole, 30)
            found_tubes, found_arrivals, found_leavings = shape.split(parameters)
            print("joint %d: %d evaluations; %.1f dB; tubes rms %s; times off %s ms; velum %s" % (
                rung, used, level(shape.residual(parameters, whole), reference),
                " ".join("%.3f" % tube_recover.tube_rms((found_tubes[k], length), (true_tubes[k], length))
                         for k in range(phones)),
                " ".join("%.2f" % (1e3 * value) for value in numpy.concatenate([
                    found_leavings - true_leavings, found_arrivals[1:] - true_arrivals[1:]])),
                " ".join("%.3f" % v for v in shape.velum(parameters)[nasals])))
            sys.stdout.flush()
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
        costs = numpy.array(pool.map(velum_cost, tasks))
        costs[~numpy.isfinite(costs)] = numpy.inf
        best = tasks[int(numpy.argmin(costs))][3][index]
    return best


if __name__ == "__main__":
    main()
