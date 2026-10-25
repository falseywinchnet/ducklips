"""REJECTED 2026-10-16 (tract/NOTES.md): its loss scored the kinematic source
(Perfect Paul) as near the humans as they are to each other. Kept as the
record; the fit is fit_blocks.py.

Behavioural fit of the engine's parameter tree (tract/parameter_tree.hpp)
to the human readings of the NATO alphabet, block by block, by damped
gradient descent.

Behaviours. Every NATO word is rendered by the main line (say with the
organ postures, speak with the fold tube and PARAMS=tree) and measured as
the human readings are, by the same code: over every glottal cycle in the
loud stretches (nato_snr.loud_cycles), the cycle-coherence HNR per band
(median: "hnr"; interquartile spread across cycles: "spread"), the
perturbation of period and amplitude (RAP, log; shimmer, dB:
"perturbation"), the voiced energy's balance across the bands (dB re its
total: "balance"), the pitch span (p90 - p10, semitones: "span"), local
clicks above 5 kHz (log(1 + per second): "clicks"), and the word's level
contours in three bands, compared by dynamic time warping ("contour").
Ours also has "pitch": how far the heard periods stray from the intended
pitch (median |cents|; the engine's CYCLES), so the tissue cannot drift
off its calibration unseen.

Loss. Each behaviour's difference from the humans' mean (the owner's two
passes and the NATO reader) is divided by how much those three readings
differ among themselves (median over words of their standard deviation,
with a floor), passed through a Huber function, averaged within its group,
then over groups and words. A contour scores the excess of our warped
distance over the humans' distance to each other. So a loss near the
humans' own (leave-one-out: one reading scored against the other two) is
"as far from them as they are from each other". `validate` prints that
floor beside our voice and two rejected sources before anything is fitted.

Descent. The free parameters (FREE: path, log or linear, unit, bounds) are
moved in normalized units u (log: value = start * exp(unit * u); linear:
start + unit * u). Each iteration takes one block of the tree in turn and
estimates the gradient of the loss within it from two renders (SPSA: every
parameter of the block nudged by +-c at once, signs at random, the same
seed for both), smoothed by Adam-like moments. Shifts are damped three
ways, so that a large change can only come from a gradient that persists:
- each parameter's step is capped by its block's trust radius (at most
  RADIUS_MAX units per accepted step; with unit 0.2 on a log scale, 2%);
- a spring (ANCHOR) pulls every parameter back toward its start: the
  loss carries ANCHOR/2 * sum(u^2);
- the step is rendered before it is kept: a step that does not lower the
  loss (or whose render fails, goes non-finite or loses its voicing) is
  undone, its block's radius halved and its momentum cleared; a kept step
  grows the radius by a quarter, up to RADIUS_MAX.
Every few iterations the current tree is also scored at a second seed it
was never fitted on; a fit that only learns one seed's noise shows there.

Files in WORKDIR: scores/ (say, once), references.pkl, state.json (resume),
current.params, best.params, fit.log (one JSON line per iteration).
usage: fit_spsa_summary.py WORKDIR validate
       fit_spsa_summary.py WORKDIR fit [ITERATIONS]"""
import json
import math
import multiprocessing
import os
import pickle
import random
import shutil
import subprocess
import sys

import librosa
import numpy
from scipy.signal import butter, sosfiltfilt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nato_snr  # noqa: E402

nato_compare = nato_snr.nato_compare
ROOT = nato_snr.ROOT
WORDS = nato_compare.WORDS
RATE = nato_snr.RATE
START = os.path.join(ROOT, "tract", "data", "voice.params")
FIT_SEED = "73021"
CHECK_SEED = "51117"
WORKERS = 10

# path, scale ("log" or "linear"), unit (per normalized step of 1), low, high.
FREE = (
    ("breath.lung_follow", "log", 0.3, 0.002, 0.03),
    ("fold.depth", "log", 0.2, 0.15, 0.6),
    ("fold.cover_mass", "log", 0.2, 0.03, 0.3),
    ("fold.body_mass", "log", 0.2, 0.1, 1.0),
    ("fold.cover_stiffness", "log", 0.2, 1.0e4, 1.2e5),
    ("fold.coupling_stiffness", "log", 0.2, 3.0e4, 3.0e5),
    ("fold.support_stiffness", "log", 0.2, 5.0e4, 5.0e5),
    ("fold.vertical_stiffness", "log", 0.2, 2.0e4, 2.0e5),
    ("fold.cover_loss", "log", 0.3, 0.02, 0.5),
    ("fold.support_loss", "log", 0.3, 0.05, 1.5),
    ("fold.contact_stiffness", "log", 0.3, 2.0e5, 4.0e6),
    ("fold.contact_softness", "log", 0.3, 1.0e-4, 3.0e-3),
    ("fold.muscle_fluctuation", "log", 0.3, 1.0e-3, 0.05),
    ("fold.muscle_corner", "log", 0.3, 20.0, 500.0),
    ("fold.natural_pitch", "log", 0.1, 120.0, 320.0),
    ("fold.pitch_per_tension", "log", 0.1, 90.0, 340.0),
    ("fold.elongation", "log", 0.3, 10.0, 300.0),
    ("fold.swing_per_pascal", "log", 0.2, 1.0e-6, 2.0e-5),
    ("posture.voicing_lower", "linear", 0.002, 0.002, 0.03),
    ("posture.voicing_upper", "linear", 0.002, -0.02, 0.02),
    ("posture.voicing_chink", "log", 0.3, 2.0e-5, 5.0e-3),
    ("posture.closure_press", "linear", 0.002, 0.0, 0.02),
    ("ear.pitch_feedback", "log", 0.3, 0.01, 0.3),
    ("variation.wall", "log", 0.3, 0.01, 0.4),
    ("variation.wander_gain", "log", 0.3, 0.1, 3.0),
    ("glottal_noise.aspiration", "log", 0.2, 1.0, 30.0),
    ("glottal_noise.closing_turbulence", "linear", 0.2, 0.0, 3.0),
    ("glottal_noise.fold_turbulence", "log", 0.3, 0.1, 3.0),
    ("oral_noise.frication", "log", 0.2, 0.3, 6.0),
    ("oral_noise.jet_attack", "log", 0.3, 1.0e-3, 0.03),
    ("walls.skin_gain", "log", 0.3, 0.1, 4.0),
    ("walls.skin_corner", "log", 0.3, 60.0, 1000.0),
    ("walls.stiffness_scale", "log", 0.3, 0.3, 3.0),
    ("walls.loss_scale", "log", 0.3, 0.3, 3.0),
    ("walls.closure_pull", "linear", 0.05, 0.0, 0.6),
    ("branches.piriform_loss", "log", 0.3, 0.1, 1.5),
    ("branches.piriform_port", "log", 0.3, 0.15, 2.0),
    ("branches.piriform_scale", "log", 0.2, 0.25, 1.5),
    ("branches.sinus_neck_scale", "log", 0.3, 0.05, 1.0),
    ("tract.hop_loss", "linear", 0.005, 0.0, 0.1),
    ("tract.nose_hop_loss", "linear", 0.05, 0.05, 0.8),
    ("mouth.release_speed", "log", 0.2, 0.5, 3.0),
)

# The descent's damping (see the docstring).
PERTURBATION = 0.1   # c, normalized units
RATE_OF_STEP = 0.05  # Adam step size, normalized units
RADIUS_START = 0.05
RADIUS_MAX = 0.1
RADIUS_MIN = 0.005
ANCHOR = 0.01
FIRST_MOMENT = 0.8
SECOND_MOMENT = 0.98
CHECK_EVERY = 10

# Behaviour groups compared feature by feature, with each feature's floor
# on the human spread (the scale never falls under it).
VECTOR_GROUPS = (("hnr", 1.0), ("spread", 1.0), ("perturbation", 0.1), ("balance", 1.0), ("span", 0.5),
                 ("clicks", 0.2))
PITCH_SCALE = 20.0  # cents: the intended pitch's tolerance
FILTERS = nato_snr.band_filters()
CLICK_FILTER = butter(4, 5000, btype="high", fs=RATE, output="sos")


# ---------------------------------------------------------------- block files

def split_words(text):
    words = []
    for line in text.splitlines():
        line = line.split("#", 1)[0].replace("{", " { ").replace("}", " } ")
        words.extend(line.split())
    return words


def read_tree(path):
    """The leaves of a block file as [dotted path, value], in file order."""
    with open(path) as handle:
        words = split_words(handle.read())
    entries = []
    open_blocks = []
    at = 0
    while at < len(words):
        word = words[at]
        if word == "}":
            open_blocks.pop()
            at += 1
        elif at + 1 < len(words) and words[at + 1] == "{":
            open_blocks.append(word)
            at += 2
        else:
            entries.append([".".join(open_blocks + [word]), float(words[at + 1])])
            at += 2
    return entries


def write_tree(path, entries):
    lines = []
    open_blocks = []
    for name, value in entries:
        parts = name.split(".")
        shared = 0
        while shared < len(open_blocks) and shared + 1 < len(parts) and open_blocks[shared] == parts[shared]:
            shared += 1
        while len(open_blocks) > shared:
            open_blocks.pop()
            lines.append("  " * len(open_blocks) + "}")
        while len(open_blocks) + 1 < len(parts):
            lines.append("  " * len(open_blocks) + parts[len(open_blocks)] + " {")
            open_blocks.append(parts[len(open_blocks)])
        lines.append("  " * len(open_blocks) + "%s %.10g" % (parts[-1], value))
    while open_blocks:
        open_blocks.pop()
        lines.append("  " * len(open_blocks) + "}")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def block_of(path):
    return path.split(".", 1)[0]


def free_blocks():
    blocks = []
    for path, _, _, _, _ in FREE:
        if block_of(path) not in blocks:
            blocks.append(block_of(path))
    return blocks


def tree_at(start, position):
    """The start tree with the free parameters at normalized position u
    (a dict path -> u), clamped to their bounds."""
    values = dict(start)
    for path, scale, unit, low, high in FREE:
        u = position.get(path, 0.0)
        if scale == "log":
            value = values[path] * math.exp(unit * u)
        else:
            value = values[path] + unit * u
        values[path] = min(max(value, low), high)
    return [[path, values[path]] for path, _ in start]


def clamp_position(start, position):
    """u pulled back inside the bounds (so a parameter at a wall does not
    keep momentum pushing into it)."""
    values = dict(start)
    clamped = {}
    for path, scale, unit, low, high in FREE:
        u = position.get(path, 0.0)
        if scale == "log":
            lowest = math.log(low / values[path]) / unit
            highest = math.log(high / values[path]) / unit
        else:
            lowest = (low - values[path]) / unit
            highest = (high - values[path]) / unit
        clamped[path] = min(max(u, lowest), highest)
    return clamped


# ---------------------------------------------------------------- behaviours

def cycle_runs(found):
    """Runs of adjacent cycles (each starting where the one before ended)."""
    runs = []
    run = []
    for start, period, match in found:
        if run and abs(start - (run[-1][0] + int(round(run[-1][1])))) > 1:
            runs.append(run)
            run = []
        run.append((start, period, match))
    if run:
        runs.append(run)
    return runs


def perturbation(signal, found):
    """RAP (log10 of the fraction) and shimmer (mean |dB| between
    neighbouring cycles' peak-to-peak amplitudes) over runs of matched
    cycles."""
    rap = []
    shimmer = []
    for run in cycle_runs(found):
        good = [cycle for cycle in run if cycle[2] > 0.5]
        if len(good) < 3:
            continue
        periods = numpy.array([cycle[1] for cycle in good])
        for index in range(1, len(periods) - 1):
            rap.append(abs(periods[index] - periods[index - 1:index + 2].mean()) / periods.mean())
        amplitudes = []
        for start, period, _ in good:
            piece = signal[start:start + int(round(period))]
            amplitudes.append(piece.max() - piece.min() + 1e-12)
        amplitudes = numpy.array(amplitudes)
        shimmer.extend(numpy.abs(20.0 * numpy.log10(amplitudes[1:] / amplitudes[:-1])))
    if len(rap) < 5:
        return None
    return [math.log10(float(numpy.mean(rap)) + 1e-6), float(numpy.mean(shimmer))]


def click_rate(signal, factor=8.0):
    high = sosfiltfilt(CLICK_FILTER, signal)
    width = int(0.02 * RATE)
    envelope = numpy.sqrt(numpy.convolve(high ** 2, numpy.ones(width) / width, "same")) + 1e-9
    level = numpy.sqrt(numpy.convolve(signal ** 2, numpy.ones(width) / width, "same"))
    loud = level > level.max() * 10.0 ** (-40.0 / 20.0)
    times = []
    for index in numpy.where((numpy.abs(high) > factor * envelope) & loud)[0]:
        if not times or index - times[-1] > 0.003 * RATE:
            times.append(index)
    seconds = max(float(loud.sum()) / RATE, 0.1)
    return math.log1p(len(times) / seconds)


def behaviours(signal):
    """The behaviour groups of one word (None when it has too little voice)."""
    found = nato_snr.loud_cycles(signal)
    if len(found) < 8:
        return None
    banded = [sosfiltfilt(sos, signal) for sos in FILTERS]
    hnr = []
    spread = []
    energy = []
    for band in banded:
        values = []
        total = 0.0
        for start, period, _ in found:
            size = int(round(period))
            following = start + size
            if following + size > len(signal):
                continue
            now = band[start:following]
            then = band[following:following + size]
            correlation = numpy.dot(now, then) / numpy.sqrt(numpy.dot(now, now) * numpy.dot(then, then) + 1e-20)
            correlation = min(max(correlation, 1e-3), 0.9999)
            values.append(10.0 * numpy.log10(correlation / (1.0 - correlation)))
            total += float(numpy.dot(now, now))
        if len(values) < 4:
            return None
        quartiles = numpy.percentile(values, [25, 50, 75])
        hnr.append(float(quartiles[1]))
        spread.append(float(quartiles[2] - quartiles[0]))
        energy.append(total + 1e-20)
    energy = numpy.array(energy)
    balance = list(10.0 * numpy.log10(energy / energy.sum()))
    shaking = perturbation(signal, found)
    if shaking is None:
        return None
    pitches = numpy.array([RATE / period for _, period, match in found if match > 0.5])
    if len(pitches) < 5:
        return None
    span = 12.0 * math.log2(numpy.percentile(pitches, 90) / numpy.percentile(pitches, 10))
    return {"hnr": hnr, "spread": spread, "perturbation": shaking, "balance": balance, "span": [span],
            "clicks": [click_rate(signal)], "contour": nato_compare.contours(signal)}


def pitch_straying(cycles_path):
    """Median |cents| between the heard periods and the intended pitch,
    over voiced cycles (lungs over 250 Pa, within half an octave)."""
    straying = []
    with open(cycles_path) as handle:
        for line in handle:
            fields = line.split()
            heard = 1000.0 / float(fields[1])
            intended = float(fields[4])
            if float(fields[6]) > 250.0 and intended > 0.0 and 0.707 < heard / intended < 1.414:
                straying.append(abs(1200.0 * math.log2(heard / intended)))
    return float(numpy.median(straying)) if len(straying) >= 5 else None


# ---------------------------------------------------------------- rendering

def prepare_scores(workdir):
    folder = os.path.join(workdir, "scores")
    os.makedirs(folder, exist_ok=True)
    for word in WORDS:
        score = os.path.join(folder, word + ".score")
        if os.path.exists(score):
            continue
        subprocess.run([os.path.join(ROOT, "voice", "build", "say"), word + ".", "/dev/null", "--score", score,
                        "--plan-only", "--pressure", "750", "--phones",
                        os.path.join(ROOT, "tract", "data", "phones_tract.txt"), "--organs",
                        os.path.join(ROOT, "tract", "data", "organs.txt")],
                       cwd=os.path.join(ROOT, "voice"), check=True, stdout=subprocess.DEVNULL)


def render_word(task):
    """Renders and measures one word: task = (workdir, word, params or None,
    seed, extra environment, tag). Returns (tag, word, behaviours, pitch
    straying); behaviours None when the render failed or lost its voice."""
    workdir, word, params, seed, extra, tag = task
    folder = os.path.join(workdir, "renders", tag)
    os.makedirs(folder, exist_ok=True)
    wave = os.path.join(folder, word + ".wav")
    cycles = os.path.join(folder, word + ".cycles")
    environment = dict(os.environ)
    environment["FOLDTUBE"] = "1"
    environment["CYCLES"] = cycles
    if params is not None:
        environment["PARAMS"] = params
    environment.update(extra)
    finished = subprocess.run([os.path.join(ROOT, "tract", "build", "speak"),
                               os.path.join(workdir, "scores", word + ".score"), wave, seed],
                              env=environment, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = (tag, word, None, None)
    if finished.returncode == 0:
        signal, _ = librosa.load(wave, sr=RATE)
        if numpy.all(numpy.isfinite(signal)) and numpy.abs(signal).max() > 0.0:
            result = (tag, word, behaviours(signal), pitch_straying(cycles))
    for path in (wave, cycles):
        if os.path.exists(path):
            os.remove(path)
    return result


def reference_behaviours(workdir):
    """[owner pass 1, owner pass 2, NATO]: per reading, a list over words."""
    path = os.path.join(workdir, "references.pkl")
    if os.path.exists(path):
        with open(path, "rb") as handle:
            return pickle.load(handle)
    readings = [[], [], []]
    for index in range(len(WORDS)):
        readings[0].append(behaviours(nato_snr.reference_owner(index, 0)))
        readings[1].append(behaviours(nato_snr.reference_owner(index, 1)))
        readings[2].append(behaviours(nato_snr.reference_nato(index)))
    with open(path, "wb") as handle:
        pickle.dump(readings, handle)
    return readings


# ---------------------------------------------------------------- loss

def huber(z):
    z = numpy.abs(numpy.asarray(z, dtype=float))
    return numpy.where(z < 1.0, 0.5 * z * z, z - 0.5)


def human_scales(readings):
    """Per group, per feature: the median over words of the three readings'
    standard deviation, floored."""
    scales = {}
    for group, floor in VECTOR_GROUPS:
        deviations = []
        for index in range(len(WORDS)):
            values = [reading[index][group] for reading in readings if reading[index] is not None]
            if len(values) >= 2:
                deviations.append(numpy.std(numpy.array(values), axis=0, ddof=1))
        scales[group] = numpy.maximum(numpy.median(numpy.array(deviations), axis=0), floor)
    return scales


def contour_distance(first, second):
    cost, path = librosa.sequence.dtw(X=first, Y=second, metric="cityblock")
    return float(cost[-1, -1]) / len(path) / first.shape[0]


def word_terms(ours, humans, scales):
    """The groups' Huber means for one word against a set of human readings."""
    terms = {}
    for group, _ in VECTOR_GROUPS:
        target = numpy.mean(numpy.array([human[group] for human in humans]), axis=0)
        terms[group] = float(numpy.mean(huber((numpy.array(ours[group]) - target) / scales[group])))
    ours_distance = numpy.mean([contour_distance(ours["contour"], human["contour"]) for human in humans])
    between = []
    for first in range(len(humans)):
        for second in range(first + 1, len(humans)):
            between.append(contour_distance(humans[first]["contour"], humans[second]["contour"]))
    if between:
        excess = (ours_distance - numpy.mean(between)) / numpy.mean(between)
        terms["contour"] = float(huber(max(excess, 0.0)))
    return terms


def loss_of(candidate, readings, scales, straying=None):
    """Mean over words of the mean over groups; (loss, per-group means).
    A word without voice scores 10 in every group."""
    per_group = {}
    word_losses = []
    for index in range(len(WORDS)):
        humans = [reading[index] for reading in readings if reading[index] is not None]
        if len(humans) < 2:
            continue
        if candidate[index] is None:
            terms = {group: 10.0 for group, _ in VECTOR_GROUPS}
            terms["contour"] = 10.0
        else:
            terms = word_terms(candidate[index], humans, scales)
        if straying is not None:
            terms["pitch"] = 10.0 if straying[index] is None else float(huber(straying[index] / PITCH_SCALE))
        for group, value in terms.items():
            per_group.setdefault(group, []).append(value)
        word_losses.append(numpy.mean(list(terms.values())))
    return float(numpy.mean(word_losses)), {group: float(numpy.mean(values)) for group, values in per_group.items()}


# ---------------------------------------------------------------- evaluation

def render_many(pool, workdir, candidates, seed):
    """candidates: list of (tag, params path or None, extra environment).
    Returns {tag: ([behaviours per word], [straying per word])}."""
    tasks = []
    for tag, params, extra in candidates:
        for word in WORDS:
            tasks.append((workdir, word, params, seed, extra, tag))
    results = {}
    for tag, _, _ in candidates:
        results[tag] = ([None] * len(WORDS), [None] * len(WORDS))
    for tag, word, measured, straying in pool.imap_unordered(render_word, tasks):
        index = WORDS.index(word)
        results[tag][0][index] = measured
        results[tag][1][index] = straying
    return results


def anchor_penalty(position):
    return 0.5 * ANCHOR * sum(u * u for u in position.values())


def score_positions(pool, workdir, start, positions, seed, readings, scales):
    """positions: list of (tag, u). Returns {tag: (total loss, behaviour
    loss, per-group)}."""
    candidates = []
    for tag, position in positions:
        path = os.path.join(workdir, "trial_%s.params" % tag)
        write_tree(path, tree_at(start, position))
        candidates.append((tag, path, {}))
    rendered = render_many(pool, workdir, candidates, seed)
    scored = {}
    for tag, position in positions:
        measured, straying = rendered[tag]
        loss, groups = loss_of(measured, readings, scales, straying)
        scored[tag] = (loss + anchor_penalty(position), loss, groups)
    return scored


def show_groups(groups):
    return " ".join("%s %.3f" % (group, value) for group, value in sorted(groups.items()))


def validate(pool, workdir):
    readings = reference_behaviours(workdir)
    scales = human_scales(readings)
    print("human scales:", {group: numpy.round(scale, 3).tolist() for group, scale in scales.items()})
    start = read_tree(START)
    clean = [[path, 0.0 if path == "glottal_noise.aspiration" else value] for path, value in start]
    write_tree(os.path.join(workdir, "clean.params"), clean)
    candidates = [("fold_tube", START, {}),
                  ("no_aspiration", os.path.join(workdir, "clean.params"), {}),
                  ("kinematic", None, {"FOLDTUBE": "0"})]
    rendered = render_many(pool, workdir, candidates, FIT_SEED)
    names = ("owner pass 1", "owner pass 2", "NATO")
    print("\nloss against the other two readings (leave one out), then against all three")
    for left_out in range(3):
        others = [readings[index] for index in range(3) if index != left_out]
        loss, groups = loss_of(readings[left_out], others, scales)
        print("%-14s vs others: %.3f  %s" % (names[left_out], loss, show_groups(groups)))
        for tag, _, _ in candidates:
            loss, groups = loss_of(rendered[tag][0], others, scales)
            print("  %-12s vs same: %.3f  %s" % (tag, loss, show_groups(groups)))
    for tag, _, _ in candidates:
        measured, straying = rendered[tag]
        loss, groups = loss_of(measured, readings, scales, straying)
        print("%-14s vs all three: %.3f  %s" % (tag, loss, show_groups(groups)))


def load_state(workdir, start):
    path = os.path.join(workdir, "state.json")
    if os.path.exists(path):
        with open(path) as handle:
            return json.load(handle)
    blocks = free_blocks()
    return {"iteration": 0, "position": {path: 0.0 for path, _, _, _, _ in FREE},
            "first": {path: 0.0 for path, _, _, _, _ in FREE}, "second": {path: 0.0 for path, _, _, _, _ in FREE},
            "steps": {path: 0 for path, _, _, _, _ in FREE},
            "radius": {block: RADIUS_START for block in blocks}, "loss": None, "best": None}


def save_state(workdir, state, start):
    path = os.path.join(workdir, "state.json")
    with open(path + ".tmp", "w") as handle:
        json.dump(state, handle)
    os.replace(path + ".tmp", path)
    write_tree(os.path.join(workdir, "current.params"), tree_at(start, state["position"]))


def fit(pool, workdir, iterations):
    readings = reference_behaviours(workdir)
    scales = human_scales(readings)
    start = read_tree(START)
    state = load_state(workdir, start)
    generator = random.Random(int(FIT_SEED) + state["iteration"])
    blocks = free_blocks()
    log = open(os.path.join(workdir, "fit.log"), "a")
    if state["loss"] is None:
        scored = score_positions(pool, workdir, start, [("now", state["position"])], FIT_SEED, readings, scales)
        state["loss"] = scored["now"][0]
        state["best"] = scored["now"][0]
        log.write(json.dumps({"iteration": 0, "loss": scored["now"][1], "groups": scored["now"][2]}) + "\n")
        log.flush()
        save_state(workdir, state, start)
        shutil.copy(os.path.join(workdir, "current.params"), os.path.join(workdir, "best.params"))
        print("start: %.4f  %s" % (scored["now"][0], show_groups(scored["now"][2])), flush=True)
    for _ in range(iterations):
        state["iteration"] += 1
        block = blocks[(state["iteration"] - 1) % len(blocks)]
        names = [path for path, _, _, _, _ in FREE if block_of(path) == block]
        position = dict(state["position"])
        signs = {path: generator.choice((-1.0, 1.0)) for path in names}
        plus = dict(position)
        minus = dict(position)
        for path in names:
            plus[path] += PERTURBATION * signs[path]
            minus[path] -= PERTURBATION * signs[path]
        scored = score_positions(pool, workdir, start, [("plus", plus), ("minus", minus)], FIT_SEED, readings, scales)
        difference = scored["plus"][0] - scored["minus"][0]
        radius = state["radius"][block]
        candidate = dict(position)
        for path in names:
            gradient = difference / (2.0 * PERTURBATION) * signs[path]
            state["steps"][path] += 1
            count = state["steps"][path]
            state["first"][path] = FIRST_MOMENT * state["first"][path] + (1.0 - FIRST_MOMENT) * gradient
            state["second"][path] = SECOND_MOMENT * state["second"][path] + (1.0 - SECOND_MOMENT) * gradient ** 2
            first = state["first"][path] / (1.0 - FIRST_MOMENT ** count)
            second = state["second"][path] / (1.0 - SECOND_MOMENT ** count)
            step = -RATE_OF_STEP * first / (math.sqrt(second) + 1e-8)
            candidate[path] += min(max(step, -radius), radius)
        candidate = clamp_position(start, candidate)
        tried = score_positions(pool, workdir, start, [("step", candidate)], FIT_SEED, readings, scales)["step"]
        kept = tried[0] < state["loss"]
        if kept:
            state["position"] = candidate
            state["loss"] = tried[0]
            state["radius"][block] = min(radius * 1.25, RADIUS_MAX)
        else:
            state["radius"][block] = max(radius * 0.5, RADIUS_MIN)
            for path in names:
                state["first"][path] = 0.0
        record = {"iteration": state["iteration"], "block": block, "kept": kept, "total": tried[0],
                  "loss": tried[1], "now": state["loss"], "plus": scored["plus"][0], "minus": scored["minus"][0],
                  "radius": state["radius"][block], "groups": tried[2],
                  "moved": {path: round(candidate[path] - position[path], 4) for path in names}}
        if state["iteration"] % CHECK_EVERY == 0:
            check = score_positions(pool, workdir, start, [("check", state["position"])], CHECK_SEED, readings,
                                    scales)["check"]
            record["check"] = check[1]
            record["check_groups"] = check[2]
        log.write(json.dumps(record) + "\n")
        log.flush()
        save_state(workdir, state, start)
        if kept and state["loss"] < state["best"]:
            state["best"] = state["loss"]
            shutil.copy(os.path.join(workdir, "current.params"), os.path.join(workdir, "best.params"))
            save_state(workdir, state, start)
        print("%4d %-14s %s %.4f (now %.4f, radius %.3f)%s" % (
            state["iteration"], block, "kept" if kept else "undone", tried[0], state["loss"],
            state["radius"][block], "  check %.4f" % record["check"] if "check" in record else ""), flush=True)


def main():
    workdir = os.path.abspath(sys.argv[1])
    mode = sys.argv[2]
    os.makedirs(workdir, exist_ok=True)
    prepare_scores(workdir)
    with multiprocessing.Pool(WORKERS) as pool:
        if mode == "validate":
            validate(pool, workdir)
        elif mode == "fit":
            fit(pool, workdir, int(sys.argv[3]) if len(sys.argv) > 3 else 100)


if __name__ == "__main__":
    main()
