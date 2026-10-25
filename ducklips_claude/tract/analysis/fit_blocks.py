"""Alternating least squares: the engine's parameter tree
(tract/parameter_tree.hpp) fitted block by block so that its renders match
real reference recordings, frame by frame: the NATO reference reader's 26
words (ref/nato).

Residuals. Every NATO word is rendered by the main line (say with the organ
postures, speak with the fold tube and PARAMS=tree). For each word, each of
the reference's loud frames (10 ms hop, 25 ms window, within 35 dB of the
word's loudest) carries seventeen features: its spectral shape, the levels of twelve
bands spaced evenly in log frequency from 200 Hz to 7 kHz re the frame's
own level (floored at -40 dB: the reference's coder empties quiet bands;
its microphone chain has little below 200 Hz); the frame's level re the
word's loudest (floored at -50: its noise gate empties pauses); and the
periodicity of four bands (0.08-1, 1-2, 2-4, 4-7 kHz: HNR from the
normalized correlation one period apart, the period found on the broadband
frame). Our frames are aligned to the reference's by dynamic time warping
on shape and level (both trimmed to their loud spans first), and each reference frame is met by the mean of our
frames aligned to it. A residual is (ours - reference) / 3 dB, averaged over each third of
the reference's loud frames (frame by frame, a nudged parameter re-draws
the folds' cycle-to-cycle course and the differences are mostly chance). One more
residual per word: the median |cents| by which the heard periods stray
from the intended pitch, / 20, so the tissue cannot drift off its
calibration unseen.

Blocks, alternated. The alignment is a block of its own: solved by DTW at
the start of each parameter block's turn and then held, so the residual
vector keeps its length and moves smoothly while that block is solved.
Each parameter block of FREE (moved in normalized units u; log: value =
start * exp(unit * u), linear: start + unit * u; bounded) takes its turn in
order, the rest held: the Jacobian by forward differences (0.3 units; k + 1 renders of
the 26 words, one seed), then a damped (Levenberg-Marquardt) least-squares
step on
    cost = mean(r^2) + ANCHOR * sum(u^2):
    (J'J/N + lambda diag(J'J/N) + ANCHOR) d = -(J'r/N + ANCHOR u).
Large shifts are slow by four dampers:
- lambda, per block, starts at LAMBDA_START: a rejected step multiplies it
  by 4 (and is tried again, up to TRIES), a kept one halves it, never under
  LAMBDA_FLOOR;
- every |d_i| is capped at the block's radius (RADIUS_START, x1.25 per
  kept step up to RADIUS_MAX, halved per rejected one): with unit 0.2 on a
  log scale, 0.2 is 4% per turn, 0.4 is 8%;
- the anchor's spring pulls every parameter back toward its start;
- a step is rendered before it is kept, and kept only if the cost falls; a
  render that fails, goes non-finite or loses its voice costs 100.
After every sweep the tree is also scored at a second seed it was never
fitted on.

`report` prints the cost of m64's tree and of the kinematic source.
Files in WORKDIR: scores/, reference.pkl, state.json (resume),
current.params, best.params, fit.log (one JSON line per block turn).
usage: fit_blocks.py WORKDIR report
       fit_blocks.py WORKDIR fit [SWEEPS]"""
import json
import math
import multiprocessing
import os
import pickle
import shutil
import subprocess
import sys

import librosa
import numpy
from scipy.signal import butter, sosfiltfilt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nato_snr  # noqa: E402
from parameter_files import read_tree, write_tree  # noqa: E402

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

# The dampers (see the docstring).
DIFFERENCE = 0.3     # forward-difference step, normalized units (over the renders' chance)
LAMBDA_START = 1.0
LAMBDA_FLOOR = 0.05
TRIES = 3
# (At 0.05 per turn the gain per step, about 0.005, sat under the renders'
# chance in the cost, 0.01-0.03: every step was undone.)
RADIUS_START = 0.2
RADIUS_MAX = 0.4
RADIUS_MIN = 0.02
ANCHOR = 0.01
FAILED_COST = 100.0

HOP = 160
WIDTH = 400
SHORTEST_LAG = RATE // 400
LONGEST_LAG = RATE // 60
LOUD_DB = 35.0
LEVEL_BANDS = 12
PERIODIC_BANDS = 4
FEATURES = LEVEL_BANDS + 1 + PERIODIC_BANDS
TOTAL = LEVEL_BANDS        # the column of the frame's level
SHAPE_FLOOR = -40.0        # dB re the frame: the reference's coder empties quiet bands (-100 dB)
LEVEL_FLOOR = -50.0        # dB re the loudest frame: its noise gate empties pauses
FEATURE_SCALE = 3.0    # dB
FFT_SIZE = 512
SEGMENTS = 3               # per word: beginning, middle, end of its loud frames
PITCH_SCALE = 20.0     # cents
FILTERS = nato_snr.band_filters()
BROAD = butter(4, [60, 4000], btype="band", fs=RATE, output="sos")
EDGES = numpy.geomspace(200.0, 7000.0, LEVEL_BANDS + 1)
BINS = numpy.fft.rfftfreq(FFT_SIZE, 1.0 / RATE)
WINDOW = numpy.hanning(WIDTH)


# ---------------------------------------------------------------- tree

def block_of(path):
    return path.split(".", 1)[0]


def free_blocks():
    blocks = []
    for path, _, _, _, _ in FREE:
        if block_of(path) not in blocks:
            blocks.append(block_of(path))
    return blocks


def free_value(start_value, scale, unit, low, high, u):
    value = start_value * math.exp(unit * u) if scale == "log" else start_value + unit * u
    return min(max(value, low), high)


def tree_at(start, position):
    """The start tree with the free parameters at normalized position u."""
    values = dict(start)
    for path, scale, unit, low, high in FREE:
        values[path] = free_value(values[path], scale, unit, low, high, position.get(path, 0.0))
    return [[path, values[path]] for path, _ in start]


def clamp_position(start, position):
    values = dict(start)
    clamped = dict(position)
    for path, scale, unit, low, high in FREE:
        if scale == "log":
            lowest = math.log(low / values[path]) / unit
            highest = math.log(high / values[path]) / unit
        else:
            lowest = (low - values[path]) / unit
            highest = (high - values[path]) / unit
        clamped[path] = min(max(position.get(path, 0.0), lowest), highest)
    return clamped


# ---------------------------------------------------------------- frames

def frame_features(signal):
    """Per frame: [twelve band levels re the frame's own level (dB, floored
    at SHAPE_FLOOR), the frame's level re the word's loudest (floored at
    LEVEL_FLOOR), four band HNRs (dB, -10..30)], and the frame's level re
    the loudest, unfloored."""
    banded = [sosfiltfilt(sos, signal) for sos in FILTERS]
    broad = sosfiltfilt(BROAD, signal)
    count = (len(signal) - WIDTH - LONGEST_LAG) // HOP
    if count < 5:
        return None
    features = numpy.zeros((count, FEATURES))
    total = numpy.zeros(count)
    for frame in range(count):
        at = frame * HOP
        spectrum = numpy.abs(numpy.fft.rfft(signal[at:at + WIDTH] * WINDOW, FFT_SIZE)) ** 2
        for band in range(LEVEL_BANDS):
            inside = (BINS >= EDGES[band]) & (BINS < EDGES[band + 1])
            features[frame, band] = 10.0 * math.log10(float(spectrum[inside].sum()) + 1e-14)
        # (Over the twelve bands only: the reference carries rumble under 200 Hz.)
        total[frame] = 10.0 * math.log10(float(spectrum[(BINS >= EDGES[0]) & (BINS < EDGES[-1])].sum()) + 1e-14)
        piece = broad[at:at + WIDTH]
        reach = broad[at:at + WIDTH + LONGEST_LAG]
        dots = numpy.correlate(reach, piece, "valid")
        squares = numpy.cumsum(numpy.concatenate(([0.0], reach ** 2)))
        energies = squares[WIDTH:WIDTH + len(dots)] - squares[:len(dots)]
        normalized = dots / numpy.sqrt(energies * dots[0] + 1e-30)
        lag = SHORTEST_LAG + int(numpy.argmax(normalized[SHORTEST_LAG:LONGEST_LAG + 1]))
        for index, band in enumerate(banded):
            now = band[at:at + WIDTH]
            then = band[at + lag:at + lag + WIDTH]
            correlation = numpy.dot(now, then) / math.sqrt(float(numpy.dot(now, now)) * float(numpy.dot(then, then)) + 1e-30)
            correlation = min(max(correlation, 0.0909), 0.999)
            features[frame, TOTAL + 1 + index] = 10.0 * math.log10(correlation / (1.0 - correlation))
    features[:, :LEVEL_BANDS] = numpy.maximum(features[:, :LEVEL_BANDS] - total[:, None], SHAPE_FLOOR)
    total = total - total.max()
    features[:, TOTAL] = numpy.maximum(total, LEVEL_FLOOR)
    return {"features": features, "total": total}


def reference_frames(workdir):
    """Per word, the NATO reference reader's frame features."""
    path = os.path.join(workdir, "reference.pkl")
    if os.path.exists(path):
        with open(path, "rb") as handle:
            return pickle.load(handle)
    frames = [frame_features(nato_snr.reference_nato(index)) for index in range(len(WORDS))]
    with open(path, "wb") as handle:
        pickle.dump(frames, handle)
    return frames


def loud_frames(reference):
    return numpy.where(reference["total"] > -LOUD_DB)[0]


def loud_span(frames, margin=3):
    """The first and one past the last frame within LOUD_DB of the loudest,
    widened by a margin."""
    loud = numpy.where(frames["total"] > -LOUD_DB)[0]
    return max(int(loud[0]) - margin, 0), min(int(loud[-1]) + 1 + margin, len(frames["total"]))


def align(reference, ours):
    """For each loud reference frame, the frames of ours aligned to it by DTW
    on the band shapes and level, both trimmed to their loud spans first (our renders
    carry silence the reference words do not)."""
    first, last = loud_span(reference)
    our_first, our_last = loud_span(ours)
    cost, path = librosa.sequence.dtw(X=reference["features"][first:last, :TOTAL + 1].T,
                                      Y=ours["features"][our_first:our_last, :TOTAL + 1].T, metric="cityblock")
    matched = {}
    for reference_frame, our_frame in path:
        matched.setdefault(int(reference_frame) + first, []).append(int(our_frame) + our_first)
    return [matched[frame] for frame in loud_frames(reference)]


def frame_residuals(reference, ours, alignment, scales):
    """(ours - reference) / scale, averaged over each of SEGMENTS equal runs
    of the reference's loud frames, flattened. (Frame by frame, a nudged
    parameter re-draws the folds' cycle-to-cycle course and the residuals
    are mostly chance; a segment's mean keeps what the body does always.)"""
    rows = []
    last = ours["features"].shape[0] - 1
    for frame, partners in zip(loud_frames(reference), alignment):
        partners = [min(partner, last) for partner in partners]
        rows.append((ours["features"][partners].mean(axis=0) - reference["features"][frame]) / scales)
    rows = numpy.array(rows)
    return numpy.concatenate([piece.mean(axis=0) for piece in numpy.array_split(rows, SEGMENTS)])


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


def pitch_straying(cycles_path):
    straying = []
    with open(cycles_path) as handle:
        for line in handle:
            fields = line.split()
            heard = 1000.0 / float(fields[1])
            intended = float(fields[4])
            if float(fields[6]) > 250.0 and intended > 0.0 and 0.707 < heard / intended < 1.414:
                straying.append(abs(1200.0 * math.log2(heard / intended)))
    return float(numpy.median(straying)) if len(straying) >= 5 else None


def render_word(task):
    """task = (workdir, word, params or None, seed, extra environment, tag).
    Returns (tag, word, frames or None, pitch straying or None)."""
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
            result = (tag, word, frame_features(signal), pitch_straying(cycles))
    for path in (wave, cycles):
        if os.path.exists(path):
            os.remove(path)
    return result


def render_many(pool, workdir, candidates, seed):
    """candidates: (tag, params path or None, extra environment). Returns
    {tag: ([frames per word], [straying per word])}."""
    tasks = []
    results = {}
    for tag, params, extra in candidates:
        results[tag] = ([None] * len(WORDS), [None] * len(WORDS))
        for word in WORDS:
            tasks.append((workdir, word, params, seed, extra, tag))
    for tag, word, frames, straying in pool.imap_unordered(render_word, tasks):
        index = WORDS.index(word)
        results[tag][0][index] = frames
        results[tag][1][index] = straying
    return results


def render_positions(pool, workdir, start, positions, seed):
    candidates = []
    for tag, position in positions:
        path = os.path.join(workdir, "trial_%s.params" % tag)
        write_tree(path, tree_at(start, position))
        candidates.append((tag, path, {}))
    return render_many(pool, workdir, candidates, seed)


# ---------------------------------------------------------------- residuals

def solve_alignment(references, rendered):
    """Per word: the alignment of our frames to the reference's (None when
    either is missing)."""
    return [align(reference, ours) if reference is not None and ours is not None else None
            for reference, ours in zip(references, rendered[0])]


def residuals(references, rendered, alignments):
    """The residual vector (fixed length for fixed alignments), and whether
    every word rendered with its voice."""
    frames, straying = rendered
    parts = []
    whole = True
    scales = numpy.full(FEATURES, FEATURE_SCALE)
    for reference, ours, alignment in zip(references, frames, alignments):
        if reference is None:
            continue
        if ours is None or alignment is None:
            whole = False
            parts.append(numpy.full(SEGMENTS * FEATURES, 10.0))
        else:
            parts.append(frame_residuals(reference, ours, alignment, scales))
    parts.append(numpy.array([10.0 if value is None else value / PITCH_SCALE for value in straying]))
    whole = whole and all(value is not None for value in straying)
    return numpy.concatenate(parts), whole


def cost_of(vector, whole, position):
    if not whole:
        return FAILED_COST
    return float(numpy.mean(vector ** 2)) + ANCHOR * sum(u * u for u in position.values())


def group_costs(vector):
    """Mean squared residual of the levels and of the periodicities."""
    frames = vector[:len(vector) - len(WORDS)].reshape(-1, FEATURES)  # rows: word segments
    return {"shape": float(numpy.mean(frames[:, :LEVEL_BANDS] ** 2)), "level": float(numpy.mean(frames[:, TOTAL] ** 2)),
            "periodicity": float(numpy.mean(frames[:, TOTAL + 1:] ** 2)),
            "pitch": float(numpy.mean(vector[len(vector) - len(WORDS):] ** 2))}


# ---------------------------------------------------------------- report

def report(pool, workdir):
    references = reference_frames(workdir)
    rendered = render_many(pool, workdir, [("m64", START, {}), ("kinematic", None, {"FOLDTUBE": "0"})], FIT_SEED)
    for tag in ("m64", "kinematic"):
        alignments = solve_alignment(references, rendered[tag])
        vector, whole = residuals(references, rendered[tag], alignments)
        print("%-10s cost %.4f  %s%s" % (tag, cost_of(vector, whole, {}), group_costs(vector),
                                        "" if whole else "  (words missing)"))
    # The renders' chance: m64 at a second seed, and a nudge too small to matter.
    nudge = [[path, value * 1.001 if path == "fold.cover_mass" else value] for path, value in read_tree(START)]
    write_tree(os.path.join(workdir, "nudge.params"), nudge)
    chance = render_many(pool, workdir, [("seed", START, {}), ("nudge", os.path.join(workdir, "nudge.params"), {})],
                         CHECK_SEED)
    first = residuals(references, rendered["m64"], solve_alignment(references, rendered["m64"]))[0]
    alignments = solve_alignment(references, rendered["m64"])
    second = residuals(references, chance["seed"], alignments)[0]
    nudged = render_many(pool, workdir, [("nudge", os.path.join(workdir, "nudge.params"), {})], FIT_SEED)["nudge"]
    third = residuals(references, nudged, alignments)[0]
    print("m64 vs itself at seed %s: rms residual change %.3f; cover mass x1.001: %.3f; m64 rms residual %.3f" % (
        CHECK_SEED, numpy.sqrt(numpy.mean((second - first) ** 2)), numpy.sqrt(numpy.mean((third - first) ** 2)),
        numpy.sqrt(numpy.mean(first ** 2))))


# ---------------------------------------------------------------- fit

def load_state(workdir):
    path = os.path.join(workdir, "state.json")
    if os.path.exists(path):
        with open(path) as handle:
            return json.load(handle)
    return {"sweep": 0, "turn": 0, "position": {path: 0.0 for path, _, _, _, _ in FREE},
            "damping": {block: LAMBDA_START for block in free_blocks()},
            "radius": {block: RADIUS_START for block in free_blocks()}, "best": None}


def save_state(workdir, state, start):
    path = os.path.join(workdir, "state.json")
    with open(path + ".tmp", "w") as handle:
        json.dump(state, handle)
    os.replace(path + ".tmp", path)
    write_tree(os.path.join(workdir, "current.params"), tree_at(start, state["position"]),
               "# fit_blocks.py: sweep %d, turn %d" % (state["sweep"], state["turn"]))


def solve_block(jacobian, vector, position, names, damping, radius):
    """The damped least-squares step for one block, each |d| capped."""
    count = len(vector)
    normal = jacobian.T @ jacobian / count
    gradient = jacobian.T @ vector / count + ANCHOR * numpy.array([position[name] for name in names])
    system = normal + damping * numpy.diag(numpy.diag(normal)) + ANCHOR * numpy.eye(len(names))
    step = numpy.linalg.solve(system, -gradient)
    return numpy.clip(step, -radius, radius)


def block_turn(pool, workdir, start, state, references, block, log):
    names = [path for path, _, _, _, _ in FREE if block_of(path) == block]
    position = dict(state["position"])
    # The alignment block: solved at the current tree, then held.
    positions = [("base", position)]
    for name in names:
        nudged = dict(position)
        nudged[name] += DIFFERENCE
        positions.append(("d_" + name, nudged))
    rendered = render_positions(pool, workdir, start, positions, FIT_SEED)
    alignments = solve_alignment(references, rendered["base"])
    base, whole = residuals(references, rendered["base"], alignments)
    base_cost = cost_of(base, whole, position)
    jacobian = numpy.zeros((len(base), len(names)))
    for column, name in enumerate(names):
        vector, nudged_whole = residuals(references, rendered["d_" + name], alignments)
        if nudged_whole:
            jacobian[:, column] = (vector - base) / DIFFERENCE
    record = {"sweep": state["sweep"], "block": block, "base": base_cost, "groups": group_costs(base), "tries": []}
    kept = False
    for _ in range(TRIES):
        damping = state["damping"][block]
        radius = state["radius"][block]
        step = solve_block(jacobian, base, position, names, damping, radius)
        candidate = dict(position)
        for name, change in zip(names, step):
            candidate[name] += float(change)
        candidate = clamp_position(start, candidate)
        moved = numpy.array([candidate[name] - position[name] for name in names])
        predicted = cost_of(base + jacobian @ moved, True, candidate)
        tried = render_positions(pool, workdir, start, [("step", candidate)], FIT_SEED)["step"]
        vector, tried_whole = residuals(references, tried, alignments)
        tried_cost = cost_of(vector, tried_whole, candidate)
        record["tries"].append({"damping": damping, "radius": radius, "cost": tried_cost, "predicted": predicted,
                                "step": {name: round(candidate[name] - position[name], 4) for name in names}})
        if tried_cost < base_cost:
            state["position"] = candidate
            state["damping"][block] = max(damping * 0.5, LAMBDA_FLOOR)
            state["radius"][block] = min(radius * 1.25, RADIUS_MAX)
            record["kept"] = tried_cost
            record["kept_groups"] = group_costs(vector)
            kept = True
            break
        state["damping"][block] = damping * 4.0
        state["radius"][block] = max(radius * 0.5, RADIUS_MIN)
    log.write(json.dumps(record) + "\n")
    log.flush()
    return base_cost, record.get("kept", base_cost), kept


def fit(pool, workdir, sweeps):
    references = reference_frames(workdir)
    start = read_tree(START)
    state = load_state(workdir)
    blocks = free_blocks()
    log = open(os.path.join(workdir, "fit.log"), "a")
    for _ in range(sweeps):
        while state["turn"] < len(blocks):
            block = blocks[state["turn"]]
            before, after, kept = block_turn(pool, workdir, start, state, references, block, log)
            state["turn"] += 1
            save_state(workdir, state, start)
            if state["best"] is None or after < state["best"]:
                state["best"] = after
                shutil.copy(os.path.join(workdir, "current.params"), os.path.join(workdir, "best.params"))
                save_state(workdir, state, start)
            print("sweep %d %-14s %.4f -> %.4f %s (lambda %.3g, radius %.3f)" % (
                state["sweep"], block, before, after, "kept" if kept else "undone", state["damping"][block],
                state["radius"][block]), flush=True)
        check = render_positions(pool, workdir, start, [("check", state["position"])], CHECK_SEED)["check"]
        alignments = solve_alignment(references, check)
        vector, whole = residuals(references, check, alignments)
        check_cost = cost_of(vector, whole, state["position"])
        log.write(json.dumps({"sweep": state["sweep"], "check": check_cost, "check_groups": group_costs(vector)}) + "\n")
        log.flush()
        print("sweep %d done; at seed %s the cost is %.4f" % (state["sweep"], CHECK_SEED, check_cost), flush=True)
        state["sweep"] += 1
        state["turn"] = 0
        save_state(workdir, state, start)


def main():
    workdir = os.path.abspath(sys.argv[1])
    mode = sys.argv[2]
    os.makedirs(workdir, exist_ok=True)
    prepare_scores(workdir)
    with multiprocessing.Pool(WORKERS) as pool:
        if mode == "report":
            report(pool, workdir)
        elif mode == "fit":
            fit(pool, workdir, int(sys.argv[3]) if len(sys.argv) > 3 else 1)


if __name__ == "__main__":
    main()
