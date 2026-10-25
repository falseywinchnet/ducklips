"""The engine as a network of parameter neurons, for policy optimization.

This is not a network that is trained to compute anything. Each neuron
holds one parameter of the engine, in three groups: "body", every leaf of
the engine's parameter tree (tract/parameter_tree.hpp, data/voice.params);
"organs", every vowel posture coordinate (data/organs.txt); "plan", every
leaf of the planner's tree (voice/plan_parameters.hpp, voice/data/
plan.params: timing, glottis, breath, release, consonant shapes, a duration
factor per phone). A network is built over the groups being fitted; the
others stay at their files' values. A neuron's state is a latent z, 0 at the
start; its output is the parameter, through a bounded, saturating
activation, so that no step anywhere can throw a parameter far:
- "log": value = exp(lerp(log low, log high, sigmoid(z0 + z)));
- "linear": value = lerp(low, high, sigmoid(z0 + z));
- "onset": value = high * tanh(max(z, 0)), for a cause that starts off
  (0) and can only be switched on.
z0 puts the start value at z = 0. The slope is steepest at the middle of a
range and falls toward its bounds: a parameter that is driven toward a wall
moves ever more slowly. The optimizer (policy_fit.py) moves the latents.

An evaluation renders the 26 NATO words through say (the network's vowel
postures and plan tree) and speak (FOLDTUBE=1, its body tree) and scores
them against the NATO reference recordings: fit_blocks.py's frame features
pooled over SEGMENTS runs of each word, plus timing (the word's length and
its inner timing against the reference's, which DTW alone forgives).

Left out of the body tree: kinematic.* (they act only on the kinematic
source, not the fold tube spoken here) and glottal_noise.flux_corner (0 is a
switch, "shaped like the jet noise", not the bottom of a range). Left out of
the plan tree: intonation (pitch.*: the cost does not measure pitch), the
area planner's mouth.* (unused by the articulatory body), the switch
glottis.onset_gap (-1: off), and lead_in and pause_pressure (silence only)."""
import math
import os
import subprocess

import librosa
import numpy

import fit_blocks
from parameter_files import read_tree, write_tree

ROOT = fit_blocks.ROOT
WORDS = fit_blocks.WORDS
RATE = fit_blocks.RATE
# (ENGINE_TREE, ENGINE_ORGANS, ENGINE_PLAN in the environment start from other files.)
TREE = os.environ.get("ENGINE_TREE", os.path.join(ROOT, "tract", "data", "voice.params"))
ORGANS = os.environ.get("ENGINE_ORGANS", os.path.join(ROOT, "tract", "data", "organs.txt"))
PLAN = os.environ.get("ENGINE_PLAN", os.path.join(ROOT, "voice", "data", "plan.params"))
SEGMENTS = 6
ORGAN_KEYS = ("jaw", "tx", "ty", "tipx", "tipy", "ap", "pro", "lw", "drop")
LEFT_OUT = ("kinematic.amplitude", "kinematic.phase_lag", "kinematic.softness", "kinematic.jitter",
            "glottal_noise.flux_corner")

# Tree leaves whose range is not the default (log, start / 2.5 .. start * 2.5):
# path -> (activation, low, high).
TREE_RANGES = {
    "fold.natural_pitch": ("log", 150.0, 260.0),
    "fold.pitch_per_tension": ("log", 110.0, 260.0),
    "fold.elongation_floor": ("linear", 0.7, 1.0),
    "posture.voicing_lower": ("linear", 0.002, 0.03),
    "posture.voicing_upper": ("linear", -0.02, 0.02),
    "posture.closure_press": ("linear", -0.002, 0.015),
    "variation.breath": ("onset", 0.0, 0.06),
    "variation.tension": ("onset", 0.0, 0.03),
    "variation.posture": ("onset", 0.0, 0.003),
    "variation.tremor": ("onset", 0.0, 0.01),
    "variation.adductor_tremor": ("onset", 0.0, 0.0015),
    "variation.breath_tremor": ("onset", 0.0, 0.03),
    "variation.wall_shared": ("linear", 0.0, 1.0),
    "glottal_noise.closing_turbulence": ("onset", 0.0, 2.0),
    "glottal_noise.jet_flux": ("onset", 0.0, 0.08),
    "oral_noise.jet_spread": ("onset", 0.0, 0.3),
    "oral_noise.jet_corner_ceiling": ("log", 6000.0, 16000.0),
    "walls.closure_pull": ("linear", 0.0, 0.6),
    "walls.closure_relaxation": ("onset", 0.0, 0.5),
    "branches.piriform_smoothing": ("onset", 0.0, 0.3),
    "tract.hop_loss": ("log", 0.005, 0.08),
    "tract.nose_hop_loss": ("linear", 0.1, 0.8),
}
DEFAULT_SPAN = 2.5

PLAN_LEFT_OUT_BLOCKS = ("pitch", "mouth")
PLAN_LEFT_OUT = ("glottis.onset_gap", "timing.lead_in", "breath.pause_pressure")
# Plan leaves whose range is not their block's default: path -> (activation, low, high).
PLAN_RANGES = {
    "timing.rate": ("log", 0.7, 1.6),
    "timing.seal_lead": ("linear", 0.0, 1.0),
    "timing.vowel_hold": ("linear", 0.2, 0.9),
    "timing.affricate_hiss_at": ("linear", 0.1, 0.6),
    "timing.vowel_arrival": ("linear", 0.4, 1.0),
    "timing.sonorant_devoicing": ("linear", 0.0, 1.0),
    "timing.fricative_hold": ("onset", 0.0, 0.5),
    "glottis.voicing_gap": ("linear", -0.015, 0.005),
    "glottis.unstressed_gap": ("linear", -0.005, 0.01),
    "glottis.stress_press": ("linear", -0.003, 0.005),
    "glottis.final_relax_gap": ("linear", -0.015, 0.01),
    "glottis.vowel_compression": ("linear", 0.2, 1.0),
    "glottis.unstressed_compression": ("linear", 0.2, 1.0),
    "glottis.adduction_lead": ("linear", -0.02, 0.03),
    "glottis.final_devoicing": ("onset", 0.0, 0.06),
    "breath.speech_pressure": ("log", 500.0, 1000.0),
    "breath.voiceless_pressure": ("linear", 0.5, 1.2),
    "breath.unstressed_pressure": ("linear", 0.6, 1.2),
    "breath.phrase_fade_from": ("linear", 0.3, 0.95),
    "release.release_peel_area": ("onset", 0.0, 0.3),
    "articulation.sibilant_groove": ("linear", 0.3, 1.2),
    "articulation.postalveolar_groove": ("linear", 0.2, 1.0),
    "articulation.affricate_closure": ("linear", 0.15, 0.6),
    "articulation.lateral_retract": ("linear", 0.0, 0.5),
    "articulation.w_round_share": ("linear", 0.0, 1.0),
    "articulation.y_raise": ("linear", -0.2, 0.3),
    "articulation.velar_press": ("linear", -0.05, 0.2),
}


def sigmoid(x):
    return 1.0 / (1.0 + math.exp(-x))


def logit(p):
    return math.log(p / (1.0 - p))


class Neuron:
    """One parameter: its name, where it goes (tree path, or vowel and key),
    activation, bounds and the bias that puts its start at z = 0."""

    def __init__(self, name, block, activation, low, high, start):
        self.name = name
        self.block = block
        self.activation = activation
        self.low = low
        self.high = high
        self.start = start
        self.bias = 0.0
        if activation == "log":
            share = (math.log(start) - math.log(low)) / (math.log(high) - math.log(low))
            self.bias = logit(share)
        elif activation == "linear":
            self.bias = logit((start - low) / (high - low))

    def value(self, z):
        if self.activation == "onset":
            return self.high * math.tanh(max(z, 0.0))
        share = sigmoid(self.bias + z)
        if self.activation == "log":
            return math.exp(math.log(self.low) + share * (math.log(self.high) - math.log(self.low)))
        return self.low + share * (self.high - self.low)


# The physical tongue's switches (never neurons).
TONGUE_SWITCHES = ("tongue.physical", "tongue.statics", "tongue.smooth_airway", "tongue.root_physical",
                   "tongue.dorsum_physical", "tongue.blade_physical", "tongue.tip_physical", "tongue.roles",
                   "tongue.shaper_pinned", "tongue.valve_pinned", "tongue.nozzle_pinned", "tongue.role_anchor",
                   "tongue.lateral_from_posture", "tongue.feedback_aim", "tongue.guard",
                   "tongue.solid_hull", "tongue.blade_body", "tongue.spine", "tongue.skeleton",
                   "tongue.skeleton_curve", "tongue.skeleton_minimal_contact", "tongue.smooth_airway",
                   "tongue.spine_base_end", "tongue.spine_tip_start")


def tongue_neuron(path, start):
    """A physical tongue setting: x/÷2 around its start (shares bounded)."""
    if path in ("tongue.tip_pliancy", "tongue.tip_mass", "tongue.tip_stretch"):
        return Neuron(path, "tongue", "log", max(start / 2.0, 0.02), min(start * 2.0, 1.0) if start < 1.0 else start * 1.5, start)
    return Neuron(path, "tongue", "log", start / 2.0, start * 2.0, start)


def tree_neuron(path, start):
    block = path.split(".", 1)[0]
    if path in TREE_RANGES:
        activation, low, high = TREE_RANGES[path]
    elif start > 0.0:
        activation, low, high = "log", start / DEFAULT_SPAN, start * DEFAULT_SPAN
    else:
        raise ValueError("no range for %s (start %g)" % (path, start))
    return Neuron(path, block, activation, low, high, start)


def organ_neuron(vowel, key, start):
    """Vowel postures move within reach of where the shape tuning left them
    (organ_vopt3 held them near MRI area functions; a free mouth flattened
    toward a tube before)."""
    name = "organs.%s.%s" % (vowel, key)
    block = "organs." + vowel
    if key == "jaw":
        return Neuron(name, block, "linear", max(start - 4.0, -0.5), start + 4.0, start)
    if key in ("tx", "ty", "tipx", "tipy"):
        return Neuron(name, block, "linear", start - 0.6, start + 0.6, start)
    if key == "ap":
        return Neuron(name, block, "log", start / 2.0, start * 2.0, start)
    if key == "pro":
        return Neuron(name, block, "linear", -0.1, start + 0.5, start)
    if key == "lw":
        return Neuron(name, block, "linear", max(start - 1.0, 1.2), start + 1.0, start)
    return Neuron(name, block, "linear", -0.2, start + 0.5, start)  # drop


def plan_neuron(path, start):
    """A planner leaf: its explicit range, or its kind's default (positions
    +-0.4 cm, presses +-0.15 cm, jaw ceilings +-2 degrees, protrusions
    -0.1 to +0.4, everything else positive x/÷2)."""
    block, name = path.split(".", 1)
    if path in PLAN_RANGES:
        activation, low, high = PLAN_RANGES[path]
    elif block == "articulation" and (name.endswith("_x") or name.endswith("_y") or name in (
            "rhotic_back", "rhotic_up", "velar_front", "velar_back")):
        activation, low, high = "linear", start - 0.4, start + 0.4
    elif block == "articulation" and name.endswith("press"):
        activation, low, high = "linear", start - 0.15, start + 0.15
    elif block == "articulation" and name.endswith("_jaw"):
        activation, low, high = "linear", max(start - 2.0, -0.5), start + 2.0
    elif block == "articulation" and name.endswith("protrusion"):
        activation, low, high = "linear", -0.1, start + 0.4
    elif start > 0.0:
        activation, low, high = "log", start / 2.0, start * 2.0
    else:
        raise ValueError("no range for %s (start %g)" % (path, start))
    return Neuron("plan." + path, "plan." + block, activation, low, high, start)


def read_organs(path):
    """[(vowel, {key: value})] in file order, and the comment lines above them."""
    postures = []
    header = []
    with open(path) as handle:
        for line in handle:
            body = line.split("#", 1)[0].split()
            if not body:
                if not postures:
                    header.append(line.rstrip("\n"))
                continue
            values = {}
            for field in body[1:]:
                key, number = field.split("=")
                values[key] = float(number)
            postures.append((body[0], values))
    return postures, header


class EngineNet:
    """The neurons of the groups being fitted, in a fixed order; the
    latents are a numpy vector."""

    def __init__(self, groups=("body", "organs")):
        self.groups = tuple(groups)
        self.tree = read_tree(TREE)
        self.postures, self.organ_header = read_organs(ORGANS)
        self.plan = read_tree(PLAN)
        self.neurons = []
        if "body" in self.groups:
            for path, value in self.tree:
                if path not in LEFT_OUT and not path.startswith("tongue."):
                    self.neurons.append(tree_neuron(path, value))
        if "tongue" in self.groups:
            # The physical tongue's own settings (not its switch); with
            # ENGINE_TONGUE_KEYS (comma-separated prefixes, e.g.
            # "tongue.skeleton_,tongue.tip_radius") only those.
            keys = tuple(key for key in os.environ.get("ENGINE_TONGUE_KEYS", "tongue.").split(",") if key)
            for path, value in self.tree:
                if path.startswith(keys) and path not in TONGUE_SWITCHES and value > 0.0:
                    self.neurons.append(tongue_neuron(path, value))
        if "organs" in self.groups:
            for vowel, values in self.postures:
                for key in ORGAN_KEYS:
                    self.neurons.append(organ_neuron(vowel, key, values[key]))
        if "plan" in self.groups:
            for path, value in self.plan:
                if path.split(".", 1)[0] not in PLAN_LEFT_OUT_BLOCKS and path not in PLAN_LEFT_OUT:
                    self.neurons.append(plan_neuron(path, value))
        self.blocks = []
        for neuron in self.neurons:
            if neuron.block not in self.blocks:
                self.blocks.append(neuron.block)

    def size(self):
        return len(self.neurons)

    def block_indices(self, block):
        return [index for index, neuron in enumerate(self.neurons) if neuron.block == block]

    def values(self, latents):
        return {neuron.name: neuron.value(float(z)) for neuron, z in zip(self.neurons, latents)}

    def write(self, latents, prefix):
        """The engine's three parameter files at these latents: PREFIX.params
        (body), PREFIX.organs (vowel postures), PREFIX.plan (planner).
        Returns their paths."""
        values = self.values(latents)
        paths = (prefix + ".params", prefix + ".organs", prefix + ".plan")
        write_tree(paths[0], [[path, values.get(path, value)] for path, value in self.tree],
                   "# engine_net.py: the network's parameter tree")
        lines = ["# engine_net.py: the network's vowel postures (from %s)" % os.path.basename(ORGANS)]
        for vowel, posture in self.postures:
            fields = []
            for key in ORGAN_KEYS:
                fields.append("%s=%.4f" % (key, values.get("organs.%s.%s" % (vowel, key), posture[key])))
            lines.append("%s %s" % (vowel, " ".join(fields)))
        with open(paths[1], "w") as handle:
            handle.write("\n".join(lines) + "\n")
        write_tree(paths[2], [[path, values.get("plan." + path, value)] for path, value in self.plan],
                   "# engine_net.py: the network's planner tree")
        return paths


def render_word(task):
    """task = (workdir, word, (body, organs, plan) paths, seed, tag): plans
    the word with these postures and plan tree, speaks it with this body
    tree, and returns (tag, word, frames or None, pitch straying or None)."""
    workdir, word, (tree_path, organs_path, plan_path), seed, tag = task
    folder = os.path.join(workdir, "renders", tag)
    os.makedirs(folder, exist_ok=True)
    score = os.path.join(folder, word + ".score")
    wave = os.path.join(folder, word + ".wav")
    cycles = os.path.join(folder, word + ".cycles")
    areas = os.path.join(folder, word + ".areas")
    result = (tag, word, None, None)
    planned = subprocess.run([os.path.join(ROOT, "voice", "build", "say"), word + ".", "/dev/null", "--score", score,
                              "--plan-only", "--pressure", "750", "--phones",
                              os.path.join(ROOT, "tract", "data", "phones_tract.txt"), "--organs", organs_path,
                              "--plan-params", plan_path],
                             cwd=os.path.join(ROOT, "voice"), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if planned.returncode == 0:
        environment = dict(os.environ)
        # ENGINE_SOURCE=kinematic: the kinematic folds (VocalTractLab 2025 glottis, 2026-10-18) instead of
        # the fold tube.
        if os.environ.get("ENGINE_SOURCE") == "kinematic":
            environment.pop("FOLDTUBE", None)
            environment["KINEMATIC"] = "1"
        else:
            environment["FOLDTUBE"] = "1"
        environment["CYCLES"] = cycles
        environment["PARAMS"] = tree_path
        environment["AREALOG"] = areas
        finished = subprocess.run([os.path.join(ROOT, "tract", "build", "speak"), score, wave, seed],
                                  env=environment, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if finished.returncode == 0:
            signal, _ = librosa.load(wave, sr=RATE)
            if numpy.all(numpy.isfinite(signal)) and numpy.abs(signal).max() > 0.0:
                frames = fit_blocks.frame_features(signal)
                if frames is not None:
                    frames["ringing"] = vowel_ringing(score, areas)
                result = (tag, word, frames, fit_blocks.pitch_straying(cycles))
    for path in (score, wave, cycles, areas):
        if os.path.exists(path):
            os.remove(path)
    return result


VOWELS = ("AA", "AE", "AH", "AO", "AW", "AY", "EH", "ER", "EY", "IH", "IY", "OW", "OY", "UH", "UW", "AX")


def moving_mean(values, width):
    return numpy.convolve(values, numpy.ones(width) / width, "same")


def vowel_ringing(score_path, areas_path):
    """The tongue's ringing in the word's vowels: the oral nodes' log areas
    (past 7 cm), their 5-30 Hz part (6 ms less 35 ms moving means), rms per
    vowel, the largest vowel's (0 without vowels of 20 ms)."""
    vowels = []
    with open(score_path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 4 and fields[0] == "#" and fields[1].split("-")[0].rstrip("012") in VOWELS:
                vowels.append((float(fields[2]), float(fields[2]) + float(fields[3])))
    rows = [line.split() for line in open(areas_path)]
    if not rows or not vowels:
        return 0.0
    width = min(len(row) for row in rows)
    table = numpy.array([[float(value) for value in row[:width]] for row in rows])
    times = table[:, 0]
    logs = numpy.log(numpy.maximum(table[:, 1 + int(7.0 / 0.3646):], 1.0e-3))
    band = numpy.array([moving_mean(column, 6) - moving_mean(column, 35) for column in logs.T]).T
    largest = 0.0
    for start, finish in vowels:
        inside = (times >= start + 0.01) & (times <= finish - 0.01)
        if inside.sum() > 20:
            largest = max(largest, float(numpy.sqrt(numpy.mean(band[inside] ** 2))))
    return largest


def render_latents(pool, workdir, net, candidates):
    """candidates: [(tag, latents, seed)]. Returns {tag: ([frames per word],
    [straying per word])}."""
    tasks = []
    results = {}
    for tag, latents, seed in candidates:
        paths = net.write(latents, os.path.join(workdir, "net_" + tag))
        results[tag] = ([None] * len(WORDS), [None] * len(WORDS))
        for word in WORDS:
            tasks.append((workdir, word, paths, seed, tag))
    for tag, word, frames, straying in pool.imap_unordered(render_word, tasks):
        index = WORDS.index(word)
        results[tag][0][index] = frames
        results[tag][1][index] = straying
    return results


TIMING_POINTS = 5
# The tongue's ringing inside a word's vowels (tongue_ringing.py: 5-30 Hz
# in the oral log areas): free below RING_LIMIT, then a residual per
# RING_SCALE (the owner heard a warble at 0.11 and 0.12, none at 0.04).
RING_LIMIT = 0.035
RING_SCALE = 0.002
# ENGINE_RINGING=0 drops the hinge: the 5-30 Hz band also measures how much
# the tongue moves inside short vowels, so it rewarded a sluggish tongue
# (NOTES.md, 2026-10-16).
RINGING_ON = os.environ.get("ENGINE_RINGING", "1") != "0"
DURATION_SCALE = 0.1   # log units of a word's length per residual unit
TIMING_SCALE = 0.03    # s of inner timing per residual unit


def timing_residuals(reference, ours, alignment):
    """The word's length (log ratio of the loud spans, / DURATION_SCALE)
    and its inner timing: at TIMING_POINTS evenly placed loud reference
    frames, where ours is against where the reference is, each as a share
    of its own span, times the reference's span (s, / TIMING_SCALE)."""
    first, last = fit_blocks.loud_span(reference, 0)
    our_first, our_last = fit_blocks.loud_span(ours, 0)
    reference_span = max(last - first, 1)
    our_span = max(our_last - our_first, 1)
    loud = fit_blocks.loud_frames(reference)
    rows = [math.log(our_span / reference_span) / DURATION_SCALE]
    for point in range(TIMING_POINTS):
        pick = int((point + 0.5) / TIMING_POINTS * len(loud))
        frame = loud[pick]
        partner = float(numpy.mean(alignment[pick]))
        drift = (partner - our_first) / our_span - (frame - first) / reference_span
        rows.append(drift * reference_span * fit_blocks.HOP / RATE / TIMING_SCALE)
    return numpy.array(rows)


def score(references, rendered, alignments=None):
    """(residual vector, whole, alignments): per word, the frame residuals
    pooled over SEGMENTS runs of the reference's loud frames, then timing;
    then the pitch straying. The alignment is solved here unless given (ALS
    holds it through a block's turn)."""
    frames, straying = rendered
    if alignments is None:
        alignments = fit_blocks.solve_alignment(references, rendered)
    scales = numpy.full(fit_blocks.FEATURES, fit_blocks.FEATURE_SCALE)
    parts = []
    whole = True
    saved = fit_blocks.SEGMENTS
    fit_blocks.SEGMENTS = SEGMENTS
    for reference, ours, alignment in zip(references, frames, alignments):
        if reference is None:
            continue
        if ours is None or alignment is None:
            whole = False
            parts.append(numpy.full(SEGMENTS * fit_blocks.FEATURES + 2 + TIMING_POINTS, 10.0))
            continue
        parts.append(fit_blocks.frame_residuals(reference, ours, alignment, scales))
        parts.append(timing_residuals(reference, ours, alignment))
        parts.append(numpy.array([max(ours.get("ringing", 0.0) - RING_LIMIT, 0.0) / RING_SCALE if RINGING_ON else 0.0]))
    fit_blocks.SEGMENTS = saved
    parts.append(numpy.array([10.0 if value is None else value / fit_blocks.PITCH_SCALE for value in straying]))
    whole = whole and all(value is not None for value in straying)
    return numpy.concatenate(parts), whole, alignments


def residual_mask(view, length):
    """Which entries of score()'s vector a view sees: "timing" (each word's
    length and inner timing), "shape" (the frame residuals and the pitch
    straying) or "all"."""
    mask = numpy.ones(length, dtype=bool)
    if view == "all":
        return mask
    per_word = SEGMENTS * fit_blocks.FEATURES + 2 + TIMING_POINTS
    words = (length - len(WORDS)) // per_word
    timing = numpy.zeros(length, dtype=bool)
    for word in range(words):
        timing[word * per_word + SEGMENTS * fit_blocks.FEATURES:(word + 1) * per_word - 1] = True
    return timing if view == "timing" else ~timing


def groups_of(vector):
    """Mean squared residual by kind: spectral shape, level, periodicity,
    timing, pitch."""
    per_word = SEGMENTS * fit_blocks.FEATURES + 2 + TIMING_POINTS
    body = vector[:len(vector) - len(WORDS)].reshape(-1, per_word)
    frames = body[:, :SEGMENTS * fit_blocks.FEATURES].reshape(-1, fit_blocks.FEATURES)
    return {"shape": float(numpy.mean(frames[:, :fit_blocks.LEVEL_BANDS] ** 2)),
            "level": float(numpy.mean(frames[:, fit_blocks.TOTAL] ** 2)),
            "periodicity": float(numpy.mean(frames[:, fit_blocks.TOTAL + 1:] ** 2)),
            "duration": float(numpy.mean(body[:, SEGMENTS * fit_blocks.FEATURES] ** 2)),
            "timing": float(numpy.mean(body[:, SEGMENTS * fit_blocks.FEATURES + 1:-1] ** 2)),
            "ringing": float(numpy.mean(body[:, -1] ** 2)),
            "pitch": float(numpy.mean(vector[len(vector) - len(WORDS):] ** 2))}


def behaviour_cost(vector, whole):
    """The reference mismatch alone (the optimizer adds its own decay)."""
    return fit_blocks.cost_of(vector, whole, {})
