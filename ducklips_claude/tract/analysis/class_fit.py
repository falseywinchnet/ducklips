"""Separable fit of class-specific settings: the vowels' organ postures and
the consonants' plan shapes, each credited by its own phone classes' share
of the error.

The combined fidelity is 10 log10(groups / sum over groups of error/spread),
and spectrum, dynamics and timing errors are sums over phone classes, so a
class's cost (its squared errors over each group's spread) adds up to that
objective. A setting moves only its own classes' costs, nearly: one render
of a candidate perturbing every setting at once ranks every setting by its
own classes. (profile_fit's single combined score ranked 66 settings by one
number whose phase group swung 5-11 dB with the render seed alone.)

Each generation: antithetic pairs around the mean, all rendered on the
training sentences with one seed (a new seed each generation, so no setting
fits one noise draw); per setting, the best four candidates by its classes'
cost recombine into its new mean. The mean is scored on the held-out
sentences every few generations.
usage: class_fit.py NAME REFDIR TRAIN HELDOUT PARAMS ORGANS PLAN GENERATIONS [SETS]
REFDIR "nato": the 26 NATO words in one render against the reader's clips
(even words train, odd held out; TRAIN and HELDOUT are then ignored).
SETS: vowels,consonants,tree (default all three). Writes
/tmp/es/class/NAME/{mean,best}_{organs.txt,plan.params,voice.params}."""
import json
import multiprocessing
import os
import pickle
import re
import shutil
import sys
import tempfile

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fidelity
import render_set
import slider_server

# Organ posture fields: half-range about the start (cm, degrees) and floor.
ORGAN_FIELDS = {"jaw": (4.0, 0.0), "tx": (0.6, None), "ty": (0.4, None), "tipx": (0.4, None), "tipy": (0.4, None),
                "ap": (0.5, 0.02), "pro": (0.4, 0.0), "lw": (1.0, 0.5), "drop": (0.5, 0.0)}
# Each posture's classes and their weights (diphthongs start or end on it;
# planner.cpp vowel_organs).
VOWEL_CLASSES = {"AA": {"aa": 1, "ay": 0.5, "aw": 0.5}, "AE": {"ae": 1}, "AH": {"ah": 1},
                 "AO": {"ao": 1, "oy": 0.5}, "EH": {"eh": 1}, "EY": {"ey": 0.6},
                 "IH": {"ih": 1, "ay": 0.5, "oy": 0.5, "ey": 0.4}, "IY": {"iy": 1, "glide_palatal": 0.5},
                 "OW": {"ow": 0.6}, "UH": {"uh": 1, "aw": 0.5, "ow": 0.4}, "UW": {"uw": 1, "glide_labial": 0.5},
                 "ER": {"er": 1, "liquid_rhotic": 0.5}, "AX": {"ax": 1}}
# Plan shapes: (block.key, low, high, classes).
SIBILANT = {"fric_alveolar": 1, "fric_alveolar_voiced": 1}
VOICELESS_SOUNDS = {"stop_voiceless": 1, "fric_labial": 1, "fric_dental": 1, "fric_alveolar": 1, "fric_postalveolar": 1,
                    "affricate": 1, "fric_glottal": 1}
POSTALVEOLAR = {"fric_postalveolar": 1, "affricate": 1, "affricate_voiced": 1}
CONSONANTS = [
    ("articulation.labiodental_aperture", 0.005, 0.3, {"fric_labial": 1, "fric_labial_voiced": 1}),
    ("articulation.labiodental_jaw", 0.0, 10.0, {"fric_labial": 1, "fric_labial_voiced": 1}),
    ("articulation.dental_press", -0.2, 0.2, {"fric_dental": 1}),
    # /dh/ its own press (bdl's moves 11-14.5 dB within its 55 ms, a near
    # stop; ours 7-12).
    ("articulation.dh_press", -0.2, 0.3, {"fric_dental_voiced": 1}),
    ("articulation.dental_jaw", 0.0, 10.0, {"fric_dental": 1, "fric_dental_voiced": 1}),
    ("articulation.sibilant_press", -0.3, 0.1, SIBILANT),
    ("articulation.sibilant_x", 6.6, 7.5, SIBILANT),
    ("articulation.sibilant_groove", 0.0, 2.0, SIBILANT),
    ("articulation.sibilant_aperture", 0.1, 1.5, SIBILANT),
    ("articulation.postalveolar_hiss_press", -0.4, 0.1, POSTALVEOLAR),
    ("articulation.postalveolar_x", 5.8, 6.8, POSTALVEOLAR),
    ("articulation.postalveolar_groove", 0.0, 2.0, POSTALVEOLAR),
    ("articulation.postalveolar_aperture", 0.1, 1.5, POSTALVEOLAR),
    ("articulation.postalveolar_protrusion", 0.0, 1.0, POSTALVEOLAR),
    ("articulation.alveolar_press", 0.0, 0.4, {"stop_voiceless": 0.4, "stop_voiced": 0.4, "nasal_alveolar": 1}),
    ("articulation.velar_press", 0.0, 0.3, {"stop_voiceless": 0.3, "stop_voiced": 0.3, "nasal_velar": 1}),
    ("articulation.lateral_press", -0.1, 0.4, {"liquid_lateral": 1}),
    ("articulation.lateral_retract", 0.0, 0.8, {"liquid_lateral": 1}),
    ("articulation.rhotic_press", -0.8, 0.0, {"liquid_rhotic": 1}),
    ("articulation.rhotic_back", 0.0, 1.2, {"liquid_rhotic": 1}),
    ("articulation.rhotic_up", -0.3, 0.6, {"liquid_rhotic": 1}),
    ("articulation.rhotic_aperture", 0.1, 1.5, {"liquid_rhotic": 1}),
    ("articulation.rhotic_protrusion", 0.0, 1.0, {"liquid_rhotic": 1}),
    ("articulation.y_raise", -0.4, 0.3, {"glide_palatal": 1}),
    ("articulation.y_round_share", 0.0, 1.0, {"glide_palatal": 1}),
    ("articulation.w_aperture", 0.05, 1.2, {"glide_labial": 1}),
    ("articulation.w_round_share", 0.0, 1.0, {"glide_labial": 1}),
    # Stop timing: bdl's /p t k/ hold the closure to about 70% of the token
    # and their aspiration carries almost nothing under 500 Hz; ours release
    # by 60% (a "keep" /k/ sealed 33 of 91 ms) and the release flow rises in
    # the low bands (bdl4).
    ("timing.burst_scale", 0.2, 2.0, {"stop_voiceless": 1, "stop_voiced": 1}),
    ("timing.seal_lead", 0.0, 1.5, {"stop_voiceless": 1, "stop_voiced": 1, "nasal_labial": 0.3, "nasal_alveolar": 0.3,
                                    "nasal_velar": 0.3, "affricate": 0.5, "affricate_voiced": 0.5}),
    ("glottis.aspirate_gap", 0.02, 0.4, {"stop_voiceless": 1}),
    # The folds part for voiceless sounds (bdl's /s/ loses its voicing band
    # within the first quarter of the token; ours through half).
    ("glottis.abduction_move", 0.01, 0.1, VOICELESS_SOUNDS),
    ("glottis.fricative_abduct_lead", 0.0, 0.08, VOICELESS_SOUNDS),
    ("articulation.dorsum_close_move", 0.02, 0.12, {"stop_voiceless": 0.3, "stop_voiced": 0.3, "nasal_velar": 1}),
    ("articulation.tip_close_move", 0.01, 0.06, {"stop_voiceless": 0.4, "stop_voiced": 0.4, "nasal_alveolar": 1}),
    # Affricates: bdl's /ch/ closes for about 40% and hisses strongly after;
    # ours closed to 75% (the hiss capped at affricate_hiss_max, 60 ms) and
    # hissed 10 dB weaker; his /dZ/ devoices through its hiss, ours re-voiced.
    ("timing.affricate_hiss_max", 0.03, 0.2, {"affricate": 1, "affricate_voiced": 1}),
    ("articulation.affricate_closure", 0.15, 0.7, {"affricate": 1, "affricate_voiced": 1}),
    ("glottis.affricate_hiss_gap", 0.0, 0.25, {"affricate_voiced": 1}),
    # bdl's /z/ before no vowel devoices through the token (under 400 Hz -5
    # -> -34 dB re speech), ours held its voicing (-19).
    ("glottis.final_devoicing", 0.0, 0.3, {"fric_alveolar_voiced": 1, "fric_labial_voiced": 1}),
    # rms's /z/ hisses like his /s/ (-5 to -10 dB re speech at 3.4-5.8 kHz);
    # ours, the folds nearly shut at 0.04, 21-30 dB short up there (rms6).
    ("glottis.voiced_fricative_gap", 0.0, 0.25, {"fric_alveolar_voiced": 1, "fric_labial_voiced": 1,
                                                 "fric_postalveolar": 0.3, "affricate_voiced": 0.5}),
    ("glottis.voiceless_gap", 0.05, 0.6, {"stop_voiceless": 1, "fric_labial": 1, "fric_dental": 1, "fric_alveolar": 1,
                                          "fric_postalveolar": 1, "affricate": 1}),
]
# Engine tree settings with class credit: (block.key, low, high, classes).
FRICATION = {"fric_labial": 1, "fric_labial_voiced": 1, "fric_dental": 1, "fric_dental_voiced": 1, "fric_alveolar": 1,
             "fric_alveolar_voiced": 1, "fric_postalveolar": 1, "affricate": 1, "affricate_voiced": 1,
             "stop_voiceless": 0.5}
VOICELESS = {"stop_voiceless": 1, "fric_labial": 1, "fric_dental": 1, "fric_alveolar": 1, "fric_postalveolar": 1,
             "affricate": 1, "fric_glottal": 1}
VOICED_CLOSURE = {"stop_voiced": 1, "affricate_voiced": 0.5, "fric_labial_voiced": 0.5, "fric_dental_voiced": 0.5,
                  "fric_alveolar_voiced": 0.5}
# (block.key, low, high, default when the voice lacks it, classes)
TREE = [
    ("oral_noise.labiodental", 0.05, 3.0, 1.2, {"fric_labial": 1, "fric_labial_voiced": 1}),
    ("oral_noise.frication", 0.3, 5.0, 0.5, FRICATION),
    ("oral_noise.obstacle_reach", 0.2, 1.6, 0.8, {"fric_alveolar": 1, "fric_alveolar_voiced": 1, "fric_postalveolar": 1,
                                                  "fric_dental": 1, "affricate": 1}),
    # (bdl's /m/ against profile_bdl_s2: +8 to +16 dB at 0.4-3 kHz.)
    ("branches.sinus_neck_scale", 0.05, 1.5, 0.25, {"nasal_labial": 1, "nasal_alveolar": 1, "nasal_velar": 1}),
    # The folds kept swinging through open /f/ (0.17-0.22 cm^2 at f0, 2200
    # dPa across them): voicing at 100-500 Hz in voiceless sounds.
    ("kinematic.abduction_fade_start", 0.0, 0.15, 0.0, VOICELESS),
    ("kinematic.abduction_fade_end", 0.0, 0.3, 0.0, VOICELESS),
    # (bdl4 round 1: /hh/ +13/+6/+7/+15 dB by quarters of the band.)
    ("vtl.aspiration_gain", 0.5, 30.0, 1.0, {"fric_glottal": 1, "stop_voiceless": 0.5}),
    # Voiced closures devoice: bdl's voice bar sits at -9 to -16 dB under
    # 500 Hz through /b d g/, ours fell to -29 (bdl4 round 1).
    ("walls.closure_pull", 0.0, 1.0, 0.25, VOICED_CLOSURE),
    ("walls.closure_relaxation", 0.0, 0.9, 0.0, VOICED_CLOSURE),
    ("posture.closure_press", 0.0, 0.02, 0.004, VOICED_CLOSURE),
    # The voice bar radiates through the walls: the folds swing all through a
    # /b/ closure, but at skin gain 0.037 its 0-500 Hz level sat at -28 dB re
    # speech (bdl -9 to -16); at 1.0, -12 (one sentence) though the whole
    # held-out set scored 6.9 -> 6.7 (voiced affricates +3.3 dB).
    ("walls.skin_gain", 0.0, 2.0, 1.0, dict(VOICED_CLOSURE, stop_voiceless=0.5)),
    ("walls.skin_corner", 100.0, 600.0, 250.0, VOICED_CLOSURE),
    # (Nasal murmurs too bright: bdl's vowel-to-murmur drop at 0.5-2 kHz is
    # 10 dB, ours 5, bdl4.)
    ("branches.nose_smoothing", 0.05, 0.8, 0.4, {"nasal_labial": 1, "nasal_alveolar": 1, "nasal_velar": 1}),
    # (Our /l/ sat 7-19 dB short at 0.25-1.3 kHz for bdl, rms and slt.)
    ("anatomy.lateral_area", 0.06, 0.5, 0.12, {"liquid_lateral": 1}),
    # The /l/ side branch over the blade (its zero near 2 kHz; off at 0).
    ("branches.lateral_pocket", 0.0, 0.8, 0.0, {"liquid_lateral": 1}),
    ("branches.lateral_pocket_length", 2.0, 6.0, 4.0, {"liquid_lateral": 1}),
]
GROUPS = ("spectrum", "dynamics", "timing")
WORK = "/tmp/es/class"
HELDOUT_SEEDS = (render_set.SEED, "11")
NATO_TEXT = ("alfa. bravo. charlie. delta. echo. foxtrot. golf. hotel. india. juliett. kilo. lima. mike. november. "
             "oscar. papa. quebec. romeo. sierra. tango. uniform. victor. whiskey. x-ray. yankee. zulu.")
NATO_DIR = os.path.expanduser("~/ducklips/voice/ref/nato")
STATE = {}


def setup(state):
    STATE.update(state)


def organ_values(text):
    out = {}
    for line in text.split("\n"):
        name = line.split(" ", 1)[0]
        if name in VOWEL_CLASSES:
            for field in ORGAN_FIELDS:
                found = re.search(r"\b%s=([-0-9.e]+)" % field, line)
                if found:
                    out["organs.%s.%s" % (name, field)] = float(found.group(1))
    return out


def organs_with(text, values):
    lines = text.split("\n")
    for index, line in enumerate(lines):
        name = line.split(" ", 1)[0]
        if name not in VOWEL_CLASSES:
            continue
        for field in ORGAN_FIELDS:
            key = "organs.%s.%s" % (name, field)
            if key in values:
                line = re.sub(r"\b%s=[-0-9.e]+" % field, "%s=%.4f" % (field, values[key]), line)
        lines[index] = line
    return "\n".join(lines)


def settings_of(organs_text, plan_text, params_text, sets):
    """[(key, low, high, {class: weight})] with bounds about the start."""
    out = []
    if "tree" in sets:
        for key, low, high, default, classes in TREE:
            out.append(("tree." + key, low, high, classes))
    if "vowels" in sets:
        for key, value in sorted(organ_values(organs_text).items()):
            _, name, field = key.split(".")
            half, floor = ORGAN_FIELDS[field]
            low = value - half if floor is None else max(value - half, floor)
            out.append((key, low, value + half, VOWEL_CLASSES[name]))
    if "consonants" in sets:
        for key, low, high, classes in CONSONANTS:
            if slider_server.read_value(plan_text, key) is not None or key in ("glottis.voiceless_gap",
                                                                                "glottis.affricate_hiss_gap",
                                                                                "articulation.dh_press"):
                out.append((key, low, high, classes))
    return out


def start_values(settings, organs_text, plan_text, params_text):
    organs = organ_values(organs_text)
    out = {}
    for key, low, high, classes in settings:
        if key.startswith("organs."):
            value = organs[key]
        elif key.startswith("tree."):
            value = slider_server.read_value(params_text, key[5:])
            if value is None:
                value = next(t[3] for t in TREE if t[0] == key[5:])
        else:
            found = slider_server.read_value(plan_text, key)
            if found is None and key == "articulation.dh_press":
                found = slider_server.read_value(plan_text, "articulation.dental_press")
            if found is None and key == "glottis.affricate_hiss_gap":
                # (Absent, the planner uses the voicing gap.)
                found = slider_server.read_value(plan_text, "glottis.voicing_gap")
            value = found if found is not None else (0.14 if key == "glottis.voiceless_gap" else low)
        out[key] = min(max(value, low), high)
    return out


def write_files(folder, values):
    organs = os.path.join(folder, "organs.txt")
    plan = os.path.join(folder, "plan.params")
    params = os.path.join(folder, "voice.params")
    open(organs, "w").write(organs_with(STATE["organs_text"], values))
    plan_values = {k: v for k, v in values.items() if not k.startswith("organs.") and not k.startswith("tree.")}
    open(plan, "w").write(slider_server.with_values(STATE["plan_text"], plan_values) if plan_values else STATE["plan_text"])
    tree_values = {k[5:]: v for k, v in values.items() if k.startswith("tree.")}
    open(params, "w").write(slider_server.with_values(STATE["params_text"], tree_values) if tree_values else STATE["params_text"])
    return organs, plan, params


def describe_nato(values, half, seed):
    """NATO mode: the 26 words as one render, the reader's clips cut by our
    phones (DTW); `half` 0 is the even words (training), 1 the odd ones.
    Returns (reference description, ours)."""
    folder = tempfile.mkdtemp(dir=STATE["work"])
    organs, plan, params = write_files(folder, values)
    render_set.SEED = seed
    done = render_set.render(("nato", NATO_TEXT, folder, params, organs, plan, {"KINEMATIC": "1"}))
    if done[1] != "ok":
        shutil.rmtree(folder, ignore_errors=True)
        return None
    reference, ours = fidelity.nato_sets(os.path.join(folder, "nato.wav"), os.path.join(folder, "nato.score"), NATO_DIR)
    shutil.rmtree(folder, ignore_errors=True)
    return fidelity.describe(reference[half::2]), fidelity.describe(ours[half::2])


def describe_values(job):
    values, ids, seed = job[:3]
    threads = job[3] if len(job) > 3 else 1
    if STATE.get("nato"):
        return describe_nato(values, ids, seed)
    folder = tempfile.mkdtemp(dir=STATE["work"])
    organs, plan, params = write_files(folder, values)
    table = render_set.texts()
    render_set.SEED = seed
    jobs = [(name, table[name], folder, params, organs, plan, {"KINEMATIC": "1"}) for name in ids]
    if any(result[1] != "ok" for result in render_set.render_many(jobs, threads)):
        shutil.rmtree(folder, ignore_errors=True)
        return None
    described = fidelity.describe(fidelity.ours_set(folder, ids, STATE.get("ref_dir")))
    shutil.rmtree(folder, ignore_errors=True)
    return described


def spreads(reference):
    out = {}
    for group in GROUPS:
        r = numpy.array(list(reference[group].values()))
        out[group] = float(numpy.sum((r - r.mean()) ** 2))
    return out


def class_costs(reference, ours, spread):
    """{class: its squared errors over each group's spread}."""
    out = {}
    for group in GROUPS:
        for key, r in reference[group].items():
            if key in ours[group]:
                c = key.rsplit("_", 1)[0]
                out[c] = out.get(c, 0.0) + (ours[group][key] - r) ** 2 / spread[group]
    return out


def credit(costs, classes, missing):
    return sum(weight * costs.get(c, missing) for c, weight in classes.items())


def main():
    name, ref_dir, train_ids, heldout_ids, params, organs_path, plan_path, generations = sys.argv[1:9]
    sets = set(sys.argv[9].split(",")) if len(sys.argv) > 9 else {"vowels", "consonants", "tree"}
    work = os.path.join(WORK, name)
    os.makedirs(work, exist_ok=True)
    cache = os.path.join(work, "reference.pkl")
    nato = ref_dir == "nato"
    if nato:
        reference_train, reference_heldout = None, None
    elif os.path.exists(cache):
        reference_train, reference_heldout = pickle.load(open(cache, "rb"))
    else:
        reference_train = fidelity.describe(fidelity.arctic_set(ref_dir, fidelity.id_range(train_ids)))
        reference_heldout = fidelity.describe(fidelity.arctic_set(ref_dir, fidelity.id_range(heldout_ids)))
        pickle.dump((reference_train, reference_heldout), open(cache, "wb"))
    organs_text, plan_text, params_text = open(organs_path).read(), open(plan_path).read(), open(params).read()
    settings = settings_of(organs_text, plan_text, params_text, sets)
    keys = [s[0] for s in settings]
    low = numpy.array([s[1] for s in settings])
    high = numpy.array([s[2] for s in settings])
    mean = numpy.array([start_values(settings, organs_text, plan_text, params_text)[k] for k in keys])
    sigma = 0.15 * (high - low)
    spread = spreads(reference_train) if not nato else None
    state = {"work": work, "params_text": params_text, "organs_text": organs_text, "plan_text": plan_text, "nato": nato,
             "ref_dir": None if nato else ref_dir}
    setup(state)
    # (NATO: the "ids" are the halves of the 26 words.)
    train, heldout = (0, 1) if nato else (fidelity.id_range(train_ids), fidelity.id_range(heldout_ids))
    pairs, parents = 5, 4
    weights = numpy.log(parents + 0.5) - numpy.log(numpy.arange(1, parents + 1))
    weights /= weights.sum()
    rng = numpy.random.default_rng(abs(hash(name)) % (2 ** 32))
    log = open(os.path.join(work, "log.jsonl"), "a")
    best_heldout = -100.0
    print("%d settings (%s)" % (len(keys), ",".join(sorted(sets))), flush=True)
    with multiprocessing.Pool(2 * pairs, initializer=setup, initargs=(state,)) as pool:
        for generation in range(int(generations) + 1):
            values = dict(zip(keys, mean))
            if generation % 3 == 0 or generation == int(generations):
                # Two render seeds, each group's dB averaged (one seed's phase
                # group alone ranged 11.3-16.0 dB over a fixed voice).
                # (Each seed's renders 5 at a time: the other workers wait.)
                helds = pool.map(describe_values, [(values, heldout, seed, 5) for seed in HELDOUT_SEEDS])
                scores = []
                for held in helds:
                    if held is None:
                        scores.append({"combined": -100.0})
                    elif nato:
                        scores.append(fidelity.score(held[0], held[1]))
                    else:
                        scores.append(fidelity.score(reference_heldout, held))
                score = {k: float(numpy.mean([one.get(k, -100.0) for one in scores])) for k in scores[0]}
                print("generation %2d mean heldout %s" % (generation, fidelity.report(score)), flush=True)
                write_files(work, values)
                for suffix in ("organs.txt", "plan.params", "voice.params"):
                    shutil.copy(os.path.join(work, suffix), os.path.join(work, "mean_" + suffix))
                if score["combined"] > best_heldout:
                    best_heldout = score["combined"]
                    for suffix in ("organs.txt", "plan.params", "voice.params"):
                        shutil.copy(os.path.join(work, suffix), os.path.join(work, "best_" + suffix))
                log.write(json.dumps({"generation": generation, "heldout": score, "values": values}) + "\n")
                log.flush()
            if generation == int(generations):
                break
            noise = rng.standard_normal((pairs, len(keys)))
            candidates = numpy.clip(numpy.concatenate([mean + sigma * noise, mean - sigma * noise]), low, high)
            seed = str(1000 + generation)
            described = pool.map(describe_values, [(dict(zip(keys, c)), train, seed) for c in candidates])
            if nato:
                costs = [class_costs(d[0], d[1], spreads(d[0])) if d is not None else {} for d in described]
            else:
                costs = [class_costs(reference_train, d, spread) if d is not None else {} for d in described]
            total = [sum(c.values()) if c else float("inf") for c in costs]
            for i, (key, _, _, classes) in enumerate(settings):
                scores = numpy.array([credit(c, classes, 1e3) if c else float("inf") for c in costs])
                order = numpy.argsort(scores)[:parents]
                chosen = candidates[order, i]
                step = numpy.sum(weights * chosen) - mean[i]
                mean[i] = min(max(mean[i] + step, low[i]), high[i])
                # The step size follows the selected spread, within limits.
                selected = numpy.sqrt(numpy.sum(weights * (chosen - mean[i]) ** 2))
                sigma[i] = float(numpy.clip(0.7 * sigma[i] + 0.3 * 1.2 * selected, 0.02 * (high[i] - low[i]),
                                            0.25 * (high[i] - low[i])))
            print("generation %2d train class cost total: best %.3f median %.3f" % (
                generation, min(total), float(numpy.median(total))), flush=True)


if __name__ == "__main__":
    main()
