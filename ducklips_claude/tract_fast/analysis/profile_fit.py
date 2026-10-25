"""Fit a voice profile to one reference speaker by the fidelity score
(fidelity.py): a diagonal evolution strategy over chosen leaves of the
engine's parameter tree, nine candidates a generation rendered in parallel,
each saying the training sentences; the reference's features computed once.
Every generation's best is also scored on held-out sentences.

A profile is a parameter tree (speak PARAMS), organ postures and plan
parameters; this fits the tree's leaves (LEAVES below, or a group of them).
usage: profile_fit.py NAME REFDIR TRAIN_IDS HELDOUT_IDS PARAMS ORGANS PLAN GENERATIONS [GROUPS]
REFDIR "nato": the NATO reader's 26 words (TRAIN_IDS and HELDOUT_IDS unused, pass -).
GROUPS: comma-separated groups of LEAVES (default all)."""
import json
import multiprocessing
import os
import pickle
import shutil
import subprocess
import sys
import tempfile

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fidelity
import render_set
import slider_server

# (group, leaf, low, high, log scale, engine default when the tree leaves it out)
LEAVES = [
    ("source", "kinematic.amplitude", 0.0001, 0.0015, True, 0.0007),  # (0.0003 was a wall, bdl r1)
    ("source", "kinematic.phase_lag", 0.3, 2.0, False, 1.22),
    ("source", "kinematic.pulse_shape", -1.0, 1.0, False, 0.5),
    ("source", "kinematic.flutter", 0.0, 120.0, False, 0.0),  # (60 was a wall, bdl5 r1)
    ("source", "kinematic.jitter", 0.0, 1.5, False, 0.17),
    ("source", "posture.voicing_lower", 0.0, 0.08, False, 0.011),
    ("source", "posture.voicing_upper", 0.0, 0.10, False, 0.006),
    ("source", "posture.voicing_chink", 0.0, 0.05, False, 0.0002),
    ("source", "kinematic.abduction_fade_start", 0.0, 0.1, False, 0.0),
    ("source", "kinematic.abduction_fade_end", 0.0, 0.2, False, 0.0),
    ("noise", "vtl.aspiration_gain", 0.2, 30.0, True, 1.0),
    ("noise", "oral_noise.frication", 0.1, 3.0, True, 2.0),
    ("noise", "oral_noise.jet_corner_ceiling", 3000.0, 12000.0, True, 12000.0),
    ("variation", "variation.breath", 0.0, 0.15, False, 0.03),
    ("variation", "variation.tension", 0.0, 0.05, False, 0.015),
    ("variation", "variation.tremor", 0.0, 0.03, False, 0.0),
    ("variation", "variation.breath_tremor", 0.0, 0.3, False, 0.0),
    ("tract", "walls.loss_scale", 0.05, 3.0, True, 1.0),
    ("tract", "tract.viscous_scale", 0.05, 5.0, True, 1.0),  # (2 was a wall, bdl r1)
    ("tract", "tract.radiation_mass_scale", 0.2, 1.5, False, 1.0),  # (0.5 was a wall, bdl r1)
    ("tract", "tract.radiation_resistance_scale", 0.5, 1.5, False, 1.0),
    ("tract", "walls.skin_gain", 0.0, 3.0, False, 1.0),
    ("tract", "tract.hop_loss", 0.0, 0.15, False, 0.02),
    ("tongue", "tongue.vtl_body_radius", 1.4, 2.2, False, 1.8),
    ("tongue", "tongue.vtl_blade_along", 0.4, 0.9, False, 0.61),
    ("tongue", "tongue.vtl_blade_lift", 0.4, 2.0, False, 1.15),
    ("tongue", "tongue.vtl_root_bulge", -0.3, 1.0, False, 0.3),
    ("cycle", "kinematic.cycle_shape", 0.0, 0.8, False, 0.0),
    # (bdl's phase scatter wanted lag 0.6, shape 0.35: source_calibrate.py.)
    ("cycle", "kinematic.cycle_lag", 0.0, 1.2, False, 0.0),
    ("cycle", "kinematic.cycle_amplitude", 0.0, 0.2, False, 0.0),
    ("cycle", "vtl.aspiration_corner", 500.0, 6000.0, True, 1000.0),
    # The planner's (plan: leaves go to the plan parameters, not the tree).
    ("plan_pitch", "plan:pitch.base_pitch", 70.0, 260.0, True, 120.0),
    ("plan_pitch", "plan:pitch.pitch_median", 70.0, 260.0, True, 122.0),
    ("plan_pitch", "plan:pitch.pitch_floor", 50.0, 220.0, True, 100.0),
    ("plan_pitch", "plan:pitch.lowest_pitch", 40.0, 200.0, True, 80.0),
    ("plan_pitch", "plan:pitch.accent", 5.0, 60.0, True, 22.0),
    ("plan_pitch", "plan:pitch.pitch_spread", 0.5, 3.0, False, 1.5),
    ("plan_pitch", "plan:pitch.final_fall", 0.7, 1.0, False, 0.9),
    ("plan_pitch", "plan:pitch.unstressed_dip", 0.0, 1.0, False, 0.5),
    ("plan_timing", "plan:timing.rate", 0.7, 1.4, False, 1.0),
    ("plan_timing", "plan:timing.final_vowel_scale", 1.0, 2.6, False, 1.8),
    ("plan_timing", "plan:timing.final_consonant_scale", 1.0, 2.4, False, 1.6),
    ("plan_timing", "plan:timing.vowel_hold", 0.3, 0.9, False, 0.6),
    ("plan_timing", "plan:timing.voice_lag", 0.0, 0.05, False, 0.02),
    ("plan_timing", "plan:breath.speech_pressure", 500.0, 1100.0, False, 750.0),
    ("plan_consonants", "plan:articulation.sibilant_x", 6.6, 7.5, False, 7.05),
    ("plan_consonants", "plan:articulation.sibilant_y", 8.7, 9.6, False, 9.15),
    ("plan_consonants", "plan:articulation.sibilant_press", -0.2, 0.05, False, -0.08),
    ("plan_consonants", "plan:articulation.postalveolar_x", 5.9, 6.8, False, 6.3),
    ("plan_consonants", "plan:articulation.postalveolar_hiss_press", -0.25, 0.0, False, -0.12),
    ("plan_consonants", "plan:articulation.alveolar_press", 0.0, 0.3, False, 0.12),
    ("plan_consonants", "plan:articulation.lateral_press", -0.1, 0.3, False, 0.1),
    ("plan_consonants", "plan:articulation.rhotic_press", -0.6, 0.0, False, -0.35),
    ("plan_consonants", "plan:articulation.rhotic_back", 0.0, 1.2, False, 0.6),
    ("plan_consonants", "plan:articulation.velar_press", -0.05, 0.2, False, 0.05),
    ("plan_consonants", "plan:articulation.labiodental_aperture", 0.0, 0.1, False, 0.03),
    # When the constrictions form (measured 2026-10-18: alfa's /f/ tightest only
    # halfway through, so its middle still sounded like the vowel: +28 dB low).
    ("plan_consonants", "plan:articulation.lead", 0.0, 0.12, False, 0.04),
    ("plan_consonants", "plan:articulation.tip_close_move", 0.01, 0.06, False, 0.028),
    ("plan_consonants", "plan:articulation.dorsum_close_move", 0.02, 0.1, False, 0.06),
    ("plan_consonants", "plan:timing.seal_lead", 0.0, 1.5, False, 0.4),
    # The glottis's own planning (measured 2026-10-18: through a /t/ closure the
    # folds stayed within 0.03-0.16 cm^2 and kept vibrating).
    ("plan_glottis", "plan:glottis.voiceless_gap", 0.08, 0.4, False, 0.14),
    ("plan_glottis", "plan:glottis.aspirate_gap", 0.02, 0.2, False, 0.06),
    ("plan_glottis", "plan:glottis.voiced_fricative_gap", 0.0, 0.1, False, 0.04),
    ("plan_glottis", "plan:glottis.fricative_abduct_lead", 0.0, 0.08, False, 0.03),
    ("plan_glottis", "plan:glottis.abduction_move", 0.015, 0.1, False, 0.05),
    ("plan_glottis", "plan:glottis.adduction_move", 0.01, 0.06, False, 0.02),
    ("plan_glottis", "plan:glottis.adduction_lead", -0.02, 0.04, False, 0.0),
    # Each syllable's breath (bdl's vowels move 3.6-6.8 dB within a token, ours
    # 0.5-2.4 at rise 0.35, decay 0.88; unstressed syllables were as loud).
    # How fast the articulators move: bdl's tokens change 1.5-2x more within
    # themselves (octave-level spread), while our vowels' loudness envelope is
    # already more peaked than his (2026-10-18): spectral change, not breath.
    ("plan_motion", "plan:timing.vowel_hold", 0.3, 0.9, False, 0.6),
    ("plan_motion", "plan:timing.release_move", 0.01, 0.08, False, 0.03),
    ("plan_motion", "plan:timing.jaw_lead", 0.0, 0.1, False, 0.04),
    ("plan_motion", "plan:articulation.tip_release_move", 0.01, 0.08, False, 0.03),
    ("plan_motion", "plan:articulation.dorsum_release_move", 0.02, 0.12, False, 0.05),
    ("plan_motion", "plan:articulation.liquid_move", 0.03, 0.15, False, 0.07),
    ("plan_motion", "plan:glottis.adduction_move", 0.01, 0.08, False, 0.02),
    # Every articulator's move time at once (bdl's transitions are sharp at
    # the phone boundaries; ours smeared through half a token).
    ("plan_motion", "motors.move_scale", 0.3, 1.5, False, 1.0),
    # The velum around nasals (bdl's vowel after a nasal is 16 dB over the
    # murmur at 0.5-2 kHz within 10-15 ms; ours rose 8-10 dB, nasalized).
    ("plan_nasal", "plan:nasal.velum", 0.2, 0.7, False, 0.45),
    ("plan_nasal", "plan:nasal.close_lead", -0.03, 0.05, False, 0.01),
    ("plan_nasal", "plan:nasal.close_move", 0.02, 0.12, False, 0.07),
    ("plan_nasal", "plan:nasal.open_lead", 0.0, 0.12, False, 0.05),
    ("plan_nasal", "plan:nasal.anticipation", 0.0, 0.95, False, 0.5),
    # The spectrum's global shape: every vowel of ours sat under bdl at
    # 2.6-3.9 kHz (-4, -8, -10, -3 dB) and over at 4.5 kHz (+4). A speech
    # microphone's presence rise fills the first; a longer epilarynx drops
    # the second (2026-10-18).
    ("shape", "recording.presence_db", 0.0, 12.0, False, 0.0),
    ("shape", "recording.presence_hz", 2000.0, 4500.0, False, 3500.0),
    ("shape", "recording.presence_q", 0.5, 3.0, False, 0.7),
    ("shape", "anatomy.epilarynx_height", 1.8, 3.5, False, 2.3),
    ("shape", "anatomy.epilarynx_area", 0.5, 1.6, False, 1.0),
    # The recording chain's low cut (bdl sits 14-19 dB under his 150-300 Hz
    # band from 10 to 90 Hz, ours 0-5).
    ("noise", "recording.highpass", 0.0, 150.0, False, 0.0),
    ("plan_breath", "plan:breath.unstressed_pressure", 0.5, 1.2, False, 1.0),
    ("plan_breath", "plan:breath.syllable_rise", 0.1, 0.7, False, 0.35),
    ("plan_breath", "plan:breath.syllable_decay", 0.4, 1.0, False, 0.88),
    ("plan_breath", "plan:breath.voiceless_pressure", 0.5, 1.2, False, 0.85),
    # (ARCTIC lead-ins are silent; ours carried the inhalation.)
    ("plan_breath", "plan:breath.inhale_pressure", -150.0, 0.0, False, -110.0),
]
STATE = {}
HELDOUT_SEEDS = (render_set.SEED, "11")
# Renders at once for an evaluation made while the other workers wait.
WIDE = int(os.environ.get("WORKERS", "10"))


def setup(state):
    STATE.update(state)


def values_of(unit, leaves):
    out = {}
    for (group, leaf, low, high, log, default), u in zip(leaves, unit):
        u = min(max(u, 0.0), 1.0)
        out[leaf] = float(numpy.exp(numpy.log(low) + u * (numpy.log(high) - numpy.log(low)))) if log else low + u * (high - low)
    return out


def unit_of(values, leaves):
    out = []
    for group, leaf, low, high, log, default in leaves:
        v = values[leaf]
        out.append((numpy.log(v) - numpy.log(low)) / (numpy.log(high) - numpy.log(low)) if log else (v - low) / (high - low))
    return numpy.clip(numpy.array(out), 0.0, 1.0)


def describe_rendered(values, ids, threads=1):
    folder = tempfile.mkdtemp(dir=STATE["scratch"])
    params = os.path.join(folder, "p.params")
    tree = {k: v for k, v in values.items() if not k.startswith("plan:")}
    plan_values = {k[5:]: v for k, v in values.items() if k.startswith("plan:")}
    open(params, "w").write(slider_server.with_values(STATE["base_text"], tree))
    plan = STATE["plan"]
    if plan_values:
        plan = os.path.join(folder, "p.plan")
        open(plan, "w").write(slider_server.with_values(STATE["plan_text"], plan_values))
    table = render_set.texts()
    jobs = [(name, table[name], folder, params, STATE["organs"], plan, STATE["extra"]) for name in ids]
    if any(result[1] != "ok" for result in render_set.render_many(jobs, threads)):
        shutil.rmtree(folder, ignore_errors=True)
        return None
    description = fidelity.describe(fidelity.ours_set(folder, ids, STATE.get("ref_dir")))
    shutil.rmtree(folder, ignore_errors=True)
    return description


NATO_TEXT = ("alfa. bravo. charlie. delta. echo. foxtrot. golf. hotel. india. juliett. kilo. lima. mike. november. "
             "oscar. papa. quebec. romeo. sierra. tango. uniform. victor. whiskey. x-ray. yankee. zulu.")


def nato_scores(values):
    """NATO mode: the 26 words rendered as one utterance; the reader's clips
    cut by our phones carried over by DTW (so the reference's features follow
    each render's mapping); even words train, odd words are held out."""
    folder = tempfile.mkdtemp(dir=STATE["scratch"])
    params = os.path.join(folder, "p.params")
    tree = {k: v for k, v in values.items() if not k.startswith("plan:")}
    plan_values = {k[5:]: v for k, v in values.items() if k.startswith("plan:")}
    open(params, "w").write(slider_server.with_values(STATE["base_text"], tree))
    plan = STATE["plan"]
    if plan_values:
        plan = os.path.join(folder, "p.plan")
        open(plan, "w").write(slider_server.with_values(STATE["plan_text"], plan_values))
    result = render_set.render(("nato", NATO_TEXT, folder, params, STATE["organs"], plan, STATE["extra"]))
    if result[1] != "ok":
        shutil.rmtree(folder, ignore_errors=True)
        return None, None
    reference, ours = fidelity.nato_sets(os.path.join(folder, "nato.wav"), os.path.join(folder, "nato.score"),
                                         os.path.expanduser("~/ducklips/voice/ref/nato"))
    shutil.rmtree(folder, ignore_errors=True)
    train = fidelity.score(fidelity.describe(reference[0::2]), fidelity.describe(ours[0::2]))
    held = fidelity.score(fidelity.describe(reference[1::2]), fidelity.describe(ours[1::2]))
    return train, held


def evaluate(unit, threads=1):
    values = values_of(unit, STATE["leaves"])
    if STATE.get("nato"):
        train, held = nato_scores(values)
        if train is None:
            return -100.0, values, {}
        train["heldout_combined"] = held["combined"]
        return float(train["combined"]), values, train
    ours = describe_rendered(values, STATE["train"], threads)
    if ours is None:
        return -100.0, values, {}
    result = fidelity.score(STATE["reference_train"], ours)
    return float(result["combined"]), values, result


def heldout(values):
    """The held-out score, each group's dB averaged over two render seeds (the
    stages' acceptance compares these; one seed's phase group alone ranged
    11.3-16.0 dB over a fixed voice)."""
    if STATE.get("nato"):
        return nato_scores(values)[1] or {}
    scores = []
    for seed in HELDOUT_SEEDS:
        render_set.SEED = seed
        ours = describe_rendered(values, STATE["heldout"], WIDE)
        if ours is None:
            return {}
        scores.append(fidelity.score(STATE["reference_heldout"], ours))
    render_set.SEED = HELDOUT_SEEDS[0]
    return {k: float(numpy.mean([one[k] for one in scores])) for k in scores[0]}


def main():
    name, ref_dir, train_ids, heldout_ids, params, organs, plan, generations = sys.argv[1:9]
    groups = set(sys.argv[9].split(",")) if len(sys.argv) > 9 else None
    leaves = [leaf for leaf in LEAVES if groups is None or leaf[0] in groups]
    workdir = os.path.join("/tmp/es/profile", name)
    os.makedirs(workdir, exist_ok=True)
    cache = os.path.join(workdir, "reference.pkl")
    if ref_dir == "nato":
        reference_train, reference_heldout = {}, {}
    elif os.path.exists(cache):
        reference_train, reference_heldout = pickle.load(open(cache, "rb"))
    else:
        reference_train = fidelity.describe(fidelity.arctic_set(ref_dir, fidelity.id_range(train_ids)))
        reference_heldout = fidelity.describe(fidelity.arctic_set(ref_dir, fidelity.id_range(heldout_ids)))
        pickle.dump((reference_train, reference_heldout), open(cache, "wb"))
    base_text = open(params).read()
    plan_text = open(plan).read()
    start = {}
    for group, leaf, low, high, log, default in leaves:
        found = slider_server.read_value(plan_text, leaf[5:]) if leaf.startswith("plan:") else slider_server.read_value(base_text, leaf)
        start[leaf] = found if found is not None else default
    state = {"scratch": workdir, "base_text": base_text, "plan_text": plan_text, "organs": os.path.abspath(organs), "plan": os.path.abspath(plan),
             "extra": {"KINEMATIC": "1"}, "leaves": leaves,
             "train": fidelity.id_range(train_ids) if ref_dir != "nato" else [],
             "heldout": fidelity.id_range(heldout_ids) if ref_dir != "nato" else [], "reference_train": reference_train,
             "reference_heldout": reference_heldout, "nato": ref_dir == "nato",
             "ref_dir": ref_dir if ref_dir != "nato" else None}
    setup(state)
    mean = unit_of(start, leaves)
    # The first step size (unit space). Near an optimum 0.15 put every
    # candidate about a unit from the mean over 40 settings: bdl4 round 2's
    # tree stage found nothing better in 10 generations (all 6.8-8.0 against
    # 8.27 at the start).
    sigma = numpy.full(len(leaves), float(os.environ.get("PROFILE_SIGMA", "0.15")))
    rng = numpy.random.default_rng(abs(hash(name)) % (2 ** 32))
    population, parents = 9, 4
    weights = numpy.log(parents + 0.5) - numpy.log(numpy.arange(1, parents + 1))
    weights /= weights.sum()
    log = open(os.path.join(workdir, "log.jsonl"), "a")
    with multiprocessing.Pool(9, initializer=setup, initargs=(state,)) as pool:
        best = pool.apply(evaluate, (mean, WIDE))
        best_heldout = pool.apply(heldout, (best[1],))
        print("start train %s | heldout %s" % (fidelity.report(best[2]), fidelity.report(best_heldout)), flush=True)
        for generation in range(int(generations)):
            candidates = numpy.clip(mean + sigma * rng.standard_normal((population, len(leaves))), 0.0, 1.0)
            results = pool.map(evaluate, list(candidates))
            order = numpy.argsort([-r[0] for r in results])
            chosen = candidates[order[:parents]]
            spread = numpy.sqrt(numpy.sum(weights[:, None] * (chosen - mean) ** 2, axis=0))
            mean = numpy.sum(weights[:, None] * chosen, axis=0)
            sigma = numpy.clip(0.7 * sigma + 0.3 * spread * 1.4, 0.02, 0.3)
            top = results[order[0]]
            improved = top[0] > best[0]
            if improved:
                best = top
                best_heldout = pool.apply(heldout, (best[1],))
                header = "# profile_fit.py %s: train %s | heldout %s\n" % (name, fidelity.report(best[2]), fidelity.report(best_heldout))
                open(os.path.join(workdir, "best.params"), "w").write(header + slider_server.with_values(
                    base_text, {k: v for k, v in best[1].items() if not k.startswith("plan:")}))
                open(os.path.join(workdir, "best.plan"), "w").write(header + slider_server.with_values(
                    plan_text, {k[5:]: v for k, v in best[1].items() if k.startswith("plan:")}))
            log.write(json.dumps({"generation": generation, "train": top[0], "best": best[0], "values": top[1],
                                  "groups": top[2], "heldout": best_heldout}) + "\n")
            log.flush()
            print("generation %2d best train %.2f (this %.2f)%s" % (generation, best[0], top[0],
                  " heldout %s" % fidelity.report(best_heldout) if improved else ""), flush=True)
    print("best values: " + " ".join("%s=%.4g" % (k, v) for k, v in best[1].items()))


if __name__ == "__main__":
    main()
