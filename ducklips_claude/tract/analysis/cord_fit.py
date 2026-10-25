"""Fit the VocalTractLab-2025 glottis (and the body's variation) to the human
held-vowel features (cord_features.py), owner 2026-10-18.

Targets (validated 2026-10-18; bdl = ARCTIC aa / iy / uw over 0.15 s, 85
stretches; owner = voice/ref/user_three.wav's held stretches, vowel unknown):
  tilt and h1h2 per vowel from bdl, and their means over the vowels (the
    source's share; the tract is ours, not bdl's, so per vowel weighs less);
  harmonics over troughs per band, owner and bdl pooled (they agree);
  pitch flutter 8-20 Hz from both; tremor 3-8 Hz, drift 0.5-3 Hz and level
    wobble 3-20 Hz from the owner (bdl's vowels are too short to resolve them).
Search: a diagonal evolution strategy over the settings scaled to [0, 1],
nine candidates a generation rendered in parallel (held ah / ee / oo each,
measured over 1.0-2.0 s).
usage: cord_fit.py BASE.params OUT.params [GENERATIONS]"""
import json
import multiprocessing
import os
import subprocess
import sys
import tempfile

import numpy
import soundfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cord_features
import slider_server

# (leaf, low, high, log scale)
SETTINGS = [
    ("posture.voicing_lower", 0.0, 0.08, False),
    ("posture.voicing_upper", 0.0, 0.10, False),
    ("posture.voicing_chink", 0.0, 0.05, False),
    ("kinematic.phase_lag", 0.3, 2.0, False),
    ("kinematic.pulse_shape", -1.0, 1.0, False),
    ("kinematic.amplitude", 0.0003, 0.0015, True),
    ("kinematic.flutter", 0.0, 60.0, False),
    ("kinematic.jitter", 0.0, 1.5, False),    # (0.6 was a wall in the first fit)
    ("vtl.aspiration_gain", 0.2, 20.0, True),
    ("variation.breath", 0.0, 0.15, False),
    ("variation.tension", 0.0, 0.05, False),
    ("variation.tremor", 0.0, 0.03, False),
    ("variation.breath_tremor", 0.0, 0.3, False),  # (0.1 was a wall in the first fit)
]
PER_VOWEL = {"tilt": ({"ah": -11.0, "ee": -0.9, "oo": -9.7}, 6.0),
             "h1h2": ({"ah": -2.5, "ee": -6.4, "oo": -5.0}, 6.0)}
MEAN_OVER_VOWELS = {"tilt": (-7.2, 2.0), "h1h2": (-4.6, 2.0)}
POOLED = {"hnr_0_1k": (20.5, 3.0), "hnr_1_2k": (8.5, 2.5), "hnr_2_4k": (7.0, 2.5), "hnr_4_6k": (5.8, 2.5),
          "flutter": (6.5, 2.0), "tremor": (14.0, 7.0), "drift": (9.0, 10.0), "wobble": (1.5, 0.7)}
VOWELS = ("ah", "ee", "oo")
ENGINE_DEFAULTS = {"vtl.aspiration_gain": 1.0}  # settings a base file may leave to the engine
HELD = (1.0, 2.0)
BASE_TEXT = ""


def set_base(text):
    # Each worker's copy (macOS starts pool workers fresh, without main's globals).
    global BASE_TEXT
    BASE_TEXT = text


def values_of(x):
    values = {}
    for (leaf, low, high, log), unit in zip(SETTINGS, x):
        u = min(max(unit, 0.0), 1.0)
        values[leaf] = float(numpy.exp(numpy.log(low) + u * (numpy.log(high) - numpy.log(low)))) if log else low + u * (high - low)
    return values


def unit_of(values):
    x = []
    for leaf, low, high, log in SETTINGS:
        v = values[leaf]
        x.append((numpy.log(v) - numpy.log(low)) / (numpy.log(high) - numpy.log(low)) if log else (v - low) / (high - low))
    return numpy.clip(numpy.array(x), 0.0, 1.0)


def measure(values):
    folder = tempfile.mkdtemp(dir="/tmp/es/cordfit")
    params = os.path.join(folder, "p.params")
    open(params, "w").write(slider_server.with_values(BASE_TEXT, values))
    found = {}
    for vowel in VOWELS:
        wave = os.path.join(folder, vowel + ".wav")
        done = subprocess.run([slider_server.SPEAK, os.path.join(slider_server.SCORES, vowel + ".score"), wave, "73021"],
                              env=dict(os.environ, KINEMATIC="1", PARAMS=params), capture_output=True, text=True)
        if done.returncode != 0:
            return None
        x, r = soundfile.read(wave)
        row = cord_features.features(x[int(HELD[0] * r):int(HELD[1] * r)], r)
        if row is None:
            return None
        found[vowel] = row
    subprocess.run(["rm", "-rf", folder])
    return found


def cost(found):
    if found is None:
        return 1.0e3, {}
    terms = {}
    for name, (targets, scale) in PER_VOWEL.items():
        for vowel in VOWELS:
            terms["%s_%s" % (name, vowel)] = ((found[vowel][name] - targets[vowel]) / scale) ** 2
    for name, (target, scale) in MEAN_OVER_VOWELS.items():
        terms["mean_" + name] = ((numpy.mean([found[v][name] for v in VOWELS]) - target) / scale) ** 2
    for name, (target, scale) in POOLED.items():
        value = numpy.nanmean([found[v][name] for v in VOWELS])
        terms[name] = ((value - target) / scale) ** 2 if numpy.isfinite(value) else 25.0
    return float(sum(terms.values()) / len(terms)), terms


def evaluate(x):
    values = values_of(x)
    found = measure(values)
    total, terms = cost(found)
    return total, values, found, terms


def summary(found):
    rows = []
    for name in cord_features.FEATURES:
        rows.append("%s %s" % (name, " ".join("%6.1f" % found[v][name] for v in VOWELS)))
    return " | ".join(rows)


def main():
    global BASE_TEXT
    BASE_TEXT = open(sys.argv[1]).read()
    out = sys.argv[2]
    generations = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    os.makedirs("/tmp/es/cordfit", exist_ok=True)
    log = open("/tmp/es/cordfit/log.jsonl", "a")
    start = {}
    for leaf, low, high, log_scale in SETTINGS:
        found = slider_server.read_value(BASE_TEXT, leaf)
        start[leaf] = found if found is not None else ENGINE_DEFAULTS.get(leaf, low)
    mean = unit_of(start)
    sigma = numpy.full(len(SETTINGS), 0.2)
    rng = numpy.random.default_rng(20261018)
    population = 9
    parents = 4
    weights = numpy.log(parents + 0.5) - numpy.log(numpy.arange(1, parents + 1))
    weights /= weights.sum()
    with multiprocessing.Pool(9, initializer=set_base, initargs=(BASE_TEXT,)) as pool:
        best = pool.apply(evaluate, (mean,))
        print("start cost %.3f | %s" % (best[0], summary(best[2]) if best[2] else "failed"), flush=True)
        for generation in range(generations):
            steps = rng.standard_normal((population, len(SETTINGS)))
            candidates = numpy.clip(mean + sigma * steps, 0.0, 1.0)
            results = pool.map(evaluate, list(candidates))
            order = numpy.argsort([result[0] for result in results])
            for index in order:
                log.write(json.dumps({"generation": generation, "cost": results[index][0], "values": results[index][1],
                                      "terms": results[index][3]}) + "\n")
            log.flush()
            chosen = candidates[order[:parents]]
            spread = numpy.sqrt(numpy.sum(weights[:, None] * (chosen - mean) ** 2, axis=0))
            mean = numpy.sum(weights[:, None] * chosen, axis=0)
            # Each setting's step follows the spread of the chosen, damped.
            sigma = numpy.clip(0.7 * sigma + 0.3 * spread * 1.4, 0.02, 0.35)
            if results[order[0]][0] < best[0]:
                best = results[order[0]]
            print("generation %2d best %.3f this %.3f sigma %.3f" % (generation, best[0], results[order[0]][0], sigma.mean()), flush=True)
    total, values, found, terms = best
    open(out, "w").write("# cord_fit.py: VTL-2025 glottis and body variation fitted to human held-vowel features, cost %.3f\n" % total +
                         slider_server.with_values(BASE_TEXT, values))
    print("best cost %.3f" % total)
    for leaf in values:
        print("  %-28s %.5g (start %.5g)" % (leaf, values[leaf], start[leaf]))
    print("features (ah ee oo): " + summary(found))
    print("terms: " + " ".join("%s %.2f" % (k, v) for k, v in sorted(terms.items(), key=lambda kv: -kv[1])))


if __name__ == "__main__":
    main()
