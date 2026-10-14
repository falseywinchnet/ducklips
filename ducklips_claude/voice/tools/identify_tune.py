"""Tune the phone inventory for four-view phone identification: each
candidate override renders ARCTIC a0001-a0060 (every view analysed) and is
judged by tools/phone_test.py's nearest-neighbour vote against the human
exemplars (bdl a0001-a0040, stft/cq/causal/wigner). The objective is our
accuracy over a0001-a0040; a0041-a0060 is reported beside it as the check.
Each candidate is scored over SEEDS (the articulators' wander and the
turbulence are random; tuned at one seed, the gains partly fitted that
realisation: 40.2% at seed 1 read 37.7-38.4% once the sequence shifted).
The best value of each step is kept before the next step, unless it costs
the check sentences more than CHECK_SLACK (one or two phones there).

With PROFILE=s,sh in the environment the objective is instead spectral
shape: the negative RMS (dB) between our and bdl's mean band profiles of
those phones at their middle slice, 1-7 kHz, constant offset removed
(constant-Q view) - what a listener hears as a lisp or a dull hiss. An
entry PHONE@SLICE (slice 0-6: -30 ms, start, 25%, 50%, 75%, end, +30 ms)
compares that slice instead and keeps its level re the utterance, over
1-7 kHz: k@4 is a /k/ about 15-30 ms after its burst in both bdl's labels
and ours (the burst-frication a listener needs to hear /k/, not /p/).

With BOTH=1 as well, the reference profile is the mean of bdl's and rms's
(ref/arctic_rms, analysed into out/arctic_rms): a target that is no single
speaker's habit.

usage: identify_tune.py PLAN.json OVERRIDES.txt
  PLAN.json: [{"phone": "L", "key": "vplace", "values": [6, 7, 8]}, ...]
  Global settings use the phone "*" with keys aspiration, frication, pressure.
  A value may be an object of several keys set together (key is then a label):
  {"phone": "L", "key": "fit", "values": [{"place": 14.6, "vplace": 4}]}
"""
import json
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
import phone_test as pt
import render_arctic as ra
import tune_phones as tp

VIEWS = ["stft", "cq", "causal", "wigner"]
WORK = os.path.join(ca.ROOT, "out", "identify_tune")
CHECK_SLACK = 0.3  # percent, or dB of shape error under PROFILE
HUMAN_TUNED = {}
HUMAN_CHECK = {}
SEEDS = (1, 2)
TUNED = ["arctic_a%04d" % n for n in range(1, 41)]
CHECK = ["arctic_a%04d" % n for n in range(41, 61)]


PROFILE = [p for p in os.environ.get("PROFILE", "").split(",") if p]
BAND_HZ = numpy.array([ca.GRID[a:b].mean() for a, b in pt.BANDS])
SHAPE_BANDS = (BAND_HZ >= 1000) & (BAND_HZ <= 7000)
MIDDLE = 3  # slice index of the phone's 50% point


def middle_profiles(samples, entries):
    profiles = {}
    for entry in entries:
        phone, _, slice_text = entry.partition("@")
        slice_index = int(slice_text) if slice_text else MIDDLE
        rows = [s[1] for s in samples if s[0] == phone]
        if rows:
            profiles[entry] = numpy.mean(rows, axis=0).reshape(len(pt.OFFSETS), -1)[slice_index]
    return profiles


def shape_score(ours, human):
    errors = []
    for entry, reference in human.items():
        if entry in ours:
            difference = (ours[entry] - reference)[SHAPE_BANDS]
            offset = 0.0 if "@" in entry else difference.mean()
            errors.append(numpy.sqrt(numpy.mean((difference - offset) ** 2)))
    return -float(numpy.mean(errors))


def second_speaker(names):
    """rms's samples (constant-Q fingerprints) for the named utterances."""
    folder = os.path.join(ca.ROOT, "out", "arctic_rms")
    labels = os.path.join(ca.ROOT, "ref", "arctic_rms", "lab")
    samples = []
    for name in names:
        fields = [pt.load_levels(os.path.join(folder, name + "_rms"), "cq")]
        start = 0.0
        for line in open(os.path.join(labels, name + ".lab")):
            parts = line.split()
            if len(parts) != 3:
                continue
            end = float(parts[0])
            phone = {"dx": "t", "ax": "ah"}.get(parts[2], parts[2])
            if phone not in ("pau", "h#", "ssil"):
                samples.append((phone, pt.fingerprint(fields, start, end), name))
            start = end
    return samples


def accuracy(results, known):
    results = [r for r in results if r[0] in known]
    return 100.0 * sum(1 for a, b in results if a == b) / max(len(results), 1)


def apply_value(settings, globals_, phone, key, value):
    target = globals_ if phone == "*" else settings.setdefault(phone, {})
    if isinstance(value, dict):
        target.update(value)
    else:
        target[key] = value


def describe(value):
    if isinstance(value, dict):
        return ",".join("%s=%g" % (k, v) for k, v in sorted(value.items()))
    return "%g" % value


def score(settings, globals_, train, known):
    os.makedirs(WORK, exist_ok=True)
    overrides_path = os.path.join(WORK, "phones.txt")
    tp.write_overrides(overrides_path, settings)
    extra = ["--phones", overrides_path]
    for key in sorted(globals_):
        extra += ["--" + key, "%g" % globals_[key]]
    tuned = 0.0
    check = 0.0
    for seed in SEEDS:
        ra.render_all(1, 60, extra + ["--seed", str(seed)], True)
        if PROFILE:
            tuned += shape_score(middle_profiles(pt.collect(TUNED, "ours", ["cq"]), PROFILE), HUMAN_TUNED) / len(SEEDS)
            check += shape_score(middle_profiles(pt.collect(CHECK, "ours", ["cq"]), PROFILE), HUMAN_CHECK) / len(SEEDS)
            continue
        tuned += accuracy(pt.classify(train, pt.collect(TUNED, "ours", VIEWS), 5), known) / len(SEEDS)
        check += accuracy(pt.classify(train, pt.collect(CHECK, "ours", VIEWS), 5), known) / len(SEEDS)
    return tuned, check


def main():
    plan = json.load(open(sys.argv[1]))
    output = sys.argv[2]
    settings = tp.read_overrides(output)
    globals_ = settings.pop("*", {})
    global HUMAN_TUNED, HUMAN_CHECK
    if PROFILE:
        HUMAN_TUNED = middle_profiles(pt.collect(TUNED, "ref", ["cq"]), PROFILE)
        HUMAN_CHECK = middle_profiles(pt.collect(CHECK, "ref", ["cq"]), PROFILE)
        if os.environ.get("BOTH"):
            second_tuned = middle_profiles(second_speaker(TUNED), PROFILE)
            second_check = middle_profiles(second_speaker(CHECK), PROFILE)
            HUMAN_TUNED = {k: 0.5 * (v + second_tuned[k]) for k, v in HUMAN_TUNED.items() if k in second_tuned}
            HUMAN_CHECK = {k: 0.5 * (v + second_check[k]) for k, v in HUMAN_CHECK.items() if k in second_check}
    train = pt.collect(["arctic_a%04d" % n for n in range(1, 41)], "ref", VIEWS)
    known = {s[0] for s in train}
    best, check = score(settings, globals_, train, known)
    print("baseline identification %.2f (check %.2f)" % (best, check), flush=True)
    for step in plan:
        phone, key = step["phone"], step["key"]
        results = []
        for value in step["values"]:
            trial_settings = {p: dict(v) for p, v in settings.items()}
            trial_globals = dict(globals_)
            apply_value(trial_settings, trial_globals, phone, key, value)
            tuned, held = score(trial_settings, trial_globals, train, known)
            results.append((tuned, held, len(results)))
            print("  %-3s %s=%-7s identification %.2f (check %.2f)" % (phone, key, describe(value), tuned, held), flush=True)
        admissible = [r for r in results if r[1] >= check - CHECK_SLACK]
        tuned, held, index = max(admissible) if admissible else (best, check, None)
        if tuned > best:
            best, check = tuned, held
            value = step["values"][index]
            apply_value(settings, globals_, phone, key, value)
            print("%s %s -> %s (identification %.2f, check %.2f)" % (phone, key, describe(value), best, check), flush=True)
        else:
            print("%s %s unchanged (best %.2f)" % (phone, key, best), flush=True)
        saved = {p: dict(v) for p, v in settings.items()}
        if globals_:
            saved["*"] = globals_
        tp.write_overrides(output, saved)
    print("final identification %.2f (check %.2f)" % (best, check), flush=True)


if __name__ == "__main__":
    main()
