"""Fit the VTL-2025 glottis and the body's variation to a human reading of
the same words: the NATO reader (voice/ref/nato/nato_[a-z].wav, the owner's
choice 2026-10-18: "id rather you fit to the nato voice"), or the owner's
sentences (voice/ref/user_three.wav 0-9.8 s: "hello my name is ducklips, the
fox jumped over the lazy dog. she sells seashells by the seashore"; the rest
is the same whispered). Ours says the same words (say, rate 0.75,
organs_skeleton); both are measured over every voiced run of at least
0.12 s (cord_features), medians over the runs pooled over the files.

Fitted: tilt, H1-H2, harmonics over troughs per band, the harmonics' phase
scatter per band (np.angle, phase distortion deviation) and pitch flutter.
Not fitted: the phase distortion's mean (the owner's +0.9 rad against bdl's
+2.6 to +3.0: it differs by speaker or recording), drift, tremor and
wobble (sentences carry intonation). A click rate above the owner's is
penalized (the held-vowel fit shut the glottis: 26 clicks/s).
usage: cord_fit_sentences.py BASE.params SCORE OUT.params GENERATIONS REFERENCE.wav[:START:END] ..."""
import json
import multiprocessing
import os
import subprocess
import sys
import tempfile

import librosa
import numpy
import soundfile
from scipy.signal import butter, sosfiltfilt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cord_features
import cord_fit
import slider_server

NAMES = ["tilt", "h1h2", "hnr_0_1k", "hnr_1_2k", "hnr_2_4k", "hnr_4_6k", "pdd_0_1k", "pdd_1_2k", "pdd_2_4k", "flutter"]
SCALES = {"tilt": 2.0, "h1h2": 2.5, "hnr_0_1k": 3.0, "hnr_1_2k": 2.5, "hnr_2_4k": 2.5, "hnr_4_6k": 2.5,
          "pdd_0_1k": 0.15, "pdd_1_2k": 0.15, "pdd_2_4k": 0.15, "flutter": 2.5}
SCORE = ""
TARGETS = {}
OWNER_CLICKS = 0.0


def setup(base_text, score, targets, owner_clicks):
    global SCORE, TARGETS, OWNER_CLICKS
    cord_fit.set_base(base_text)
    SCORE = score
    TARGETS = targets
    OWNER_CLICKS = owner_clicks


def clicks_per_second(x, r, k=8.0):
    # As analysis/clicks2.py: the 5 kHz high-pass over 8 times its 20 ms rms, while loud.
    h = sosfiltfilt(butter(4, 5000, btype="high", fs=r, output="sos"), x)
    w = int(0.02 * r)
    e = numpy.sqrt(numpy.convolve(h ** 2, numpy.ones(w) / w, "same")) + 1e-9
    level = numpy.sqrt(numpy.convolve(x ** 2, numpy.ones(w) / w, "same"))
    loud = level > level.max() * 10 ** (-40 / 20)
    times = []
    for i in numpy.where((numpy.abs(h) > k * e) & loud)[0]:
        t = i / r
        if not times or t - times[-1] > 0.003:
            times.append(t)
    return len(times) / (len(x) / r)


def voiced_runs(x, r):
    hop = 256
    f0, voiced, _ = librosa.pyin(x, fmin=60, fmax=320, sr=r, frame_length=2048, hop_length=hop)
    rows = []
    start = None
    for i, on in enumerate(list(voiced) + [False]):
        if on and start is None:
            start = i
        elif not on and start is not None:
            if (i - start) * hop / r >= 0.12:
                found = cord_features.features(x[start * hop:i * hop], r)
                if found:
                    rows.append(found)
            start = None
    return rows


def pooled(rows):
    if len(rows) < 5:
        return None
    out = {}
    for name in NAMES:
        values = [row[name] for row in rows if numpy.isfinite(row[name])]
        out[name] = float(numpy.median(values)) if values else float("nan")
    out["runs"] = len(rows)
    return out


def evaluate(unit):
    values = cord_fit.values_of(unit)
    folder = tempfile.mkdtemp(dir="/tmp/es/cordfit_s")
    params = os.path.join(folder, "p.params")
    open(params, "w").write(slider_server.with_values(cord_fit.BASE_TEXT, values))
    wave = os.path.join(folder, "s.wav")
    done = subprocess.run([slider_server.SPEAK, SCORE, wave, "73021"], env=dict(os.environ, KINEMATIC="1", PARAMS=params),
                          capture_output=True, text=True)
    if done.returncode != 0:
        return 1.0e3, values, None, {}
    x, r = soundfile.read(wave)
    subprocess.run(["rm", "-rf", folder])
    found = pooled(voiced_runs(x, r))
    if found is None:
        return 1.0e3, values, None, {}
    found["clicks"] = clicks_per_second(x, r)
    terms = {}
    for name in NAMES:
        value = found[name]
        terms[name] = ((value - TARGETS[name]) / SCALES[name]) ** 2 if numpy.isfinite(value) else 25.0
    terms["clicks"] = max(found["clicks"] - OWNER_CLICKS, 0.0) ** 2
    return float(sum(terms.values()) / len(terms)), values, found, terms


def main():
    base_text = open(sys.argv[1]).read()
    score = sys.argv[2]
    out = sys.argv[3]
    generations = int(sys.argv[4])
    os.makedirs("/tmp/es/cordfit_s", exist_ok=True)
    rows = []
    click_count = 0.0
    seconds = 0.0
    for spec in sys.argv[5:]:
        parts = spec.split(":")
        x, r = soundfile.read(parts[0])
        x = x if x.ndim == 1 else x.mean(axis=1)
        if len(parts) == 3:
            x = x[int(float(parts[1]) * r):int(float(parts[2]) * r)]
        rows += voiced_runs(x, r)
        click_count += clicks_per_second(x, r) * len(x) / r
        seconds += len(x) / r
    targets = pooled(rows)
    owner_clicks = click_count / seconds
    print("targets (%d runs): %s | clicks %.1f/s" % (targets["runs"], " ".join("%s %.2f" % (n, targets[n]) for n in NAMES),
                                                     owner_clicks), flush=True)
    setup(base_text, score, targets, owner_clicks)
    log = open("/tmp/es/cordfit_s/log.jsonl", "a")
    start = {}
    for leaf, low, high, log_scale in cord_fit.SETTINGS:
        found = slider_server.read_value(base_text, leaf)
        start[leaf] = found if found is not None else cord_fit.ENGINE_DEFAULTS.get(leaf, low)
    mean = cord_fit.unit_of(start)
    sigma = numpy.full(len(cord_fit.SETTINGS), 0.2)
    rng = numpy.random.default_rng(20261019)
    population = 9
    parents = 4
    weights = numpy.log(parents + 0.5) - numpy.log(numpy.arange(1, parents + 1))
    weights /= weights.sum()

    def show(result):
        found = result[2]
        return "cost %.3f | %s | clicks %.1f/s" % (result[0], " ".join("%s %.2f" % (n, found[n]) for n in NAMES), found["clicks"]) if found else "failed"

    with multiprocessing.Pool(9, initializer=setup, initargs=(base_text, score, targets, owner_clicks)) as pool:
        best = pool.apply(evaluate, (mean,))
        print("start " + show(best), flush=True)
        for generation in range(generations):
            candidates = numpy.clip(mean + sigma * rng.standard_normal((population, len(cord_fit.SETTINGS))), 0.0, 1.0)
            results = pool.map(evaluate, list(candidates))
            order = numpy.argsort([result[0] for result in results])
            for index in order:
                log.write(json.dumps({"generation": generation, "cost": results[index][0], "values": results[index][1],
                                      "terms": results[index][3]}) + "\n")
            log.flush()
            chosen = candidates[order[:parents]]
            spread = numpy.sqrt(numpy.sum(weights[:, None] * (chosen - mean) ** 2, axis=0))
            mean = numpy.sum(weights[:, None] * chosen, axis=0)
            sigma = numpy.clip(0.7 * sigma + 0.3 * spread * 1.4, 0.02, 0.35)
            if results[order[0]][0] < best[0]:
                best = results[order[0]]
            print("generation %2d best %.3f this %.3f" % (generation, best[0], results[order[0]][0]), flush=True)
    total, values, found, terms = best
    open(out, "w").write("# cord_fit_sentences.py: VTL-2025 glottis and variation fitted to the owner's sentences, cost %.3f\n" % total +
                         slider_server.with_values(base_text, values))
    print("best " + show(best))
    for leaf in values:
        print("  %-28s %.5g (start %.5g)" % (leaf, values[leaf], start[leaf]))
    print("terms: " + " ".join("%s %.2f" % (k, v) for k, v in sorted(terms.items(), key=lambda kv: -kv[1])))


if __name__ == "__main__":
    main()
