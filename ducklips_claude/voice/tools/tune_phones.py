"""Analysis by synthesis: tune each phone's articulation so that the
automatic phone test identifies it as intended, judged against human
exemplars (CMU ARCTIC bdl). No listener in the loop.

For a phone and one parameter, each candidate value is written to a phone
override file, the utterances containing that phone are rendered and
analysed in parallel, and the candidate scores

  score = accuracy of that phone + 0.5 x its mean vote share
          - penalty if overall accuracy on those utterances falls

The best value is kept and the next phone is tuned on top of it.

usage: tune_phones.py PLAN.json OVERRIDES.txt
  PLAN.json: [{"phone": "T", "key": "place", "values": [14.5, 15.0, 15.5]}, ...]
"""
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
import phone_test as pt

ROOT = ca.ROOT
TUNE = os.path.join(ROOT, "out", "tune")
VIEWS = ["cq", "causal"]


def training_set():
    cache = os.path.join(TUNE, "training_%s.npz" % "_".join(VIEWS))
    if os.path.exists(cache):
        data = numpy.load(cache)
        return data["labels"], data["normal"], data["mean"], data["spread"]
    names = ["arctic_a%04d" % n for n in range(1, 41)]
    train = pt.collect(names, "ref", VIEWS)
    labels = numpy.array([s[0] for s in train])
    matrix = numpy.array([s[1] for s in train])
    mean = matrix.mean(axis=0)
    spread = matrix.std(axis=0) + 1.0
    normal = (matrix - mean) / spread
    normal /= numpy.linalg.norm(normal, axis=1, keepdims=True) + 1e-9
    os.makedirs(TUNE, exist_ok=True)
    numpy.savez(cache, labels=labels, normal=normal, mean=mean, spread=spread)
    return labels, normal, mean, spread


def vote(model, vector, k=5):
    labels, normal, mean, spread = model
    v = (vector - mean) / spread
    v /= numpy.linalg.norm(v) + 1e-9
    similarity = normal @ v
    nearest = numpy.argsort(-similarity)[:k]
    votes = {}
    for index in nearest:
        votes[labels[index]] = votes.get(labels[index], 0.0) + max(similarity[index], 0.0)
    total = sum(votes.values()) + 1e-9
    guess = max(votes.items(), key=lambda item: item[1])[0]
    return guess, votes


def utterances_with(phone, limit=8):
    counts = []
    for number in range(1, 61):
        name = "arctic_a%04d" % number
        phones = [p for p, s, e in ca.our_phones(os.path.join(ca.OUT, name + ".score"))]
        hits = sum(1 for p in phones if p == phone.lower())
        if hits:
            counts.append((hits, name))
    counts.sort(reverse=True)
    return [name for hits, name in counts[:limit]]


def render_and_analyse(args):
    name, text, overrides, folder = args
    os.makedirs(folder, exist_ok=True)
    wave48 = os.path.join(folder, name + "_48.wav")
    score = os.path.join(folder, name + ".score")
    subprocess.run([os.path.join(ROOT, "build", "say"), text, wave48, "--score", score, "--phones", overrides],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    data, rate = ca.read_wave(wave48)
    wave16 = os.path.join(folder, name + "_ours.wav")
    ca.write_wave(wave16, ca.resample(data, rate, 16000), 16000)
    subprocess.run([os.path.join(ROOT, "build", "tf_analyze"), wave16, os.path.join(folder, name + "_ours"),
                    "--views", ",".join(VIEWS)],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return name


def score_candidate(model, phone, names, folder):
    target = phone.lower()
    hits, shares, overall = [], [], []
    for name in names:
        prefix = os.path.join(folder, name + "_ours")
        fields = [pt.load_levels(prefix, view) for view in VIEWS]
        for p, start, end in ca.our_phones(os.path.join(folder, name + ".score")):
            p = {"dx": "t", "ax": "ah"}.get(p, p)
            guess, votes = vote(model, pt.fingerprint(fields, start, end))
            overall.append(guess == p)
            if p == target:
                hits.append(guess == p)
                shares.append(votes.get(p, 0.0) / (sum(votes.values()) + 1e-9))
    return float(numpy.mean(hits)), float(numpy.mean(shares)), float(numpy.mean(overall))


def write_overrides(path, settings):
    with open(path, "w") as handle:
        handle.write("# phone overrides found by tools/tune_phones.py against CMU ARCTIC bdl\n")
        for phone, values in sorted(settings.items()):
            handle.write(phone + "".join(" %s %g" % (k, v) for k, v in sorted(values.items())) + "\n")


def read_overrides(path):
    settings = {}
    if os.path.exists(path):
        for line in open(path):
            parts = line.split()
            if not parts or parts[0].startswith("#"):
                continue
            settings[parts[0]] = {parts[i]: float(parts[i + 1]) for i in range(1, len(parts) - 1, 2)}
    return settings


def main():
    plan = json.load(open(sys.argv[1]))
    output = sys.argv[2]
    texts = ca.transcripts()
    model = training_set()
    settings = read_overrides(output)
    for step in plan:
        phone, key = step["phone"], step["key"]
        names = utterances_with(phone)
        if not names:
            print("%s: no utterances" % phone)
            continue
        # Every candidate of this step renders and analyses at once, so all
        # cores stay busy; scoring follows.
        folders = []
        jobs = []
        for value in step["values"]:
            trial = {p: dict(v) for p, v in settings.items()}
            trial.setdefault(phone, {})[key] = value
            folder = os.path.join(TUNE, "%s_%s_%g" % (phone, key, value))
            overrides = os.path.join(folder, "phones.txt")
            os.makedirs(folder, exist_ok=True)
            write_overrides(overrides, trial)
            folders.append((value, folder))
            jobs.extend((name, texts[name], overrides, folder) for name in names)
        with ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as pool:
            list(pool.map(render_and_analyse, jobs))
        results = []
        for value, folder in folders:
            accuracy, share, overall = score_candidate(model, phone, names, folder)
            results.append((accuracy + 0.5 * share + overall, value, accuracy, share, overall))
            print("  %-3s %s=%-6g acc %.2f share %.2f overall %.3f" % (phone, key, value, accuracy, share, overall), flush=True)
            shutil.rmtree(folder, ignore_errors=True)
        best = max(results)
        settings.setdefault(phone, {})[key] = best[1]
        write_overrides(output, settings)
        print("%s %s -> %g (acc %.2f, overall %.3f)" % (phone, key, best[1], best[2], best[4]), flush=True)


if __name__ == "__main__":
    main()
