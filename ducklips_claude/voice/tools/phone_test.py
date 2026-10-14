"""Automatic phone-identification test: are our phones recognisable as the
intended phones by the standard of a human speaker?

Each phone becomes a fingerprint: 7 time slices (30 ms before, start, 25%,
50%, 75%, end, 30 ms after) x 38 bands (100 Hz - 7.6 kHz, ~1/6 octave), in
dB relative to the utterance's loud frames, from the constant-Q and causal
Gabor views. A k-nearest-neighbour vote over human exemplars (CMU ARCTIC bdl
a0001-a0040, from their phone labels) names each test phone. It is a
template match, not a recogniser: no model of words or language.

Reported:
  human held-out (a0041-a0060)  - ceiling of this test
  ours (a0001-a0060)            - our phones judged against the human ones
Both as phone accuracy, per-class accuracy, and our worst confusions.

usage: phone_test.py [--views cq,causal] [--k 5] [--json OUT]
"""
import json
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca

OFFSETS = (-0.03, 0.0, 0.25, 0.5, 0.75, 1.0, 0.03)  # seconds for the outer two, fractions otherwise
FLOOR = -60.0


def band_matrix():
    groups = []
    for start in range(0, len(ca.GRID) - 4 + 1, 4):
        groups.append((start, start + 4))
    return groups


BANDS = band_matrix()


def load_levels(prefix, view):
    mapped, seconds = ca.load_view(prefix, view)
    power = 10 ** (mapped / 10.0)
    reference = ca.utterance_level(mapped)
    banded = numpy.stack([power[a:b].mean(axis=0) for a, b in BANDS])
    level = 10 * numpy.log10(banded + 1e-14) - reference
    return numpy.maximum(level, FLOOR), seconds


def fingerprint(fields, start, end):
    parts = []
    for level, seconds in fields:
        for offset in OFFSETS:
            if offset == -0.03:
                time = start - 0.03
            elif offset == 0.03 and offset not in (0.0,):
                time = end + 0.03
            else:
                time = start + offset * (end - start)
            column = int(numpy.clip(numpy.searchsorted(seconds, time), 0, len(seconds) - 1))
            lo = max(column - 1, 0)
            parts.append(level[:, lo:column + 2].mean(axis=1))
    return numpy.concatenate(parts)


def collect(names, side, views):
    samples = []
    for name in names:
        prefix = os.path.join(ca.OUT, name + "_" + side)
        fields = [load_levels(prefix, view) for view in views]
        if side == "ref":
            phones = ca.reference_phones(name)
        else:
            phones = ca.our_phones(os.path.join(ca.OUT, name + ".score"))
        for phone, start, end in phones:
            phone = {"dx": "t", "ax": "ah"}.get(phone, phone)
            samples.append((phone, fingerprint(fields, start, end), name))
    return samples


def classify(train, test, k):
    labels = numpy.array([s[0] for s in train])
    matrix = numpy.array([s[1] for s in train])
    mean = matrix.mean(axis=0)
    spread = matrix.std(axis=0) + 1.0
    normal = (matrix - mean) / spread
    normal /= numpy.linalg.norm(normal, axis=1, keepdims=True) + 1e-9
    results = []
    for phone, vector, name in test:
        v = (vector - mean) / spread
        v /= numpy.linalg.norm(v) + 1e-9
        similarity = normal @ v
        nearest = numpy.argsort(-similarity)[:k]
        votes = {}
        for index in nearest:
            votes[labels[index]] = votes.get(labels[index], 0.0) + similarity[index]
        guess = max(votes.items(), key=lambda item: item[1])[0]
        results.append((phone, guess))
    return results


def summarise(label, results, known):
    results = [r for r in results if r[0] in known]
    correct = sum(1 for a, b in results if a == b)
    print("\n%s: %d phones, accuracy %.1f%%" % (label, len(results), 100.0 * correct / max(len(results), 1)))
    by_class = {}
    for a, b in results:
        by_class.setdefault(ca.phone_class(a), []).append(a == b)
    print("  " + "  ".join("%s %.0f%%" % (c, 100.0 * numpy.mean(v)) for c, v in sorted(by_class.items())))
    return 100.0 * correct / max(len(results), 1)


def main():
    views = ["cq", "causal"]
    k = 5
    for index, argument in enumerate(sys.argv):
        if argument == "--views":
            views = sys.argv[index + 1].split(",")
        if argument == "--k":
            k = int(sys.argv[index + 1])
    train_names = ["arctic_a%04d" % n for n in range(1, 41)]
    held_names = ["arctic_a%04d" % n for n in range(41, 61)]
    our_names = ["arctic_a%04d" % n for n in range(1, 61)]
    train = collect(train_names, "ref", views)
    known = {s[0] for s in train}
    human = classify(train, collect(held_names, "ref", views), k)
    ours = classify(train, collect(our_names, "ours", views), k)
    human_accuracy = summarise("human held-out (ceiling)", human, known)
    our_accuracy = summarise("ours", ours, known)
    confusions = {}
    for a, b in ours:
        if a != b and a in known:
            confusions[(a, b)] = confusions.get((a, b), 0) + 1
    print("\nour most frequent confusions (intended -> identified):")
    for (a, b), count in sorted(confusions.items(), key=lambda item: -item[1])[:20]:
        print("  %-3s -> %-3s %d" % (a, b, count))
    per_phone = {}
    for a, b in ours:
        per_phone.setdefault(a, []).append(a == b)
    print("\nour per-phone accuracy (n>=8):")
    line = []
    for phone, values in sorted(per_phone.items(), key=lambda item: numpy.mean(item[1])):
        if len(values) >= 8:
            line.append("%s %.0f%%" % (phone, 100.0 * numpy.mean(values)))
    print("  " + ", ".join(line))
    # ROWS=dh,v,z prints those phones' full rows (what they were heard as).
    for phone in [p for p in os.environ.get("ROWS", "").split(",") if p]:
        row = {}
        for a, b in ours:
            if a == phone:
                row[b] = row.get(b, 0) + 1
        print("  %-3s -> %s" % (phone, " ".join("%s %d" % item for item in sorted(row.items(), key=lambda item: -item[1]))))
    print("\nPHONE_TEST human %.1f ours %.1f" % (human_accuracy, our_accuracy))


if __name__ == "__main__":
    main()
