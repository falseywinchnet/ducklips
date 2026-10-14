"""Learns phone durations from a human reading instead of hand rules.

For every ARCTIC sentence the planner's own phones (with stress, accent,
word and phrase position: say --plan-only --features) are aligned to the
human's labelled phones; log duration is fitted as a sum of factors
(phone identity plus context) with a per-sentence tempo offset, by
Huber-weighted least squares. Human labels under 20 ms (automatic-label
fragments) are dropped. One sentence in five is held out, and the learned
model and the planner's current durations are scored there on the same
footing (log error after removing each sentence's mean tempo).

Writes the factors as "factor NAME VALUE" lines (natural-log units, VALUE
added to log seconds) for the planner to read.
usage: duration_fit.py OUT.txt   (SPEAKER=rms by default)"""
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca

SILENCE = ("pau", "h#", "ssil")
VOWELS = set(ca.VOWELS)
VOICED_OBSTRUENTS = {"b", "d", "g", "v", "dh", "z", "zh", "jh"}
VOICELESS_OBSTRUENTS = {"p", "t", "k", "f", "th", "s", "sh", "ch", "hh"}


def planned(name, text):
    features = "/tmp/duration_fit_%s.txt" % name
    subprocess.run([os.path.join(ca.ROOT, "build", "say"), text, "/tmp/unused.wav", "--plan-only", "--features",
                    features], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    rows = []
    for line in open(features):
        if line.startswith("#"):
            continue
        parts = line.split()
        rows.append({"name": parts[0].lower(), "stress": int(parts[2]), "start": float(parts[3]),
                     "duration": float(parts[4]), "pitch_hz": float(parts[5]), "word_start": parts[6] == "1", "word_end": parts[7] == "1",
                     "accent": parts[8] == "1", "nuclear": parts[9] == "1", "phrase_final": parts[10] == "1",
                     "word": parts[11]})
    return rows


def human(directory, name):
    phones = []
    start = 0.0
    for line in open(os.path.join(directory, "lab", name + ".lab")):
        parts = line.split()
        if len(parts) != 3:
            continue
        end = float(parts[0])
        if parts[2] not in SILENCE:
            phones.append((parts[2], start, end))
        start = end
    return phones


def consonant_class(name):
    if name in VOICED_OBSTRUENTS:
        return "voiced"
    if name in VOICELESS_OBSTRUENTS:
        return "voiceless"
    if name in VOWELS:
        return "vowel"
    return "sonorant"


def context(rows, index):
    """Named factors for one planned phone (the planner knows all of them)."""
    row = rows[index]
    name = row["name"]
    factors = ["phone:" + name]
    vowel = name in VOWELS
    # Boundary: a pause after this phone's syllable, and which kind.
    ends_phrase = index + 1 == len(rows) or rows[index + 1]["start"] - (row["start"] + row["duration"]) > 0.05
    last_phrase = all(not (rows[later + 1]["start"] - (rows[later]["start"] + rows[later]["duration"]) > 0.05)
                      for later in range(index, len(rows) - 1))
    if row["phrase_final"]:
        boundary = "end" if last_phrase else "comma"
        factors.append(("final_vowel:" if vowel else "final_consonant:") + boundary)
        if ends_phrase and not vowel:
            factors.append("last_consonant:" + boundary)
    # Word shape.
    first = index
    while first > 0 and not rows[first]["word_start"]:
        first -= 1
    last = index
    while last + 1 < len(rows) and not rows[last]["word_end"]:
        last += 1
    syllables = sum(1 for k in range(first, last + 1) if rows[k]["name"] in VOWELS)
    factors.append("syllables:%d" % min(syllables, 3))
    if vowel:
        factors.append("stress:%d" % row["stress"])
        if row["nuclear"]:
            factors.append("nuclear")
        elif row["accent"]:
            factors.append("accent")
        later_vowel = any(rows[k]["name"] in VOWELS for k in range(index + 1, last + 1))
        earlier_vowel = any(rows[k]["name"] in VOWELS for k in range(first, index))
        factors.append("syllable_position:" + ("only" if not later_vowel and not earlier_vowel else
                                              "first" if not earlier_vowel else "last" if not later_vowel else "middle"))
        following = rows[index + 1]["name"] if index < last else "#"
        factors.append("vowel_before:" + (consonant_class(following) if following != "#" else "word_end"))
    else:
        position = "initial" if row["word_start"] else "final" if row["word_end"] else "medial"
        factors.append("consonant_position:" + position)
        before = index > first and rows[index - 1]["name"] not in VOWELS
        after = index < last and rows[index + 1]["name"] not in VOWELS
        if before or after:
            factors.append("cluster")
        if index + 1 < len(rows) and rows[index + 1]["name"] in VOWELS and rows[index + 1]["stress"] == 1:
            factors.append("before_stress")
    return factors


def sentence_records(item):
    name, text, directory = item
    try:
        rows = planned(name, text)
        reference = human(directory, name)
    except (subprocess.CalledProcessError, OSError):
        return []
    ours = [(row["name"], row["start"], row["start"] + row["duration"]) for row in rows]
    records = []
    for i, j in ca.align(reference, ours):
        seconds = reference[i][2] - reference[i][1]
        if 0.02 <= seconds <= 0.5:
            records.append((name, context(rows, j), numpy.log(seconds), numpy.log(rows[j]["duration"])))
    return records


def design(records, names, sentences):
    matrix = numpy.zeros((len(records), len(names) + len(sentences)))
    for row, (sentence, factors, target, planned_log) in enumerate(records):
        for factor in factors:
            if factor in names:
                matrix[row, names[factor]] = 1.0
        matrix[row, len(names) + sentences[sentence]] = 1.0
    return matrix


def huber_fit(matrix, target, ridge=1.0, rounds=6):
    weights = numpy.ones(len(target))
    regulariser = ridge * numpy.eye(matrix.shape[1])
    for _ in range(rounds):
        weighted = matrix * weights[:, None]
        solution = numpy.linalg.solve(matrix.T @ weighted + regulariser, weighted.T @ target)
        residual = target - matrix @ solution
        scale = 1.4826 * numpy.median(numpy.abs(residual)) + 1e-9
        weights = numpy.minimum(1.0, 1.345 * scale / (numpy.abs(residual) + 1e-12))
    return solution


def centred_error(records, predicted):
    """RMS log error after removing each sentence's mean (tempo)."""
    by_sentence = {}
    for (sentence, factors, target, planned_log), guess in zip(records, predicted):
        by_sentence.setdefault(sentence, []).append(target - guess)
    errors = numpy.concatenate([numpy.array(v) - numpy.mean(v) for v in by_sentence.values()])
    return numpy.sqrt(numpy.mean(errors ** 2)), numpy.median(numpy.abs(errors))


def main():
    out = sys.argv[1]
    directory = os.path.join(ca.ROOT, "ref", "arctic_" + os.environ.get("SPEAKER", "rms"))
    texts = ca.transcripts()
    items = [(name, texts[name], directory) for name in sorted(texts)
             if os.path.exists(os.path.join(directory, "lab", name + ".lab"))]
    with ThreadPoolExecutor(8) as pool:
        records = [record for batch in pool.map(sentence_records, items) for record in batch]
    held = {name for index, (name, text, d) in enumerate(items) if index % 5 == 0}
    train = [r for r in records if r[0] not in held]
    test = [r for r in records if r[0] in held]
    counts = {}
    for record in train:
        for factor in record[1]:
            counts[factor] = counts.get(factor, 0) + 1
    names = {factor: k for k, factor in enumerate(sorted(f for f, c in counts.items() if c >= 8))}
    sentences = {s: k for k, s in enumerate(sorted({r[0] for r in records}))}
    solution = huber_fit(design(train, names, sentences), numpy.array([r[2] for r in train]))
    test_matrix = design(test, names, sentences)
    test_matrix[:, len(names):] = 0.0  # tempo removed by centring
    learned = test_matrix @ solution
    current = numpy.array([r[3] for r in test])
    print("records %d train, %d held out; factors %d" % (len(train), len(test), len(names)))
    print("held-out log error (rms, median |.|): learned %.3f %.3f   planner now %.3f %.3f" % (
        centred_error(test, learned) + centred_error(test, current)))
    tempo = numpy.median(solution[len(names):])
    with open(out, "w") as handle:
        handle.write("# learned log-duration factors (natural log; sum + tempo = log seconds)\n")
        handle.write("tempo %.5f\n" % tempo)
        for factor, k in sorted(names.items(), key=lambda item: item[0]):
            handle.write("factor %s %.5f\n" % (factor, solution[k]))
    contexts = sorted(((solution[k], f) for f, k in names.items() if not f.startswith("phone:")), reverse=True)
    for value, factor in contexts:
        print("  %-34s x%.2f" % (factor, numpy.exp(value)))


if __name__ == "__main__":
    main()
