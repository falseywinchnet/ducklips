"""Learns vowel pitch from a human reading instead of hand rules.

For every ARCTIC sentence the planner's phones (say --plan-only --features)
are aligned to the human's labelled phones; each aligned vowel's F0 at its
middle is taken from the recording (autocorrelation, 40 ms) and expressed in
semitones re the sentence's median. It is fitted as a sum of factors (stress,
accent, nuclear accent, phrase-final boundary type, position in the phrase
and sentence, syllable position in the word, the neighbouring consonants'
voicing) by Huber-weighted least squares, holding one sentence in five out.
The planner's own vowel pitches are scored on the same held-out vowels.

Writes "factor NAME VALUE" lines (semitones).
usage: pitch_fit.py OUT.txt   (SPEAKER=rms by default)"""
import os
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
import duration_fit as df

VOWELS = df.VOWELS


def f0_at(data, rate, time):
    size = int(0.04 * rate)
    start = int(time * rate) - size // 2
    segment = data[max(start, 0):start + size]
    if len(segment) < size:
        return numpy.nan
    segment = segment - segment.mean()
    correlation = numpy.correlate(segment, segment, "full")[size - 1:]
    if correlation[0] <= 0.0:
        return numpy.nan
    low, high = int(rate / 300), int(rate / 60)
    lag = low + int(numpy.argmax(correlation[low:high]))
    if correlation[lag] / correlation[0] < 0.5:
        return numpy.nan
    # Parabolic refinement of the lag.
    if 0 < lag < len(correlation) - 1:
        a, b, c = correlation[lag - 1], correlation[lag], correlation[lag + 1]
        denominator = a - 2.0 * b + c
        if denominator != 0.0:
            lag = lag + 0.5 * (a - c) / denominator
    return rate / lag


def vowel_factors(rows, index):
    """Named pitch factors for one planned vowel."""
    row = rows[index]
    factors = ["stress:%d" % row["stress"]]
    if row["nuclear"]:
        factors.append("nuclear")
    elif row["accent"]:
        factors.append("accent")
    # Phrase extent from the planner's own pauses.
    starts = [0]
    for k in range(len(rows) - 1):
        if rows[k + 1]["start"] - (rows[k]["start"] + rows[k]["duration"]) > 0.05:
            starts.append(k + 1)
    phrase = max(s for s in starts if s <= index)
    phrase_end = min([s for s in starts if s > index] + [len(rows)])
    vowels_in_phrase = [k for k in range(phrase, phrase_end) if rows[k]["name"] in VOWELS]
    position = vowels_in_phrase.index(index) / max(len(vowels_in_phrase) - 1, 1)
    factors.append("phrase_progress:%d" % min(int(position * 4), 3))
    factors.append("phrase_number:%d" % min(starts.index(phrase), 2))
    last_phrase = phrase_end == len(rows)
    if row["phrase_final"]:
        factors.append("final:" + ("end" if last_phrase else "comma") + (":nuclear" if row["nuclear"] else ""))
    if vowels_in_phrase and vowels_in_phrase[0] == index:
        factors.append("phrase_first_vowel")
    before = rows[index - 1]["name"] if index > 0 else "#"
    factors.append("after:" + (df.consonant_class(before) if before != "#" else "start"))
    return factors


def sentence_records(item):
    name, text, directory = item
    try:
        rows = df.planned(name, text)
        reference = df.human(directory, name)
        recording, rate = ca.read_wave(os.path.join(directory, "wav", name + ".wav"))
    except (OSError, Exception):
        return []
    ours = [(row["name"], row["start"], row["start"] + row["duration"]) for row in rows]
    pairs = []
    for i, j in ca.align(reference, ours):
        if rows[j]["name"] not in VOWELS or reference[i][2] - reference[i][1] < 0.04:
            continue
        hz = f0_at(recording, rate, 0.5 * (reference[i][1] + reference[i][2]))
        if numpy.isfinite(hz) and rows[j]["pitch_hz"] > 0.0:
            pairs.append((j, 12.0 * numpy.log2(hz), 12.0 * numpy.log2(rows[j]["pitch_hz"])))
    if len(pairs) < 4:
        return []
    human_median = numpy.median([p[1] for p in pairs])
    our_median = numpy.median([p[2] for p in pairs])
    return [(name, vowel_factors(rows, j), human - human_median, mine - our_median) for j, human, mine in pairs]


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

    def matrix(rows):
        result = numpy.zeros((len(rows), len(names) + 1))
        result[:, -1] = 1.0
        for row, record in enumerate(rows):
            for factor in record[1]:
                if factor in names:
                    result[row, names[factor]] = 1.0
        return result

    solution = df.huber_fit(matrix(train), numpy.array([r[2] for r in train]), ridge=0.5)
    learned = matrix(test) @ solution
    target = numpy.array([r[2] for r in test])
    current = numpy.array([r[3] for r in test])
    flat = numpy.zeros(len(test))
    for label, guess in (("learned", learned), ("planner now", current), ("flat (no movement)", flat)):
        error = target - guess
        print("held-out %-19s rms %.2f st  median |e| %.2f st  correlation %.2f" % (
            label, numpy.sqrt(numpy.mean(error ** 2)), numpy.median(numpy.abs(error)),
            numpy.corrcoef(target, guess)[0, 1] if numpy.std(guess) > 0 else 0.0))
    with open(out, "w") as handle:
        handle.write("# learned vowel pitch factors, semitones re the sentence median\n")
        handle.write("intercept %.4f\n" % solution[-1])
        for factor, k in sorted(names.items()):
            handle.write("factor %s %.4f\n" % (factor, solution[k]))
    for value, factor in sorted(((solution[k], f) for f, k in names.items()), reverse=True):
        print("  %-28s %+.2f st" % (factor, value))


if __name__ == "__main__":
    main()
