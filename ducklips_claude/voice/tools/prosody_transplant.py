"""Renders an ARCTIC sentence with a human reading's prosody: each of our
phones takes the duration of the human phone it aligns to (edit distance on
phone names), each vowel the human's F0 at its middle, and each phrase break
the human's pause there. Our sound source and tract stay ours, so this
separates prosody from voice quality by ear.

Writes OUT_planned.wav (our own prosody) and OUT_human_prosody.wav, 16 kHz.
usage: prosody_transplant.py OUT UTTERANCE   (SPEAKER=rms by default;
ENGINE=tract renders through the segment-module body)"""
import os
import subprocess
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca

SILENCE = ("pau", "h#", "ssil")


def human_phones(directory, name):
    """(phone, start, end) including pauses."""
    phones = []
    start = 0.0
    for line in open(os.path.join(directory, "lab", name + ".lab")):
        parts = line.split()
        if len(parts) != 3:
            continue
        end = float(parts[0])
        phones.append((parts[2], start, end))
        start = end
    return phones


def pitch_at(data, rate, time):
    """Autocorrelation F0 over 40 ms around time, or 0 if unvoiced."""
    size = int(0.04 * rate)
    start = max(int(time * rate) - size // 2, 0)
    segment = data[start:start + size]
    if len(segment) < size:
        return 0.0
    segment = segment - segment.mean()
    correlation = numpy.correlate(segment, segment, "full")[size - 1:]
    low, high = int(rate / 300), int(rate / 60)
    lag = low + int(numpy.argmax(correlation[low:high]))
    if correlation[0] <= 0.0 or correlation[lag] / correlation[0] < 0.45:
        return 0.0
    return rate / lag


def phrase_starts(score):
    """Indexes (into our phone list) of each phrase's first phone."""
    marks = []
    phone_starts = []
    for line in open(score):
        parts = line.split()
        if line.startswith("# ") and len(parts) == 5 and parts[1] != "phones:":
            phone_starts.append(float(parts[2]))
        elif len(parts) >= 3 and parts[1] == "mark" and parts[2] == "phrase":
            marks.append(float(parts[0]))
    return [int(numpy.argmin(numpy.abs(numpy.array(phone_starts) - mark))) for mark in marks]


def main():
    out, name = sys.argv[1], sys.argv[2]
    extra = sys.argv[3:]
    directory = os.path.join(ca.ROOT, "ref", "arctic_" + os.environ.get("SPEAKER", "rms"))
    text = ca.transcripts()[name]
    say = os.path.join(ca.ROOT, "build", "say")
    planned = "/tmp/" + name + "_planned48.wav"
    score = planned + ".score"
    tract = os.environ.get("ENGINE") == "tract"
    speak = os.path.join(ca.ROOT, "..", "tract", "build", "speak")
    if tract:
        # The segment-module body speaks the front end's score.
        extra = ["--final-stop-breath", "--pressure", "600", "--phones",
                 os.path.join(ca.ROOT, "..", "tract", "data", "phones_tract.txt")] + extra
    subprocess.run([say, text, planned, "--score", score] + extra + (["--plan-only"] if tract else []),
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    if tract:
        subprocess.run([speak, score, planned], cwd=os.path.join(ca.ROOT, "..", "tract"), stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=True)
    ours = ca.our_phones(score)
    human_all = human_phones(directory, name)
    human = [phone for phone in human_all if phone[0] not in SILENCE]
    recording, rate = ca.read_wave(os.path.join(directory, "wav", name + ".wav"))
    pairs = dict((j, i) for i, j in ca.align(human, ours))
    lines = []
    for index, (phone, start, end) in enumerate(ours):
        if index not in pairs:
            lines.append("0 0")
            continue
        match = human[pairs[index]]
        pitch = pitch_at(recording, rate, 0.5 * (match[1] + match[2])) if phone in ca.VOWELS else 0.0
        lines.append("%.4f %.1f" % (match[2] - match[1], pitch))
    # Pauses: the human silence (if any) just before the phone that begins
    # each of our later phrases, and the ending silence after the last.
    starts = phrase_starts(score)
    for begin in starts[1:]:
        pause = 0.05
        if begin in pairs:
            match = human[pairs[begin]]
            position = human_all.index(match)
            if position > 0 and human_all[position - 1][0] in SILENCE:
                pause = human_all[position - 1][2] - human_all[position - 1][1]
        lines.append("pause %.4f" % pause)
    timing = "/tmp/" + name + ".timing"
    with open(timing, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    transplanted = "/tmp/" + name + "_human48.wav"
    transplanted_score = transplanted + ".score"
    subprocess.run([say, text, transplanted, "--timing", timing, "--score", transplanted_score] + extra
                   + (["--plan-only"] if tract else []), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   check=True)
    if tract:
        subprocess.run([speak, transplanted_score, transplanted], cwd=os.path.join(ca.ROOT, "..", "tract"),
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    for suffix, path in (("planned", planned), ("human_prosody", transplanted)):
        data, wave_rate = ca.read_wave(path)
        data = ca.resample(data, wave_rate, 16000)
        ca.write_wave(out + "_" + suffix + ".wav", 0.5 * data / (numpy.abs(data).max() + 1e-12), 16000)
    print("matched %d of %d phones; timing in %s" % (len(pairs), len(ours), timing))


if __name__ == "__main__":
    main()
