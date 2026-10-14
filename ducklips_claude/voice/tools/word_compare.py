"""The same word in the same ARCTIC sentence, bdl beside ours: renders our
version, finds the word's phone sequence in both timelines and draws their
spectrograms (with 60 ms either side) one above the other.
usage: word_compare.py OUT.png UTTERANCE "PH ON ES" [UTTERANCE "PH ON ES" ...]
       phones in ARCTIC label spelling, e.g. arctic_a0036 "sh iy".
SPEAKER=rms takes the reference from ref/arctic_rms (all 1132 a-set
sentences; bdl here has only a0001-a0060)."""
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca

PAD = 0.06
FOLD = {"ax": "ah", "axr": "er", "ix": "ih", "dx": "t"}


SPEAKER = os.environ.get("SPEAKER", "bdl")
SPEAKER_DIR = os.path.join(ca.ROOT, "ref", "arctic_" + SPEAKER)


def speaker_phones(name):
    phones = []
    start = 0.0
    for line in open(os.path.join(SPEAKER_DIR, "lab", name + ".lab")):
        parts = line.split()
        if len(parts) != 3:
            continue
        end = float(parts[0])
        if parts[2] not in ("pau", "h#", "ssil"):
            phones.append((parts[2], start, end))
        start = end
    return phones


def find_span(phones, sequence):
    names = [FOLD.get(p, p) for p, s, e in phones]
    target = [FOLD.get(p, p) for p in sequence]
    for start in range(len(names) - len(target) + 1):
        if names[start:start + len(target)] == target:
            return phones[start][1], phones[start + len(target) - 1][2]
    return None


def spectrogram(axis, data, rate, start, end, title):
    a0, a1 = int(max(start - PAD, 0) * rate), int((end + PAD) * rate)
    segment = data[a0:a1]
    size, step = int(0.012 * rate), int(0.001 * rate)
    window = numpy.hanning(size)
    frames = [numpy.abs(numpy.fft.rfft(segment[i:i + size] * window)) for i in range(0, len(segment) - size, step)]
    level = 20 * numpy.log10(numpy.array(frames).T + 1e-9)
    level -= level.max()
    axis.imshow(numpy.maximum(level, -70), origin="lower", aspect="auto", cmap="magma",
                extent=[0, len(segment) / rate, 0, rate / 2])
    axis.axvline(PAD, color="w", lw=0.5)
    axis.axvline(PAD + end - start, color="w", lw=0.5)
    axis.set_title(title, fontsize=9)
    axis.set_ylim(0, 8000)


def main():
    out = sys.argv[1]
    pairs = list(zip(sys.argv[2::2], sys.argv[3::2]))
    texts = ca.transcripts()
    figure, axes = plt.subplots(2, len(pairs), figsize=(3.2 * len(pairs), 6), squeeze=False)
    for column, (name, phones_text) in enumerate(pairs):
        sequence = phones_text.split()
        wave48 = os.path.join("/tmp", name + "_word48.wav")
        score = os.path.join("/tmp", name + "_word.score")
        subprocess.run([os.path.join(ca.ROOT, "build", "say"), texts[name], wave48, "--score", score],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        ours, our_rate = ca.read_wave(wave48)
        ours = ca.resample(ours, our_rate, 16000)
        ref, ref_rate = ca.read_wave(os.path.join(SPEAKER_DIR, "wav", name + ".wav"))
        ref_span = find_span(speaker_phones(name), sequence)
        our_span = find_span(ca.our_phones(score), sequence)
        if ref_span is None or our_span is None:
            print("not found:", name, phones_text, ref_span, our_span)
            continue
        spectrogram(axes[0][column], ref, ref_rate, ref_span[0], ref_span[1], "%s %s  %s" % (SPEAKER, name[-4:], phones_text))
        spectrogram(axes[1][column], ours, 16000, our_span[0], our_span[1], "ours  %s" % phones_text)
        print("%s %-10s ref %.0f ms  ours %.0f ms" % (name, phones_text, 1000 * (ref_span[1] - ref_span[0]),
                                                      1000 * (our_span[1] - our_span[0])))
    figure.tight_layout()
    figure.savefig(out, dpi=90)


if __name__ == "__main__":
    main()
