"""Spectrograms of NATO words: the NATO reference, the owner's reading
(ref/user_nato.wav, first pass) and ours, one row per word, each trimmed as
in nato_compare.py.
usage: nato_views.py OUT.png WORD [WORD ...] [-- say arguments...]"""
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
import nato_compare as nc


def trimmed(data, span=None):
    frame = int(0.01 * nc.RATE)
    level = numpy.array([10 * numpy.log10(numpy.mean(data[i:i + frame] ** 2) + 1e-12)
                         for i in range(0, len(data) - frame, frame)])
    if span is None:
        keep = numpy.where(level > level.max() - 35)[0]
        first, last = max(keep[0] - 3, 0), keep[-1] + 4
    else:
        first, last = int(span[0] * 100), int(span[1] * 100) + 1
    return data[first * frame:last * frame]


def draw(axis, data, title):
    size, step = int(0.016 * nc.RATE), int(0.002 * nc.RATE)
    window = numpy.hanning(size)
    frames = [numpy.abs(numpy.fft.rfft(data[i:i + size] * window)) for i in range(0, len(data) - size, step)]
    level = 20 * numpy.log10(numpy.array(frames).T + 1e-9)
    level -= level.max()
    axis.imshow(numpy.maximum(level, -70), origin="lower", aspect="auto", cmap="magma",
                extent=[0, len(data) / nc.RATE, 0, nc.RATE / 2])
    axis.set_title(title, fontsize=9)


def render(word, extra):
    wave = os.path.join("/tmp", "nato_view_%s.wav" % word[0].lower())
    score = wave + ".score"
    subprocess.run([os.path.join(ca.ROOT, "build", "say"), word + ".", wave, "--score", score] + extra,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    phones = ca.our_phones(score)
    ours, rate = ca.read_wave(wave)
    return trimmed(ca.resample(ours, rate, nc.RATE), (phones[0][1] - 0.03, phones[-1][2] + 0.15))


def main():
    out = sys.argv[1]
    arguments = sys.argv[2:]
    extra = []
    if "--" in arguments:
        extra = arguments[arguments.index("--") + 1:]
        arguments = arguments[:arguments.index("--")]
    figure, axes = plt.subplots(len(arguments), 3, figsize=(13, 2.6 * len(arguments)), squeeze=False)
    for row, word in enumerate(arguments):
        index = nc.WORDS.index(word)
        os.environ.pop("REF", None)
        draw(axes[row][0], trimmed(nc.reference_word(index)), "NATO " + word)
        os.environ["REF"] = "owner"
        draw(axes[row][1], trimmed(nc.reference_word(index)), "owner " + word)
        draw(axes[row][2], render(word, extra), "ours " + word)
    os.environ.pop("REF", None)
    figure.tight_layout()
    figure.savefig(out, dpi=80)


if __name__ == "__main__":
    main()
