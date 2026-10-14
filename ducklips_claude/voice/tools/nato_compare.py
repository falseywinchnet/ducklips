"""Every NATO word, ours against the NATO reference recording, by level
contours in three bands (0.2-0.7, 0.7-2.5, 2.5-7 kHz; 10 ms; dB re
the word's loudest frame). Both are trimmed to speech (the reference: frames
within 35 dB of the peak, 30 ms margin; ours: its phone timeline ± 30 ms) and aligned by dynamic time warping on the three contours, so
duration differences do not count; the mean absolute difference along the
path is reported per band, worst words first. A figure shows each word's
aligned contours (reference solid, ours dashed).
usage: nato_compare.py OUT.png [say arguments...]
ENGINE=tract renders through the segment-module engine (../tract/build/speak);
ENGINE=sheet through the sheet tube (../tract/build/sheet_speak).
REF=owner compares against the owner's own reading (ref/user_nato.wav, two
passes through the alphabet; PASS=2 takes the second)."""
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy
from scipy.signal import butter, lfilter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca

WORDS = ["Alfa", "Bravo", "Charlie", "Delta", "Echo", "Foxtrot", "Golf", "Hotel", "India", "Juliett", "Kilo",
         "Lima", "Mike", "November", "Oscar", "Papa", "Quebec", "Romeo", "Sierra", "Tango", "Uniform", "Victor",
         "Whiskey", "X-ray", "Yankee", "Zulu"]
BANDS = ((200, 700), (700, 2500), (2500, 7000))  # NATO's microphone chain has little below 200 Hz
RATE = 16000
# Word onsets and ends (s) in ref/user_nato.wav, found by level (35 dB under
# the peak, gaps under 150 ms bridged): two passes, Alfa to Zulu.
OWNER_SPANS = (
    ((1.29, 1.98), (2.67, 3.37), (3.93, 4.54), (5.06, 5.74), (6.32, 6.91), (7.61, 8.31), (8.78, 9.42),
     (9.85, 10.53), (11.22, 11.86), (12.22, 12.89), (13.67, 14.30), (14.76, 15.41), (15.95, 16.63),
     (16.98, 17.76), (18.14, 18.79), (19.12, 19.87), (20.29, 20.97), (21.51, 22.27), (22.45, 23.16),
     (23.76, 24.36), (24.83, 25.61), (26.02, 26.61), (27.23, 28.04), (28.45, 29.16), (29.71, 30.24),
     (30.79, 31.46)),
    ((31.94, 32.49), (33.22, 33.67), (34.34, 34.89), (35.62, 36.19), (36.83, 37.33), (37.93, 38.61),
     (39.16, 39.75), (40.46, 40.94), (41.58, 42.23), (42.80, 43.48), (44.24, 44.70), (45.27, 45.85),
     (46.42, 47.03), (47.62, 48.18), (48.93, 49.42), (50.06, 50.53), (51.13, 51.72), (52.52, 53.00),
     (53.62, 54.14), (54.59, 55.22), (55.93, 56.38), (57.10, 57.58), (58.29, 58.76), (59.46, 59.91),
     (60.62, 61.05), (61.73, 62.24)))


def reference_word(index):
    letter = WORDS[index][0].lower()
    if os.environ.get("REF") != "owner":
        reference, rate = ca.read_wave(os.path.join(ca.ROOT, "ref", "nato16", "nato_%s.wav" % letter))
        return ca.resample(reference, rate, RATE)
    span = OWNER_SPANS[int(os.environ.get("PASS", "1")) - 1][index]
    recording, rate = ca.read_wave(os.path.join(ca.ROOT, "ref", "user_nato.wav"))
    piece = recording[int((span[0] - 0.1) * rate):int((span[1] + 0.1) * rate)]
    return ca.resample(piece, rate, RATE)


def contours(data, span=None):
    rows = []
    frame = int(0.01 * RATE)
    for low, high in BANDS:
        b, a = butter(2, [low / (RATE / 2), high / (RATE / 2)], "band")
        x = lfilter(b, a, data)
        rows.append([10 * numpy.log10(numpy.mean(x[i:i + frame] ** 2) + 1e-12) for i in range(0, len(x) - frame, frame)])
    level = numpy.array(rows)
    if span is not None:
        level = level[:, int(span[0] * 100):int(span[1] * 100) + 1]
    else:
        total = 10 * numpy.log10(numpy.sum(10 ** (level / 10), axis=0))
        keep = numpy.where(total > total.max() - 35)[0]
        level = level[:, max(keep[0] - 3, 0):keep[-1] + 4]
    return level - level.max()


def warp(a, b):
    n, m = a.shape[1], b.shape[1]
    cost = numpy.abs(a[:, :, None] - b[:, None, :]).sum(axis=0)
    total = numpy.full((n + 1, m + 1), numpy.inf)
    total[0, 0] = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            total[i, j] = cost[i - 1, j - 1] + min(total[i - 1, j], total[i, j - 1], total[i - 1, j - 1])
    path = []
    i, j = n, m
    while i > 0 and j > 0:
        path.append((i - 1, j - 1))
        step = numpy.argmin([total[i - 1, j - 1], total[i - 1, j], total[i, j - 1]])
        i, j = (i - 1, j - 1) if step == 0 else (i - 1, j) if step == 1 else (i, j - 1)
    return path[::-1]


def main():
    out = sys.argv[1]
    extra = sys.argv[2:]
    results = []
    figure, axes = plt.subplots(6, 5, figsize=(20, 18))
    for index, word in enumerate(WORDS):
        letter = word[0].lower()
        reference = reference_word(index)
        wave = os.path.join("/tmp", "nato_ours_%s.wav" % letter)
        score = wave + ".score"
        if os.environ.get("ENGINE") == "tract":
            # The segment-module engine (../tract) speaks the front end's score.
            subprocess.run([os.path.join(ca.ROOT, "build", "say"), word + ".", wave, "--score", score, "--plan-only",
                            "--final-stop-breath", "--pressure", "600", "--phones",
                            os.path.join(ca.ROOT, "..", "tract", "data", "phones_tract.txt")] + extra,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            subprocess.run([os.path.join(ca.ROOT, "..", "tract", "build", "speak"), score, wave],
                           cwd=os.path.join(ca.ROOT, "..", "tract"), stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, check=True)
        elif os.environ.get("ENGINE") == "sheet":
            # The sheet tube speaks the front end's score (base pitch in its cord's range).
            subprocess.run([os.path.join(ca.ROOT, "build", "say"), word + ".", wave, "--score", score, "--plan-only",
                            "--pitch", os.environ.get("PITCH", "150"), "--phones",
                            os.path.join(ca.ROOT, "..", "tract", "data", "phones_sheet.txt")] + extra,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            subprocess.run([os.path.join(ca.ROOT, "..", "tract", "build", "sheet_speak"), score, wave] +
                           os.environ.get("SHEET", "").split(),
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        else:
            subprocess.run([os.path.join(ca.ROOT, "build", "say"), word + ".", wave, "--score", score] + extra,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        phones = ca.our_phones(score)
        ours, our_rate = ca.read_wave(wave)
        ours = ca.resample(ours, our_rate, RATE)
        ref_level = contours(reference)
        our_level = contours(ours, (phones[0][1] - 0.03, phones[-1][2] + 0.03))
        path = warp(ref_level, our_level)
        ref_path = numpy.array([ref_level[:, i] for i, j in path]).T
        our_path = numpy.array([our_level[:, j] for i, j in path]).T
        error = numpy.mean(numpy.abs(ref_path - our_path), axis=1)
        results.append((error.mean(), word, error, ref_level.shape[1] * 10, our_level.shape[1] * 10))
        axis = axes[index // 5][index % 5]
        for band, color in zip(range(3), ("tab:blue", "tab:green", "tab:red")):
            axis.plot(ref_path[band], color=color, lw=1.2)
            axis.plot(our_path[band], color=color, lw=1.2, ls="--")
        axis.set_title("%s  %.1f dB" % (word, error.mean()), fontsize=9)
        axis.set_ylim(-60, 2)
    for spare in range(len(WORDS), 30):
        axes[spare // 5][spare % 5].axis("off")
    figure.tight_layout()
    figure.savefig(out, dpi=70)
    print("word        mean   low   mid  high   ref ms  ours ms")
    for mean, word, error, ref_ms, our_ms in sorted(results, reverse=True):
        print("%-10s %5.1f %5.1f %5.1f %5.1f   %5d   %5d" % (word, mean, error[0], error[1], error[2], ref_ms, our_ms))
    print("ALL %.2f" % numpy.mean([r[0] for r in results]))


if __name__ == "__main__":
    main()
