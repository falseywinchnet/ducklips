"""Band signal-to-noise per NATO word: ours against the NATO reading and the
owner's. Each word is rendered by the main line (say with the organ
postures, then speak with the fold tube; extra environment passes through)
into OUTDIR. Every glottal cycle in the loud stretches (within 20 dB of the
word's peak; found by waveform matching, none dropped for irregularity) is
compared with the next band by band: the normalized correlation r gives
HNR = 10 log10(r / (1 - r)). Medians per word and band, the cycle count and
the share of cycles whose broadband match passes 0.9 (clean); the table
lists ours minus each reference. (Gated to clean cycles only, as
bandhnr.py is, ours scored 2-10 dB over both references: the gate kept
only our best cycles.)
usage: nato_snr.py OUTDIR [render|measure]   (default: both)"""
import json
import os
import subprocess
import sys

import librosa
import numpy
from scipy.signal import butter, sosfiltfilt

sys.path.insert(0, os.path.expanduser("~/ducklips/voice/tools"))
import nato_compare  # noqa: E402

ROOT = os.path.expanduser("~/ducklips")
BANDS = ((80, 1000), (1000, 2000), (2000, 4000), (4000, 7000))
RATE = 16000
SEED = "73021"


def render_word(word, outdir):
    score = os.path.join(outdir, word + ".score")
    wave = os.path.join(outdir, word + ".wav")
    subprocess.run([os.path.join(ROOT, "voice", "build", "say"), word + ".", "/dev/null", "--score", score,
                    "--plan-only", "--pressure", "750", "--phones", os.path.join(ROOT, "tract", "data", "phones_tract.txt"),
                    "--organs", os.path.join(ROOT, "tract", "data", "organs.txt")],
                   cwd=os.path.join(ROOT, "voice"), check=True, stdout=subprocess.DEVNULL)
    environment = dict(os.environ)
    environment.setdefault("FOLDTUBE", "1")
    subprocess.run([os.path.join(ROOT, "tract", "build", "speak"), score, wave, SEED], env=environment, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def band_filters():
    filters = []
    for low, high in BANDS:
        filters.append(butter(4, [low, min(high, 7900)], btype="band", fs=RATE, output="sos"))
    return filters


def period_guess(matched, at):
    """The period (samples) by autocorrelation of 40 ms from `at`, 60-320 Hz,
    and its normalized peak."""
    piece = matched[at:at + int(0.04 * RATE)]
    piece = piece - piece.mean()
    correlation = numpy.correlate(piece, piece, "full")[len(piece) - 1:]
    shortest = RATE // 320
    longest = min(RATE // 60, len(piece) - 1)
    if correlation[0] <= 0.0 or longest <= shortest:
        return 0, 0.0
    lag = shortest + int(numpy.argmax(correlation[shortest:longest]))
    return lag, correlation[lag] / correlation[0]


def loud_cycles(signal):
    """Every glottal cycle in the loud stretches (within 20 dB of the peak):
    from a period guessed by autocorrelation, each next cycle starts where one
    period best matches the next (0.85-1.15 of the guess, parabolic
    refinement), on the 60 Hz-4 kHz signal. Returns (start, period, match)
    in samples; no cycle is dropped for being irregular."""
    matched = sosfiltfilt(butter(4, [60, 4000], btype="band", fs=RATE, output="sos"), signal)
    window = numpy.ones(320) / 320.0
    level = numpy.sqrt(numpy.convolve(signal ** 2, window, "same"))
    loud = level > level.max() * 10.0 ** (-20.0 / 20.0)
    found = []
    at = 0
    guess = 0.0
    hop = int(0.005 * RATE)
    while at + int(0.05 * RATE) < len(signal):
        if not loud[at]:
            at += hop
            guess = 0.0
            continue
        if guess <= 0.0:
            lag, strength = period_guess(matched, at)
            if strength < 0.3:
                at += hop
                continue
            guess = float(lag)
        size = int(round(guess))
        lowest = int(0.85 * guess)
        highest = int(1.15 * guess) + 1
        if at + highest + size + 1 >= len(signal):
            break
        now = matched[at:at + size]
        energy = numpy.dot(now, now)
        scores = []
        for shift in range(lowest - 1, highest + 2):
            then = matched[at + shift:at + shift + size]
            scores.append(numpy.dot(now, then) / numpy.sqrt(energy * numpy.dot(then, then) + 1e-30))
        scores = numpy.array(scores)
        best = int(numpy.argmax(scores[1:-1])) + 1
        curvature = scores[best - 1] - 2.0 * scores[best] + scores[best + 1]
        offset = 0.5 * (scores[best - 1] - scores[best + 1]) / curvature if curvature < 0.0 else 0.0
        period = lowest - 1 + best + offset
        found.append((at, period, scores[best]))
        if scores[best] < 0.3:
            guess = 0.0
        else:
            guess = period
        at += int(round(period))
    return found


def band_hnr(signal, filters):
    """Median cycle-coherence HNR per band (dB) over every loud cycle, the
    number of cycles, and the share whose broadband match passes 0.9."""
    found = loud_cycles(signal)
    medians = []
    for sos in filters:
        banded = sosfiltfilt(sos, signal)
        values = []
        for start, period, _ in found:
            size = int(round(period))
            following = start + size
            if following + size > len(signal):
                continue
            now = banded[start:following]
            then = banded[following:following + size]
            correlation = numpy.dot(now, then) / numpy.sqrt(numpy.dot(now, now) * numpy.dot(then, then) + 1e-20)
            correlation = min(max(correlation, 1e-3), 0.9999)
            values.append(10.0 * numpy.log10(correlation / (1.0 - correlation)))
        medians.append(float(numpy.median(values)) if values else float("nan"))
    clean = float(numpy.mean([match > 0.9 for _, _, match in found])) if found else float("nan")
    return medians, len(found), clean


def reference_nato(index):
    letter = nato_compare.WORDS[index][0].lower()
    signal, _ = librosa.load(os.path.join(ROOT, "voice", "ref", "nato", "nato_%s.wav" % letter), sr=RATE)
    return signal


def reference_owner(index, passage):
    span = nato_compare.OWNER_SPANS[passage][index]
    signal, _ = librosa.load(os.path.join(ROOT, "voice", "ref", "user_nato.wav"), sr=RATE,
                             offset=span[0] - 0.05, duration=span[1] - span[0] + 0.1)
    return signal


def measure(outdir):
    filters = band_filters()
    table = {}
    for index, word in enumerate(nato_compare.WORDS):
        ours, _ = librosa.load(os.path.join(outdir, word + ".wav"), sr=RATE)
        row = {"ours": band_hnr(ours, filters), "nato": band_hnr(reference_nato(index), filters)}
        owner_first = band_hnr(reference_owner(index, 0), filters)
        owner_second = band_hnr(reference_owner(index, 1), filters)
        row["owner"] = ([float(numpy.nanmean([a, b])) for a, b in zip(owner_first[0], owner_second[0])],
                        owner_first[1] + owner_second[1], float(numpy.nanmean([owner_first[2], owner_second[2]])))
        table[word] = row
    with open(os.path.join(outdir, "snr.json"), "w") as handle:
        json.dump(table, handle)
    print("band HNR dB (0.08-1 / 1-2 / 2-4 / 4-7 kHz), [cycles, clean share]")
    print("%-9s %-36s %-36s %-36s" % ("word", "ours", "NATO", "owner"))
    for word in nato_compare.WORDS:
        row = table[word]
        cells = []
        for key in ("ours", "nato", "owner"):
            values, count, clean = row[key]
            cells.append(" ".join("%5.1f" % value for value in values) + " [%3d %3.0f%%]" % (count, 100.0 * clean))
        print("%-9s %-36s %-36s %-36s" % (word, cells[0], cells[1], cells[2]))
    for key in ("nato", "owner"):
        ours = numpy.array([table[word]["ours"][0] for word in nato_compare.WORDS])
        theirs = numpy.array([table[word][key][0] for word in nato_compare.WORDS])
        print("median ours - %-5s: %s" % (key, "  ".join("%+5.1f" % value for value in numpy.nanmedian(ours - theirs, 0))))
    for key in ("ours", "nato", "owner"):
        values = numpy.array([table[word][key][0] for word in nato_compare.WORDS])
        clean = numpy.nanmedian([table[word][key][2] for word in nato_compare.WORDS])
        print("median %-5s       : %s   clean %3.0f%%" % (key, "  ".join("%5.1f" % value for value in numpy.nanmedian(values, 0)),
                                                         100.0 * clean))


def main():
    outdir = sys.argv[1]
    step = sys.argv[2] if len(sys.argv) > 2 else "both"
    os.makedirs(outdir, exist_ok=True)
    if step in ("render", "both"):
        for word in nato_compare.WORDS:
            render_word(word, outdir)
    if step in ("measure", "both"):
        measure(outdir)


if __name__ == "__main__":
    main()
