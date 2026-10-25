"""Pitch shape of vowel onsets after silence or a voiceless consonant: the
first 30 ms of voicing (F0 every 5 ms) against the vowel's own median F0, in
cents. Prints the median first-frame offset and the median of 0-15 ms and
15-30 ms means, plus counts.
usage: onset.py bdl LABDIR WAVDIR [N] | onset.py ours SCOREDIR WAVDIR"""
import sys, os, glob, numpy
sys.path.insert(0, "/tmp")
from vw import load, pitch_at, VOWELS
from dips import phones_bdl, phones_ours

VOICELESS = {"pau", "sil", "h#", "hh", "p", "t", "k", "f", "th", "s", "sh", "ch", "ssil"}
kind = sys.argv[1]
source = phones_bdl(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 29) if kind == "bdl" \
    else phones_ours(sys.argv[2], sys.argv[3])
first, early, late = [], [], []
for wav, phones in source:
    if not os.path.exists(wav):
        continue
    x = load(wav)
    for a, b in zip(phones, phones[1:]):
        if b[0] not in VOWELS or a[0] not in VOICELESS or b[2] - b[1] < 0.08:
            continue
        track = [pitch_at(x, int((t - 0.015) * 16000))[0] for t in numpy.arange(b[1] - 0.02, b[2] - 0.015, 0.005)]
        voiced = [v for v in track if v > 0]
        if len(voiced) < 8:
            continue
        med = numpy.median(voiced)
        k = next(i for i, v in enumerate(track) if v > 0)
        seg = [1200 * numpy.log2(v / med) for v in track[k:k + 6] if v > 0]
        if len(seg) < 4:
            continue
        first.append(seg[0]); early.append(numpy.mean(seg[:3])); late.append(numpy.mean(seg[3:]))
print("%-12s n %3d  first %+4.0f c  0-15ms %+4.0f  15-30ms %+4.0f   |first| %3.0f" % (
    sys.argv[3].rstrip("/").split("/")[-1], len(first), numpy.median(first), numpy.median(early), numpy.median(late),
    numpy.median(numpy.abs(first))))
