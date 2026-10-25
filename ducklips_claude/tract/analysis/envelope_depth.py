"""How far the voice's loudness falls between pulses: the rattle the owner
hears in the pinch tube at low pitch (2026-10-17, "higher frequency voices,
no rattle"), measured against real voices.

For each held voiced stretch: the sound band-passed 80 Hz - 6 kHz, its
Hilbert envelope smoothed over 0.3 ms, the pulse peaks (at least 0.75 of
the period apart), and per cycle the drop from each peak to the trough
before the next (dB). Prints F0, the median and 90th percentile drop, and
the trough's share of the cycle spent within 6 dB of it.

usage:
  envelope_depth.py wav PATH START END [LABEL]   one stretch
  envelope_depth.py voiced PATH                  the voiced stretches over 0.3 s
  envelope_depth.py arctic LABDIR WAVDIR [N]     ARCTIC vowels over 0.15 s"""
import glob
import os
import sys

import librosa
import numpy
import soundfile
from scipy.signal import butter, find_peaks, hilbert, sosfiltfilt

VOWELS = {"aa", "ae", "ah", "ao", "iy", "uw", "eh", "er", "ow", "ey", "ay"}


def depth(x, rate):
    if len(x) < int(0.08 * rate):
        return None
    band = sosfiltfilt(butter(4, [80, 6000], "bandpass", fs=rate, output="sos"), x)
    f0, voiced, _ = librosa.pyin(band, fmin=50, fmax=400, sr=rate, frame_length=2048)
    f0 = f0[voiced]
    if len(f0) < 3:
        return None
    pitch = float(numpy.median(f0))
    envelope = numpy.abs(hilbert(band))
    width = max(int(0.0003 * rate), 1)
    envelope = numpy.convolve(envelope, numpy.ones(width) / width, "same")
    period = rate / pitch
    peaks, _ = find_peaks(envelope, distance=int(0.75 * period))
    drops = []
    for left, right in zip(peaks[:-1], peaks[1:]):
        if right - left > 1.5 * period:
            continue
        trough = envelope[left:right].min()
        drops.append(20.0 * numpy.log10(envelope[left] / max(trough, 1e-9)))
    if len(drops) < 3:
        return None
    drops = numpy.array(drops)
    return pitch, float(numpy.median(drops)), float(numpy.percentile(drops, 90)), len(drops)


def report(label, x, rate):
    found = depth(x, rate)
    if found is None:
        return None
    pitch, median, high, count = found
    print("%-28s F0 %5.1f  drop between pulses: median %5.1f dB, 90%% %5.1f dB  (%d cycles)" %
          (label, pitch, median, high, count))
    return found


def voiced_stretches(path):
    x, rate = soundfile.read(path)
    if x.ndim > 1:
        x = x.mean(axis=1)
    hop = 512
    f0, voiced, _ = librosa.pyin(x, fmin=50, fmax=400, sr=rate, frame_length=2048, hop_length=hop)
    start = None
    for index, on in enumerate(list(voiced) + [False]):
        if on and start is None:
            start = index
        elif not on and start is not None:
            if (index - start) * hop / rate > 0.3:
                a = start * hop + int(0.05 * rate)
                b = index * hop - int(0.05 * rate)
                report("%s %.2f-%.2f" % (os.path.basename(path), a / rate, b / rate), x[a:b], rate)
            start = None


def arctic(lab_dir, wav_dir, limit):
    by_vowel = {}
    for lab in sorted(glob.glob(os.path.join(lab_dir, "*.lab")))[:limit]:
        wav = os.path.join(wav_dir, os.path.basename(lab).replace(".lab", ".wav"))
        if not os.path.exists(wav):
            continue
        x, rate = soundfile.read(wav)
        previous = 0.0
        for line in open(lab):
            parts = line.split()
            if len(parts) != 3 or parts[0] == "#":
                continue
            end = float(parts[0])
            phone = parts[2].lower()
            if phone in VOWELS and end - previous > 0.15:
                trim = 0.2 * (end - previous)
                found = depth(x[int((previous + trim) * rate):int((end - trim) * rate)], rate)
                if found is not None:
                    by_vowel.setdefault(phone, []).append(found)
            previous = end
    everything = []
    for phone in sorted(by_vowel):
        rows = numpy.array([row[:3] for row in by_vowel[phone]])
        everything.append(rows)
        print("%-4s n %3d  F0 %5.1f  median drop %5.1f dB  90%% %5.1f dB" %
              (phone, len(rows), rows[:, 0].mean(), numpy.median(rows[:, 1]), numpy.median(rows[:, 2])))
    rows = numpy.concatenate(everything)
    for low, high in ((50, 95), (95, 110), (110, 130), (130, 400)):
        pick = rows[(rows[:, 0] >= low) & (rows[:, 0] < high)]
        if len(pick):
            print("F0 %3d-%3d  n %3d  median drop %5.1f dB" % (low, high, len(pick), numpy.median(pick[:, 1])))


def main():
    mode = sys.argv[1]
    if mode == "wav":
        x, rate = soundfile.read(sys.argv[2])
        a, b = float(sys.argv[3]), float(sys.argv[4])
        label = sys.argv[5] if len(sys.argv) > 5 else os.path.basename(sys.argv[2])
        report(label, x[int(a * rate):int(b * rate)], rate)
    elif mode == "voiced":
        voiced_stretches(sys.argv[2])
    elif mode == "arctic":
        arctic(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 200)


if __name__ == "__main__":
    main()
