"""Level of sibilant frication relative to vowels in a recording.

usage: python3 fricative_level.py FILE.wav [FILE.wav ...]

Sibilant frames: 20 ms frames whose 4-12 kHz energy exceeds four times their
0-2 kHz energy and lies within 40 dB of the loudest frame. Vowel frames:
strongly periodic frames within 30 dB of the loudest. Reports the median
sibilant frame level against the median vowel frame level (dB, RMS), and the
sibilant spectrum's band shares.
"""
import sys

import numpy

exec(open(__file__.replace("fricative_level.py", "voice_compare.py")).read().split("def levinson")[0])


def main():
    for path in sys.argv[1:]:
        signal, rate = read_wav(path)
        frame = int(0.02 * rate)
        frequency = numpy.fft.rfftfreq(frame, 1.0 / rate)
        window = numpy.hanning(frame)
        levels = []
        kinds = []
        spectra = []
        for start in range(0, len(signal) - frame, frame // 2):
            piece = signal[start:start + frame]
            level = 10.0 * numpy.log10(numpy.mean(piece * piece) + 1.0e-20)
            power = numpy.abs(numpy.fft.rfft(piece * window)) ** 2
            low = numpy.sum(power[frequency < 2000.0])
            high = numpy.sum(power[(frequency > 4000.0) & (frequency < 12000.0)])
            f0, strength = periodicity(piece, rate)
            kind = ""
            if high > 4.0 * low:
                kind = "sibilant"
            elif strength > 0.6:
                kind = "vowel"
            levels.append(level)
            kinds.append(kind)
            spectra.append(power)
        levels = numpy.array(levels)
        top = numpy.max(levels)
        vowel = [levels[i] for i in range(len(levels)) if kinds[i] == "vowel" and levels[i] > top - 30.0]
        sibilant_index = [i for i in range(len(levels)) if kinds[i] == "sibilant" and levels[i] > top - 40.0]
        sibilant = [levels[i] for i in sibilant_index]
        print("== %s: vowel frames %d, sibilant frames %d" % (path, len(vowel), len(sibilant)))
        if vowel and sibilant:
            print("  sibilant re vowel: median %.1f dB, loudest sibilant %.1f dB"
                  % (numpy.median(sibilant) - numpy.median(vowel), numpy.max(sibilant) - numpy.median(vowel)))
            mean = numpy.mean([spectra[i] for i in sibilant_index], axis=0)
            total = numpy.sum(mean)
            bands = ((0, 2000), (2000, 4000), (4000, 6000), (6000, 8000), (8000, 12000), (12000, 16000))
            print("  sibilant spectrum share:", " ".join("%d-%dk %.0f%%" % (a / 1000, b / 1000,
                  100.0 * numpy.sum(mean[(frequency >= a) & (frequency < b)]) / total) for a, b in bands))


main()
