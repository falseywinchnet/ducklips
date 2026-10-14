"""Short-transient strength: the loudest 3 ms windows relative to the median
vowel level, for comparing release bursts between recordings.

usage: python3 burst_level.py FILE.wav [FILE.wav ...]
"""
import sys

import numpy

exec(open(__file__.replace("burst_level.py", "voice_compare.py")).read().split("def levinson")[0])


def main():
    for path in sys.argv[1:]:
        signal, rate = read_wav(path)
        short = int(0.003 * rate)
        long = int(0.02 * rate)
        energy = numpy.convolve(signal * signal, numpy.ones(short) / short, mode="same")
        context = numpy.convolve(signal * signal, numpy.ones(long * 5) / (long * 5), mode="same")
        # Vowel reference: median of the 20 ms levels in the loudest half.
        levels = numpy.convolve(signal * signal, numpy.ones(long) / long, mode="same")[::long]
        levels = numpy.sort(levels[levels > 0.0])
        vowel = numpy.median(levels[len(levels) // 2:])
        # Transients: short energy far above its 100 ms context.
        ratio = energy / (context + 1.0e-20)
        events = []
        index = 0
        while index < len(signal):
            if ratio[index] > 6.0 and energy[index] > 0.01 * vowel:
                end = min(len(signal), index + int(0.02 * rate))
                peak = index + int(numpy.argmax(energy[index:end]))
                events.append((peak / rate, 10.0 * numpy.log10(energy[peak] / vowel)))
                index = end
            else:
                index += 1
        print("== %s: %d transients; strongest re vowel level:" % (path, len(events)))
        strongest = sorted(events, key=lambda item: -item[1])[:8]
        print("  " + "  ".join("%.3fs %+.1fdB" % item for item in sorted(strongest)))
        if events:
            print("  median transient %+.1f dB" % numpy.median([item[1] for item in events]))


main()
