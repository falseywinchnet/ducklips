"""Lists impulsive transients: first differences far above their 5 ms local
RMS, in active (non-silent) stretches.

usage: python3 transients.py FILE.wav [FILE.wav ...]
"""
import sys

import numpy

exec(open(__file__.replace("transients.py", "voice_compare.py")).read().split("def periodicity")[0])


def transients(signal, rate):
    difference = numpy.diff(signal)
    window = int(0.005 * rate)
    local = numpy.sqrt(numpy.convolve(difference * difference, numpy.ones(window) / window, mode="same")) + 1.0e-12
    ratio = numpy.abs(difference) / local
    activity = numpy.sqrt(numpy.convolve(signal * signal, numpy.ones(int(0.02 * rate)) / int(0.02 * rate), "same"))
    floor = numpy.percentile(activity, 95) * 0.03
    events = []
    index = 0
    while index < len(ratio):
        if ratio[index] > 7.0 and activity[index] > floor:
            events.append(index / rate)
            index += int(0.003 * rate)
        else:
            index += 1
    return events


for path in sys.argv[1:]:
    signal, rate = read_wav(path)
    events = transients(signal, rate)
    print("%s: %d transients: %s" % (path, len(events), " ".join("%.3f" % time for time in events)))
