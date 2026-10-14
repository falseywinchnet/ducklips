"""Pop rate: short (2 ms) bursts more than 2.2 times the surrounding 30 ms
level, in audible stretches, for our ARCTIC renders and the human recordings.
usage: pops.py [FIRST LAST]"""
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca


def pops(data, rate):
    short_size = int(0.002 * rate)
    long_size = int(0.03 * rate)
    short = numpy.sqrt(numpy.convolve(data ** 2, numpy.ones(short_size) / short_size, "same")) + 1e-9
    level = numpy.sqrt(numpy.convolve(data ** 2, numpy.ones(long_size) / long_size, "same")) + 1e-9
    hits = numpy.flatnonzero((short / level > 2.2) & (level > 0.05 * level.max()))
    events = 0
    last = -10 ** 9
    for index in hits:
        if index - last > int(0.02 * rate):
            events += 1
        last = index
    return events


if __name__ == "__main__":
    first = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    last = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    for side in ("ours", "human"):
        total = 0
        seconds = 0.0
        for number in range(first, last + 1):
            name = "arctic_a%04d" % number
            path = os.path.join(ca.OUT, name + "_ours48.wav") if side == "ours" else os.path.join(ca.REF, "wav", name + ".wav")
            data, rate = ca.read_wave(path)
            total += pops(data, rate)
            seconds += len(data) / rate
        print("%-5s pops: %d in %.0f s (%.2f per second)" % (side, total, seconds, total / seconds))
