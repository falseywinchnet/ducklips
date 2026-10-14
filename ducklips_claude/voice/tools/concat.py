"""Joins WAV files (same rate) after trimming leading/trailing silence.
usage: concat.py OUT.wav IN1.wav [IN2.wav ...]"""
import sys
import wave

import numpy

pieces = []
rate = 48000
for path in sys.argv[2:]:
    h = wave.open(path)
    rate = h.getframerate()
    d = numpy.frombuffer(h.readframes(h.getnframes()), dtype="<i2").astype(float)
    level = numpy.abs(d)
    loud = numpy.where(level > level.max() * 0.01)[0]
    if len(loud):
        d = d[max(loud[0] - int(0.05 * rate), 0):loud[-1] + int(0.05 * rate)]
    pieces.append(d)
    pieces.append(numpy.zeros(int(0.12 * rate)))
joined = numpy.concatenate(pieces)
out = wave.open(sys.argv[1], "wb")
out.setnchannels(1)
out.setsampwidth(2)
out.setframerate(rate)
out.writeframes(joined.astype("<i2").tobytes())
