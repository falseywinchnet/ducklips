"""Resonances of an impulse response: peaks of its spectrum, with -3 dB widths.
usage: python3 resonances.py RESPONSE.wav [max_hz]"""
import sys
import wave

import numpy

handle = wave.open(sys.argv[1], "rb")
rate = handle.getframerate()
data = numpy.frombuffer(handle.readframes(handle.getnframes()), dtype="<i2").astype(float)[int(rate * 0.0015):]
limit = float(sys.argv[2]) if len(sys.argv) > 2 else 6000.0
size = 1 << 16
spectrum = numpy.abs(numpy.fft.rfft(data, size))
frequencies = numpy.fft.rfftfreq(size, 1.0 / rate)
level = 20 * numpy.log10(spectrum / numpy.maximum(numpy.fft.rfftfreq(size, 1.0 / rate), 50.0) ** 2 + 1e-30)
found = []
for index in range(2, len(level) - 2):
    if frequencies[index] > limit:
        break
    if frequencies[index] < 100:
        continue
    if level[index] > level[index - 1] and level[index] >= level[index + 1] and level[index] > level[index - 2] and level[index] >= level[index + 2]:
        low = index
        while low > 0 and level[low] > level[index] - 3:
            low -= 1
        high = index
        while high < len(level) - 1 and level[high] > level[index] - 3:
            high += 1
        width = frequencies[high] - frequencies[low]
        found.append((frequencies[index], width, level[index]))
top = max(item[2] for item in found)
print(" ".join("%.0f(%.0f,%+.0f)" % (f, w, l - top) for f, w, l in found if l > top - 40))
