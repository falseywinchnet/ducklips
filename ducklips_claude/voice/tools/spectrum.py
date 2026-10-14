"""Average spectra (dB) of segments, overlaid.
usage: python3 spectrum.py OUT.png FILE.wav:start:duration [...]"""
import sys
import wave

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pyplot
import numpy


def read(spec):
    parts = spec.split(":")
    handle = wave.open(parts[0], "rb")
    rate = handle.getframerate()
    data = numpy.frombuffer(handle.readframes(handle.getnframes()), dtype="<i2").astype(float) / 32768.0
    data = data.reshape(-1, handle.getnchannels()).mean(axis=1)
    start = float(parts[1])
    duration = float(parts[2])
    return rate, data[int(start * rate):int((start + duration) * rate)]


figure = pyplot.figure(figsize=(11, 5))
for spec in sys.argv[2:]:
    rate, data = read(spec)
    size = 4096
    window = numpy.hanning(size)
    total = numpy.zeros(size // 2 + 1)
    count = 0
    for start in range(0, len(data) - size, size // 4):
        total += numpy.abs(numpy.fft.rfft(data[start:start + size] * window)) ** 2
        count += 1
    if count == 0:
        padded = numpy.pad(data, (0, size - len(data)))
        total = numpy.abs(numpy.fft.rfft(padded * window)) ** 2
        count = 1
    level = 10 * numpy.log10(total / count + 1e-20)
    frequencies = numpy.fft.rfftfreq(size, 1.0 / rate)
    level -= level[(frequencies > 80) & (frequencies < 6000)].max()
    pyplot.plot(frequencies, level, linewidth=0.7, label=spec.split("/")[-1])
pyplot.xlim(0, 8000)
pyplot.ylim(-80, 5)
pyplot.grid(alpha=0.3)
pyplot.legend(fontsize=8)
pyplot.xlabel("Hz")
pyplot.ylabel("dB")
figure.tight_layout()
figure.savefig(sys.argv[1], dpi=90)
