"""Voice source and output balance of a vowel_probe TRACE file: harmonic
levels of the glottal flow derivative and of the output, averaged per band,
in dB re the strongest harmonic. usage: source_slope.py TRACE F0"""
import sys

import numpy

RATE = 48000
trace = numpy.loadtxt(sys.argv[1])
f0 = float(sys.argv[2])
segment = slice(int(0.6 * RATE), int(0.95 * RATE))


def harmonics(signal):
    x = signal[segment] - signal[segment].mean()
    spectrum = numpy.abs(numpy.fft.rfft(x * numpy.hanning(len(x))))
    hz = numpy.fft.rfftfreq(len(x), 1.0 / RATE)
    # refine f0 on the strongest low harmonic region
    levels = []
    for k in range(1, int(7000 / f0)):
        band = (hz > k * f0 - 0.4 * f0) & (hz < k * f0 + 0.4 * f0)
        levels.append(20 * numpy.log10(spectrum[band].max() + 1e-12))
    levels = numpy.array(levels)
    return numpy.arange(1, len(levels) + 1) * f0, levels - levels.max()


flow = trace[:, 1]
hz, source = harmonics(numpy.diff(flow, prepend=flow[0]))
_, output = harmonics(trace[:, 3])
for low, high in [(0, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, 4000), (4000, 5000), (5000, 7000)]:
    band = (hz >= low) & (hz < high)
    print("%5d-%5d  source %6.1f  output %6.1f" % (low, high, source[band].mean(), output[band].mean()))
