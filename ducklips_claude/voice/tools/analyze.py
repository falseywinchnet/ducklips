"""Quick acoustic summary of a rendered or recorded vowel segment.

usage: python3 analyze.py FILE.wav [start_s] [duration_s]
Prints F0 (autocorrelation), LPC formants and bandwidths, spectral tilt,
harmonic-to-noise ratio estimate, and band energies.
"""
import sys
import wave

import numpy


def read(path):
    handle = wave.open(path, "rb")
    rate = handle.getframerate()
    width = handle.getsampwidth()
    channels = handle.getnchannels()
    frames = handle.readframes(handle.getnframes())
    if width == 2:
        data = numpy.frombuffer(frames, dtype="<i2").astype(numpy.float64) / 32768.0
    elif width == 4:
        data = numpy.frombuffer(frames, dtype="<i4").astype(numpy.float64) / 2147483648.0
    else:
        raise SystemExit("unsupported width")
    data = data.reshape(-1, channels).mean(axis=1)
    return rate, data


def pitch(segment, rate):
    segment = segment - segment.mean()
    correlation = numpy.correlate(segment, segment, "full")[len(segment) - 1:]
    low = int(rate / 400)
    high = int(rate / 60)
    lag = low + int(numpy.argmax(correlation[low:high]))
    return rate / lag, correlation[lag] / correlation[0]


def lpc(segment, order):
    correlation = numpy.correlate(segment, segment, "full")[len(segment) - 1:len(segment) + order]
    coefficients = numpy.zeros(order + 1)
    coefficients[0] = 1.0
    error = correlation[0]
    for index in range(1, order + 1):
        reflection = -(correlation[index] + numpy.dot(coefficients[1:index], correlation[index - 1:0:-1])) / error
        previous = coefficients.copy()
        for inner in range(1, index):
            coefficients[inner] = previous[inner] + reflection * previous[index - inner]
        coefficients[index] = reflection
        error *= 1.0 - reflection * reflection
    return coefficients


def formants(segment, rate):
    target = 11025
    step = max(1, int(round(rate / target)))
    smoothed = numpy.convolve(segment, numpy.ones(step) / step, "same")[::step]
    used = rate / step
    smoothed = numpy.append(smoothed[0], smoothed[1:] - 0.94 * smoothed[:-1])
    smoothed = smoothed * numpy.hamming(len(smoothed))
    roots = numpy.roots(lpc(smoothed, 12))
    found = []
    for root in roots:
        if root.imag <= 0:
            continue
        frequency = numpy.angle(root) * used / (2 * numpy.pi)
        bandwidth = -numpy.log(abs(root)) * used / numpy.pi
        if 150 < frequency < used / 2 - 200 and bandwidth < 600:
            found.append((frequency, bandwidth))
    found.sort()
    return found


def main():
    rate, data = read(sys.argv[1])
    start = float(sys.argv[2]) if len(sys.argv) > 2 else len(data) / rate * 0.5 - 0.1
    duration = float(sys.argv[3]) if len(sys.argv) > 3 else 0.2
    segment = data[int(start * rate):int((start + duration) * rate)]
    f0, clarity = pitch(segment, rate)
    spectrum = numpy.abs(numpy.fft.rfft(segment * numpy.hanning(len(segment))))
    frequencies = numpy.fft.rfftfreq(len(segment), 1.0 / rate)
    bands = [(0, 500), (500, 1000), (1000, 2000), (2000, 4000), (4000, 8000), (8000, 12000)]
    total = numpy.sum(spectrum ** 2) + 1e-30
    text = []
    for low, high in bands:
        mask = (frequencies >= low) & (frequencies < high)
        text.append("%d-%d:%.0f" % (low, high, 10 * numpy.log10(numpy.sum(spectrum[mask] ** 2) / total + 1e-12)))
    harmonics = []
    for number in range(1, 9):
        index = int(round(number * f0 * len(segment) / rate))
        window = spectrum[max(index - 2, 0):index + 3]
        harmonics.append(20 * numpy.log10(window.max() + 1e-12))
    relative = [value - harmonics[0] for value in harmonics]
    print("%s @%.2f: F0 %.1f Hz clarity %.2f rms %.4f" % (sys.argv[1], start, f0, clarity, numpy.sqrt(numpy.mean(segment ** 2))))
    print("  formants: " + ", ".join("%.0f(%.0f)" % pair for pair in formants(segment, rate)))
    print("  H1..H8 rel dB: " + " ".join("%+.0f" % value for value in relative))
    print("  band dB: " + " ".join(text))


if __name__ == "__main__":
    main()
