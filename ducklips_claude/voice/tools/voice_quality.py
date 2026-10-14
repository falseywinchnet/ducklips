"""Whole-voice quality, ours against the human, over ARCTIC vowels of at
least 80 ms (middle 60 ms):
  jitter   mean |period difference| / mean period (%)
  shimmer  mean |cycle peak difference| / mean peak (%)
  H1-H2    first minus second harmonic (dB): open phase, breathiness
  HNR      harmonics-to-noise ratio from the normalized autocorrelation (dB)
usage: voice_quality.py [FIRST LAST]"""
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca


def lowpass(signal, rate, hz):
    size = max(int(rate / hz), 1)
    return numpy.convolve(signal, numpy.ones(size) / size, "same")


def measures(segment, rate):
    segment = segment - segment.mean()
    lags = numpy.arange(int(rate / 250), int(rate / 70))
    energy = numpy.dot(segment, segment) + 1e-12
    correlation = numpy.array([numpy.dot(segment[:-lag], segment[lag:]) for lag in lags]) / energy
    best = int(numpy.argmax(correlation))
    period = lags[best]
    r = min(max(correlation[best] * len(segment) / (len(segment) - period), 1e-6), 0.999999)
    hnr = 10 * numpy.log10(r / (1 - r))
    spectrum = numpy.abs(numpy.fft.rfft(segment * numpy.hanning(len(segment)), 8192))
    hz = numpy.fft.rfftfreq(8192, 1.0 / rate)
    f0 = rate / period

    def harmonic(k):
        band = (hz > k * f0 * 0.85) & (hz < k * f0 * 1.15)
        return 20 * numpy.log10(spectrum[band].max() + 1e-12)

    h1h2 = harmonic(1) - harmonic(2)
    # cycle-by-cycle: positive peaks of the low-passed signal, one per period
    smooth = lowpass(segment, rate, 1000)
    peaks = []
    index = 0
    while index + period <= len(smooth):
        window = smooth[index:index + period]
        peaks.append(index + int(numpy.argmax(window)))
        index += period
    peaks = numpy.array(peaks)
    if len(peaks) < 4:
        return None
    periods = numpy.diff(peaks).astype(float)
    good = numpy.abs(periods - period) < 0.3 * period
    periods = periods[good]
    amplitudes = smooth[peaks]
    if len(periods) < 3:
        return None
    jitter = numpy.mean(numpy.abs(numpy.diff(periods))) / numpy.mean(periods) * 100
    shimmer = numpy.mean(numpy.abs(numpy.diff(amplitudes))) / (numpy.mean(numpy.abs(amplitudes)) + 1e-12) * 100
    return jitter, shimmer, h1h2, hnr


def survey(side, first, last):
    rows = []
    for number in range(first, last + 1):
        name = "arctic_a%04d" % number
        if side == "ref":
            data, rate = ca.read_wave(os.path.join(ca.REF, "wav", name + ".wav"))
            phones = ca.reference_phones(name)
        else:
            data, rate = ca.read_wave(os.path.join(ca.OUT, name + "_ours.wav"))
            phones = ca.our_phones(os.path.join(ca.OUT, name + ".score"))
        for phone, start, end in phones:
            if phone not in ca.VOWELS or end - start < 0.08:
                continue
            middle = 0.5 * (start + end)
            segment = data[int((middle - 0.03) * rate):int((middle + 0.03) * rate)]
            result = measures(segment, rate)
            if result is not None:
                rows.append(result)
    return numpy.array(rows)


if __name__ == "__main__":
    first = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    last = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    print("%-6s %5s  %8s %8s %7s %6s" % ("", "n", "jitter%", "shimmer%", "H1-H2", "HNR"))
    for side in ("ref", "ours"):
        rows = survey(side, first, last)
        median = numpy.median(rows, axis=0)
        print("%-6s %5d  %8.2f %8.2f %7.1f %6.1f" % ("human" if side == "ref" else "ours", len(rows), *median))
