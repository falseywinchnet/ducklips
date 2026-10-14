"""Voice-quality measures for a sustained vowel render.

usage: python3 voice_quality.py VOICE.wav [FLOW.wav]

For the voice: median F0, harmonic-to-noise ratio per band (harmonic energy
against the energy between harmonics, 100 ms frames), and the level of each
band relative to 0-1 kHz. For the glottal flow (optional): open quotient,
closing quotient, and the slope of the flow's harmonic levels (dB/octave).
"""
import struct
import sys

import numpy

BANDS = ((0.0, 1000.0), (1000.0, 2000.0), (2000.0, 4000.0), (4000.0, 8000.0))


def read_wav(path):
    with open(path, "rb") as handle:
        data = handle.read()
    position = 12
    rate = 0
    channels = 1
    bits = 16
    format_code = 1
    samples = None
    while position + 8 <= len(data):
        chunk = data[position:position + 4]
        size = struct.unpack("<I", data[position + 4:position + 8])[0]
        body = data[position + 8:position + 8 + size]
        if chunk == b"fmt ":
            format_code, channels, rate = struct.unpack("<HHI", body[:8])
            bits = struct.unpack("<H", body[14:16])[0]
        elif chunk == b"data":
            if format_code == 3 and bits == 32:
                samples = numpy.frombuffer(body, dtype="<f4").astype(numpy.float64)
            else:
                samples = numpy.frombuffer(body, dtype="<i2").astype(numpy.float64) / 32768.0
        position += 8 + size + (size & 1)
    mono = samples.reshape(-1, channels)[:, 0]
    return mono, float(rate)


def frame_f0(frame, rate):
    centred = frame - numpy.mean(frame)
    spectrum = numpy.fft.rfft(centred, 2 * len(centred))
    correlation = numpy.fft.irfft(numpy.abs(spectrum) ** 2)[: len(centred)]
    low = int(rate / 400.0)
    high = int(rate / 60.0)
    lag = low + int(numpy.argmax(correlation[low:high]))
    if correlation[0] <= 0.0 or correlation[lag] / correlation[0] < 0.3:
        return 0.0
    return rate / lag


def band_measures(signal, rate):
    frame_length = int(0.1 * rate)
    window = numpy.hanning(frame_length)
    harmonic = numpy.zeros(len(BANDS))
    between = numpy.zeros(len(BANDS))
    total = numpy.zeros(len(BANDS))
    pitches = []
    start = int(0.15 * rate)
    stop = len(signal) - int(0.15 * rate) - frame_length
    for begin in range(start, max(stop, start + 1), frame_length // 2):
        frame = signal[begin:begin + frame_length]
        if len(frame) < frame_length:
            break
        f0 = frame_f0(frame, rate)
        if f0 <= 0.0:
            continue
        pitches.append(f0)
        power = numpy.abs(numpy.fft.rfft(frame * window, 8 * frame_length)) ** 2
        frequency = numpy.fft.rfftfreq(8 * frame_length, 1.0 / rate)
        phase = frequency / f0 - numpy.round(frequency / f0)
        near = numpy.abs(phase) < 0.2
        far = numpy.abs(phase) > 0.3
        for index, (low, high) in enumerate(BANDS):
            inside = (frequency >= low) & (frequency < high)
            harmonic[index] += numpy.sum(power[inside & near])
            between[index] += numpy.sum(power[inside & far])
            total[index] += numpy.sum(power[inside])
    return pitches, harmonic, between, total


def flow_measures(flow, rate, f0):
    period = int(round(rate / f0))
    begin = int(0.3 * rate)
    cycles = []
    for cycle in range(10):
        piece = flow[begin + cycle * period: begin + (cycle + 1) * period]
        if len(piece) == period:
            cycles.append(piece)
    open_quotients = []
    closing_quotients = []
    for piece in cycles:
        floor = numpy.min(piece)
        peak = numpy.max(piece)
        threshold = floor + 0.05 * (peak - floor)
        open_quotients.append(numpy.mean(piece > threshold))
        derivative = numpy.diff(piece)
        trough = int(numpy.argmin(derivative))
        top = int(numpy.argmax(piece))
        closing_quotients.append(((trough - top) % period) / period)
    spectrum = numpy.abs(numpy.fft.rfft(flow[begin:begin + 20 * period] * numpy.hanning(20 * period)))
    levels = []
    for harmonic in range(1, 41):
        bin_index = harmonic * 20
        levels.append(20.0 * numpy.log10(numpy.max(spectrum[bin_index - 2:bin_index + 3]) + 1.0e-30))
    octaves = numpy.log2(numpy.arange(1, 41))
    slope = numpy.polyfit(octaves[1:], numpy.array(levels[1:]), 1)[0]
    return numpy.median(open_quotients), numpy.median(closing_quotients), slope, levels[0] - levels[1]


def cycle_measures(flow, rate, f0):
    """Local jitter and shimmer (percent) from the flow's cycle peaks."""
    period = rate / f0
    begin = int(0.2 * rate)
    stop = len(flow) - int(0.2 * rate)
    peaks = []
    position = begin
    while position + int(1.5 * period) < stop:
        window = flow[position:position + int(period)]
        top = position + int(numpy.argmax(window))
        peaks.append(top)
        position = top + int(0.5 * period)
    peaks = numpy.array(peaks, dtype=float)
    refined = []
    for top in peaks.astype(int):
        left, middle, right = flow[top - 1], flow[top], flow[top + 1]
        curvature = left - 2.0 * middle + right
        offset = 0.5 * (left - right) / curvature if curvature != 0.0 else 0.0
        refined.append((top + offset, middle))
    times = numpy.array([item[0] for item in refined])
    heights = numpy.array([item[1] for item in refined])
    periods = numpy.diff(times)
    jitter = 100.0 * numpy.mean(numpy.abs(numpy.diff(periods))) / numpy.mean(periods)
    shimmer = 100.0 * numpy.mean(numpy.abs(numpy.diff(heights))) / numpy.mean(heights)
    return jitter, shimmer


def main():
    voice, rate = read_wav(sys.argv[1])
    pitches, harmonic, between, total = band_measures(voice, rate)
    f0 = float(numpy.median(pitches)) if pitches else 0.0
    print("F0 median %.1f Hz (frames %d)" % (f0, len(pitches)))
    reference = total[0]
    for index, (low, high) in enumerate(BANDS):
        hnr = 10.0 * numpy.log10(harmonic[index] / max(between[index], 1.0e-30))
        level = 10.0 * numpy.log10(total[index] / max(reference, 1.0e-30))
        print("  %5.0f-%5.0f Hz  level %6.1f dB  HNR %5.1f dB" % (low, high, level, hnr))
    if len(sys.argv) > 2 and f0 > 0.0:
        flow, flow_rate = read_wav(sys.argv[2])
        open_quotient, closing_quotient, slope, h1h2 = flow_measures(flow, flow_rate, f0)
        print("flow: open quotient %.2f  closing quotient %.2f  harmonic slope %.1f dB/oct  H1-H2 %.1f dB"
              % (open_quotient, closing_quotient, slope, h1h2))
        jitter, shimmer = cycle_measures(flow, flow_rate, f0)
        print("flow: jitter %.2f %%  shimmer %.2f %%" % (jitter, shimmer))


main()
