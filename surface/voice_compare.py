"""Compares voiced speech between recordings with identical measures.

usage: python3 voice_compare.py FILE.wav [FILE.wav ...]

Per file, over strongly voiced frames only (40 ms, autocorrelation peak
above 0.6, level within 30 dB of the loudest frame):
- band levels of the long-term spectrum relative to 0-1 kHz, to 22 kHz;
- harmonic-to-noise ratio per band (energy near harmonics against energy
  between them), short frames so pitch motion smears less;
- F0 median, 10th-90th percentile, and median frame-to-frame F0 change;
- formant bandwidths from LPC roots on the voiced frames (median of the
  first three resonances' bandwidths);
- voiced-segment onset rise time (10 % to 90 % of the segment's peak level)
  and the share of segments whose onset has noise before periodicity.
"""
import struct
import sys

import numpy

FRAME_SECONDS = 0.04
HOP_SECONDS = 0.01
BAND_EDGES = (0.0, 1000.0, 2000.0, 4000.0, 6000.0, 8000.0, 12000.0, 16000.0, 22000.0)


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
                usable = len(body) - (len(body) % 2)
                samples = numpy.frombuffer(body[:usable], dtype="<i2").astype(numpy.float64) / 32768.0
        position += 8 + size + (size & 1)
    usable = len(samples) - (len(samples) % channels)
    mono = samples[:usable].reshape(-1, channels)[:, 0]
    return mono, float(rate)


def periodicity(frame, rate):
    centred = frame - numpy.mean(frame)
    spectrum = numpy.fft.rfft(centred, 2 * len(centred))
    correlation = numpy.fft.irfft(numpy.abs(spectrum) ** 2)[: len(centred)]
    low = int(rate / 400.0)
    high = int(rate / 60.0)
    if correlation[0] <= 0.0:
        return 0.0, 0.0
    lag = low + int(numpy.argmax(correlation[low:high]))
    # Normalize by the overlap so long lags are not penalized.
    strength = correlation[lag] / correlation[0] * len(centred) / (len(centred) - lag)
    left = correlation[lag - 1]
    middle = correlation[lag]
    right = correlation[lag + 1]
    curvature = left - 2.0 * middle + right
    offset = 0.5 * (left - right) / curvature if curvature != 0.0 else 0.0
    return rate / (lag + offset), strength


def levinson(autocorrelation, order):
    coefficients = numpy.zeros(order + 1)
    coefficients[0] = 1.0
    error = autocorrelation[0]
    for index in range(1, order + 1):
        reflection = -numpy.dot(coefficients[:index], autocorrelation[index:0:-1]) / error
        previous = coefficients.copy()
        for inner in range(1, index + 1):
            coefficients[inner] = previous[inner] + reflection * previous[index - inner]
        error *= 1.0 - reflection * reflection
    return coefficients


def formant_bandwidths(frame, rate):
    """First three LPC resonances (frequency, bandwidth) at a 10 kHz analysis rate."""
    step = max(1, int(round(rate / 10000.0)))
    smoothed = numpy.convolve(frame, numpy.ones(step) / step, mode="same")[::step]
    analysis_rate = rate / step
    emphasized = numpy.append(smoothed[0], smoothed[1:] - 0.97 * smoothed[:-1])
    windowed = emphasized * numpy.hamming(len(emphasized))
    order = 12
    full = numpy.correlate(windowed, windowed, mode="full")
    autocorrelation = full[len(windowed) - 1: len(windowed) + order]
    if autocorrelation[0] <= 0.0:
        return []
    roots = numpy.roots(levinson(autocorrelation, order))
    resonances = []
    for root in roots:
        if numpy.imag(root) <= 0.0:
            continue
        frequency = numpy.angle(root) * analysis_rate / (2.0 * numpy.pi)
        bandwidth = -numpy.log(numpy.abs(root)) * analysis_rate / numpy.pi
        if 200.0 < frequency < 4500.0 and bandwidth < 600.0:
            resonances.append((frequency, bandwidth))
    resonances.sort()
    return resonances[:3]


def analyse(path):
    signal, rate = read_wav(path)
    frame_length = int(FRAME_SECONDS * rate)
    hop = int(HOP_SECONDS * rate)
    starts = range(0, len(signal) - frame_length, hop)
    levels = numpy.array([numpy.sqrt(numpy.mean(signal[s:s + frame_length] ** 2)) + 1.0e-12 for s in starts])
    loudest = numpy.max(levels)
    window = numpy.hanning(frame_length)
    fft_length = 4 * frame_length
    frequency = numpy.fft.rfftfreq(fft_length, 1.0 / rate)
    band_count = len(BAND_EDGES) - 1
    long_term = numpy.zeros(band_count)
    near_total = numpy.zeros(band_count)
    far_total = numpy.zeros(band_count)
    pitches = []
    voiced = numpy.zeros(len(levels), dtype=bool)
    frame_pitch = numpy.zeros(len(levels))
    bandwidths = [[], [], []]
    for index, start in enumerate(starts):
        if levels[index] < loudest * 10.0 ** (-30.0 / 20.0):
            continue
        frame = signal[start:start + frame_length]
        f0, strength = periodicity(frame, rate)
        if strength < 0.6:
            continue
        voiced[index] = True
        frame_pitch[index] = f0
        pitches.append(f0)
        power = numpy.abs(numpy.fft.rfft(frame * window, fft_length)) ** 2
        phase = frequency / f0 - numpy.round(frequency / f0)
        near = numpy.abs(phase) < 0.2
        far = numpy.abs(phase) > 0.3
        for band in range(band_count):
            inside = (frequency >= BAND_EDGES[band]) & (frequency < BAND_EDGES[band + 1])
            long_term[band] += numpy.sum(power[inside])
            near_total[band] += numpy.sum(power[inside & near])
            far_total[band] += numpy.sum(power[inside & far])
        if index % 3 == 0:
            for slot, (centre, bandwidth) in enumerate(formant_bandwidths(frame, rate)):
                bandwidths[slot].append(bandwidth)
    # Frame-to-frame pitch motion within voiced runs (semitones per 10 ms).
    motion = []
    for index in range(1, len(voiced)):
        if voiced[index] and voiced[index - 1]:
            motion.append(abs(12.0 * numpy.log2(frame_pitch[index] / frame_pitch[index - 1])))
    # Onsets: rise time of each voiced run's level and noise before periodicity.
    rises = []
    breathy = 0
    runs = 0
    index = 0
    while index < len(voiced):
        if not voiced[index]:
            index += 1
            continue
        end = index
        while end < len(voiced) and voiced[end]:
            end += 1
        if end - index >= 10:
            runs += 1
            peak = numpy.max(levels[index:end])
            look = max(0, index - 30)
            segment = levels[look:end]
            above_ten = look + int(numpy.argmax(segment > 0.1 * peak))
            above_ninety = look + int(numpy.argmax(segment > 0.9 * peak))
            rises.append((above_ninety - above_ten) * HOP_SECONDS * 1000.0)
            if index - above_ten >= 3:
                breathy += 1
        index = end
    reference = long_term[0]
    print("== %s (%.1f s, %.0f Hz, voiced frames %d)" % (path, len(signal) / rate, rate, len(pitches)))
    for band in range(band_count):
        if BAND_EDGES[band + 1] > rate / 2.0:
            break
        level = 10.0 * numpy.log10(long_term[band] / reference)
        hnr = 10.0 * numpy.log10(near_total[band] / max(far_total[band], 1.0e-30))
        print("  %5.0f-%5.0f Hz  level %6.1f dB  HNR %5.1f dB"
              % (BAND_EDGES[band], BAND_EDGES[band + 1], level, hnr))
    pitch_array = numpy.sort(numpy.array(pitches))
    print("  F0 median %.1f Hz (10th %.1f, 90th %.1f); median motion %.3f semitone/10 ms"
          % (numpy.median(pitch_array), pitch_array[len(pitch_array) // 10],
             pitch_array[(len(pitch_array) * 9) // 10], numpy.median(motion) if motion else 0.0))
    print("  LPC bandwidths (median): B1 %.0f  B2 %.0f  B3 %.0f Hz"
          % tuple(numpy.median(values) if values else 0.0 for values in bandwidths))
    if rises:
        print("  voiced runs %d: onset rise median %.0f ms; noise-led onsets %d"
              % (runs, numpy.median(rises), breathy))


for argument in sys.argv[1:]:
    analyse(argument)
