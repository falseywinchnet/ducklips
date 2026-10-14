"""Finds narrow spectral lines that do not sit on the voice's harmonics.

usage: python3 inharmonic.py FILE.wav [from_s to_s]

A tone that does not follow the pitch (a ringing mode, a numerical
artifact) shows as a narrow peak between harmonics. Each 100 ms frame is
searched for peaks more than 12 dB above the local median spectrum and at
least 0.3 harmonic spacings from any harmonic; peaks that recur across
frames at the same frequency are reported.
"""
import sys

import numpy

exec(open(__file__.replace("inharmonic.py", "voice_compare.py")).read().split("def levinson")[0])


def main():
    signal, rate = read_wav(sys.argv[1])
    begin = float(sys.argv[2]) if len(sys.argv) > 3 else 0.2
    end = float(sys.argv[3]) if len(sys.argv) > 3 else len(signal) / rate - 0.2
    frame_length = int(0.1 * rate)
    fft_length = 4 * frame_length
    window = numpy.blackman(frame_length)
    frequency = numpy.fft.rfftfreq(fft_length, 1.0 / rate)
    counts = {}
    frames = 0
    for start in range(int(begin * rate), int(end * rate) - frame_length, frame_length // 2):
        frame = signal[start:start + frame_length]
        f0, strength = periodicity(frame[: int(0.04 * rate)], rate)
        if strength < 0.5:
            continue
        frames += 1
        level = 10.0 * numpy.log10(numpy.abs(numpy.fft.rfft(frame * window, fft_length)) ** 2 + 1.0e-20)
        kernel = int(300.0 / (rate / fft_length))
        for index in range(kernel, len(level) - kernel):
            if frequency[index] < 150.0:
                continue
            if level[index] < level[index - 1] or level[index] < level[index + 1]:
                continue
            local = numpy.median(level[index - kernel:index + kernel])
            if level[index] - local < 12.0:
                continue
            phase = abs(frequency[index] / f0 - round(frequency[index] / f0))
            if phase < 0.3:
                continue
            key = int(round(frequency[index] / 25.0)) * 25
            counts[key] = counts.get(key, 0) + 1
    recurring = sorted((count, key) for key, count in counts.items() if count >= max(3, frames // 4))
    print("%s: %d voiced frames; recurring inharmonic lines:" % (sys.argv[1], frames))
    for count, key in reversed(recurring):
        print("  %6d Hz in %d frames" % (key, count))


main()
