"""Where the burst sits inside each labelled stop before a vowel: the frame
(2 ms) of steepest rise in energy above 1.5 kHz within the label. Prints, per
stop, the median closure (label start to burst) and post-burst (burst to
label end) durations in ms, for the human labels and for our score.
usage: burst_position.py [FIRST LAST]"""
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca

STOPS = ("p", "t", "k", "b", "d", "g")


def high_energy(data, rate):
    spectrum_step = int(0.002 * rate)
    size = int(0.008 * rate)
    window = numpy.hanning(size)
    hz = numpy.fft.rfftfreq(size, 1.0 / rate)
    band = hz > 1500
    frames = []
    for begin in range(0, len(data) - size, spectrum_step):
        spectrum = numpy.abs(numpy.fft.rfft(data[begin:begin + size] * window)) ** 2
        frames.append(10 * numpy.log10(spectrum[band].sum() + 1e-12))
    return numpy.array(frames), spectrum_step / rate, size / (2.0 * rate)


def survey(side, first, last):
    rows = {}
    for number in range(first, last + 1):
        name = "arctic_a%04d" % number
        if side == "ref":
            data, rate = ca.read_wave(os.path.join(ca.REF, "wav", name + ".wav"))
            phones = ca.reference_phones(name)
        else:
            data, rate = ca.read_wave(os.path.join(ca.OUT, name + "_ours.wav"))
            phones = ca.our_phones(os.path.join(ca.OUT, name + ".score"))
        energy, step, centre = high_energy(data, rate)
        for index in range(len(phones) - 1):
            phone, start, end = phones[index]
            if phone not in STOPS or phones[index + 1][0] not in ca.VOWELS:
                continue
            a = int((start - centre) / step)
            b = int((end + 0.01 - centre) / step)
            if b - a < 6 or a < 1:
                continue
            rise = energy[a + 2:b] - energy[a:b - 2]
            burst = start + (int(numpy.argmax(rise)) + 1) * step
            rows.setdefault(phone, []).append((1000 * (burst - start), 1000 * (end - burst)))
    return rows


if __name__ == "__main__":
    first = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    last = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    human = survey("ref", first, last)
    ours = survey("ours", first, last)
    print("stop   n   human closure/after   ours closure/after (median ms)")
    for phone in STOPS:
        h = numpy.median(numpy.array(human.get(phone, [(0, 0)])), axis=0)
        o = numpy.median(numpy.array(ours.get(phone, [(0, 0)])), axis=0)
        print("%-4s %3d   %6.0f %6.0f        %6.0f %6.0f" % (phone, len(human.get(phone, [])), h[0], h[1], o[0], o[1]))
