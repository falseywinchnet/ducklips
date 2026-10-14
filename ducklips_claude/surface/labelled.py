"""Wideband spectrogram with phone labels from a surface_say timeline line.

usage: python3 labelled.py OUT.png FILE.wav TIMELINE.txt [from_s to_s]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy

exec(open(__file__.replace("labelled.py", "voice_compare.py")).read().split("def periodicity")[0])


def main():
    signal, rate = read_wav(sys.argv[2])
    labels = []
    for item in open(sys.argv[3]).read().split()[1:]:
        name, time = item.split("@")
        labels.append((float(time), name))
    begin = float(sys.argv[4]) if len(sys.argv) > 5 else 0.0
    end = float(sys.argv[5]) if len(sys.argv) > 5 else len(signal) / rate
    piece = signal[int(begin * rate):int(end * rate)]
    frame = int(0.006 * rate)
    hop = int(0.001 * rate)
    window = numpy.hanning(frame)
    columns = [10.0 * numpy.log10(numpy.abs(numpy.fft.rfft(piece[s:s + frame] * window, 2048)) ** 2 + 1.0e-14)
               for s in range(0, len(piece) - frame, hop)]
    image = numpy.array(columns).T
    frequency = numpy.fft.rfftfreq(2048, 1.0 / rate)
    keep = frequency < 8000.0
    image = image[keep]
    figure, axis = plt.subplots(1, 1, figsize=(16, 4.5))
    axis.imshow(image, origin="lower", aspect="auto", cmap="magma", vmin=image.max() - 65, vmax=image.max(),
                extent=[begin, end, 0.0, 8000.0])
    for time, name in labels:
        if begin <= time <= end:
            axis.axvline(time, color="c", lw=0.6)
            axis.text(time + 0.003, 7400.0, name, color="c", fontsize=9)
    figure.tight_layout()
    figure.savefig(sys.argv[1], dpi=80)


main()
