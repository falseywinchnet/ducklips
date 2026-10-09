"""Writes a spectrogram PNG of one or more WAV files, stacked.

usage: python3 spectrogram.py OUT.png FILE.wav [FILE.wav ...] [--max-hz 8000] [--range 70]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy

exec(open(__file__.replace("spectrogram.py", "voice_compare.py")).read().split("def periodicity")[0])


def main():
    arguments = sys.argv[1:]
    maximum = 8000.0
    dynamic = 70.0
    files = []
    index = 1
    while index < len(arguments):
        if arguments[index] == "--max-hz":
            maximum = float(arguments[index + 1])
            index += 2
        elif arguments[index] == "--range":
            dynamic = float(arguments[index + 1])
            index += 2
        else:
            files.append(arguments[index])
            index += 1
    figure, axes = plt.subplots(len(files), 1, figsize=(12, 3.2 * len(files)), squeeze=False)
    for row, path in enumerate(files):
        signal, rate = read_wav(path)
        frame = int(0.025 * rate)
        hop = int(0.004 * rate)
        window = numpy.hanning(frame)
        columns = []
        for start in range(0, len(signal) - frame, hop):
            spectrum = numpy.abs(numpy.fft.rfft(signal[start:start + frame] * window, 2048)) ** 2
            columns.append(10.0 * numpy.log10(spectrum + 1.0e-14))
        image = numpy.array(columns).T
        frequency = numpy.fft.rfftfreq(2048, 1.0 / rate)
        keep = frequency <= maximum
        image = image[keep]
        top = numpy.max(image)
        axis = axes[row][0]
        axis.imshow(image, origin="lower", aspect="auto", cmap="magma", vmin=top - dynamic, vmax=top,
                    extent=[0.0, len(signal) / rate, 0.0, frequency[keep][-1]])
        axis.set_title(path, fontsize=9)
        axis.set_ylabel("Hz")
    axes[-1][0].set_xlabel("s")
    figure.tight_layout()
    figure.savefig(arguments[0], dpi=90)


main()
