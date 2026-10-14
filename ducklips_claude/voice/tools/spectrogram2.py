"""Side-by-side wideband spectrograms with waveforms.

usage: python3 spectrogram.py OUT.png FILE.wav[:start:duration] [FILE2.wav[:start:duration] ...]
"""
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
    width = handle.getsampwidth()
    channels = handle.getnchannels()
    frames = handle.readframes(handle.getnframes())
    if width == 2:
        data = numpy.frombuffer(frames, dtype="<i2").astype(numpy.float64) / 32768.0
    else:
        data = numpy.frombuffer(frames, dtype="<i4").astype(numpy.float64) / 2147483648.0
    data = data.reshape(-1, channels).mean(axis=1)
    if len(parts) >= 3:
        start = float(parts[1])
        duration = float(parts[2])
        data = data[int(start * rate):int((start + duration) * rate)]
    return parts[0], rate, data


def spectrogram(data, rate, window_seconds):
    size = int(rate * window_seconds)
    hop = max(1, int(rate * 0.002))
    window = numpy.hanning(size)
    fft_size = 1 << int(numpy.ceil(numpy.log2(size * 2)))
    frames = []
    for start in range(0, max(len(data) - size, 1), hop):
        segment = data[start:start + size]
        if len(segment) < size:
            segment = numpy.pad(segment, (0, size - len(segment)))
        frames.append(numpy.abs(numpy.fft.rfft(segment * window, fft_size)))
    image = 20 * numpy.log10(numpy.array(frames).T + 1e-9)
    frequencies = numpy.fft.rfftfreq(fft_size, 1.0 / rate)
    return image, frequencies, hop / rate


def main():
    output = sys.argv[1]
    specs = sys.argv[2:]
    figure, axes = pyplot.subplots(2, len(specs), figsize=(6 * len(specs), 7), squeeze=False,
                                   gridspec_kw={"height_ratios": [1, 4]})
    for column, spec in enumerate(specs):
        name, rate, data = read(spec)
        times = numpy.arange(len(data)) / rate
        axes[0][column].plot(times, data, linewidth=0.4)
        axes[0][column].set_xlim(0, times[-1])
        axes[0][column].set_title(spec.split("/")[-1], fontsize=9)
        image, frequencies, step = spectrogram(data, rate, float(__import__("os").environ.get("WIN","0.006")))
        top = image.max()
        limit = numpy.searchsorted(frequencies, float(__import__("os").environ.get("TOP","8000")))
        axes[1][column].imshow(image[:limit], origin="lower", aspect="auto", cmap="magma", vmin=top - 70, vmax=top,
                               extent=[0, image.shape[1] * step, 0, frequencies[limit - 1]])
        axes[1][column].set_ylabel("Hz")
        axes[1][column].set_xlabel("s")
    figure.tight_layout()
    figure.savefig(output, dpi=90)


main()
