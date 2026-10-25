"""The owner's phase-unrolled spectrogram (2026-10-18): "instead of specgram
which is abs of fft, take the irfft, then slice the top half. you wind up
with a phase unrolled representation."

Each Hann-windowed frame (N samples) goes through numpy.fft.irfft as if it
were a half spectrum: the output (2(N - 1) real values) is mirror-symmetric,
and its top half is the frame's cosine transform, frequency running from the
Nyquist rate down to 0 (here flipped to run upward). Each value is a
frequency's magnitude times the cosine of its phase: signed, real, the phase
kept as sign and swing instead of thrown away.

Shown with a signed log scale, sign(v) log10(1 + |v| / floor), floor 60 dB
below the largest value.
usage: real_spectrogram.py OUT.png TOP_HZ LABEL:WAV[:START:END] ..."""
import sys

import numpy
import soundfile


def real_spectrogram(x, r, frame=1024, hop=128):
    window = numpy.hanning(frame)
    count = 1 + (len(x) - frame) // hop
    rows = []
    for index in range(count):
        piece = x[index * hop:index * hop + frame] * window
        unrolled = numpy.fft.irfft(piece)       # 2 (frame - 1) values, mirror-symmetric
        top = unrolled[frame - 1:]               # Nyquist ... 0
        rows.append(top[::-1])                   # 0 ... Nyquist
    matrix = numpy.array(rows).T                 # frequency x time
    hz = numpy.linspace(0.0, r / 2.0, matrix.shape[0])
    seconds = (numpy.arange(count) * hop + frame / 2) / r
    return matrix, hz, seconds


def signed_log(matrix):
    floor = numpy.abs(matrix).max() * 1e-3
    return numpy.sign(matrix) * numpy.log10(1.0 + numpy.abs(matrix) / floor)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = sys.argv[1]
    top_hz = float(sys.argv[2])
    specs = sys.argv[3:]
    fig, axes = plt.subplots(len(specs), 1, figsize=(14, 3.2 * len(specs)))
    axes = numpy.atleast_1d(axes)
    for ax, spec in zip(axes, specs):
        label, path, *span = spec.split(":")
        x, r = soundfile.read(path)
        x = x if x.ndim == 1 else x.mean(axis=1)
        if span:
            x = x[int(float(span[0]) * r):int(float(span[1]) * r)]
        matrix, hz, seconds = real_spectrogram(x, r)
        keep = hz <= top_hz
        shown = signed_log(matrix[keep])
        limit = numpy.abs(shown).max()
        ax.imshow(shown, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit,
                  extent=[seconds[0], seconds[-1], 0, hz[keep][-1]])
        ax.set_title(label, fontsize=10)
        ax.set_ylabel("Hz")
    axes[-1].set_xlabel("s")
    plt.tight_layout()
    plt.savefig(out, dpi=80)


if __name__ == "__main__":
    main()
