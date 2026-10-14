"""Source/filter swap between our rendering and a human reading of the same
ARCTIC sentence, to find which layer makes ours sound wrong.

Each signal is split by LPC (order 18, 25 ms frames every 5 ms, 16 kHz)
into an excitation (the inverse-filtered residual, carrying the glottal
source, its noise and its timing) and a filter (the all-pole envelope:
vocal tract, articulation and radiation). The two are aligned by dynamic
time warping on log-spectral frames. Written, all at 16 kHz:

  OUT_human.wav           the human reading
  OUT_ours.wav            ours
  OUT_oursource.wav       our excitation through the human's filters
  OUT_humansource.wav     the human's excitation through our filters

Each hybrid keeps its excitation's timeline and frame energy: if
"oursource" sounds wrong and "humansource" does not, the source is the
problem; the reverse puts it in the tract and articulation.
usage: cross_synth.py OUT UTTERANCE   (SPEAKER=rms by default)"""
import os
import subprocess
import sys

import numpy
from scipy.linalg import solve_toeplitz

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca

RATE = 16000
ORDER = 18
FRAME = int(0.025 * RATE)
HOP = int(0.005 * RATE)
EMPHASIS = 0.94


def lpc_frames(data):
    """Per-hop LPC polynomials (a[0] = 1) and residual frame energies."""
    window = numpy.hanning(FRAME)
    emphasised = numpy.append(data[0], data[1:] - EMPHASIS * data[:-1])
    padded = numpy.concatenate([numpy.zeros(FRAME // 2), emphasised, numpy.zeros(FRAME)])
    count = len(data) // HOP
    polynomials = numpy.zeros((count, ORDER + 1))
    for index in range(count):
        segment = padded[index * HOP:index * HOP + FRAME] * window
        correlation = numpy.correlate(segment, segment, "full")[FRAME - 1:FRAME + ORDER]
        correlation[0] *= 1.0 + 1e-6
        correlation[0] += 1e-9
        coefficients = solve_toeplitz(correlation[:ORDER], -correlation[1:ORDER + 1])
        polynomials[index] = numpy.concatenate([[1.0], coefficients])
    return polynomials


def log_spectra(data):
    window = numpy.hanning(FRAME)
    padded = numpy.concatenate([numpy.zeros(FRAME // 2), data, numpy.zeros(FRAME)])
    count = len(data) // HOP
    rows = [numpy.log(numpy.abs(numpy.fft.rfft(padded[i * HOP:i * HOP + FRAME] * window, 512))[4:200] + 1e-6)
            for i in range(count)]
    rows = numpy.array(rows)
    return rows - rows.mean(axis=1, keepdims=True)


def warp_map(source, target):
    """For each frame of source, the DTW-matched frame of target."""
    n, m = len(source), len(target)
    cost = numpy.sqrt(((source[:, None, :] - target[None, :, :]) ** 2).sum(axis=2))
    total = numpy.full((n + 1, m + 1), numpy.inf)
    total[0, 0] = 0.0
    for i in range(1, n + 1):
        row = total[i - 1]
        for j in range(1, m + 1):
            total[i, j] = cost[i - 1, j - 1] + min(row[j], total[i, j - 1], row[j - 1])
    mapping = numpy.zeros(n, dtype=int)
    i, j = n, m
    while i > 0 and j > 0:
        mapping[i - 1] = j - 1
        step = numpy.argmin([total[i - 1, j - 1], total[i - 1, j], total[i, j - 1]])
        if step == 0:
            i, j = i - 1, j - 1
        elif step == 1:
            i -= 1
        else:
            j -= 1
    return mapping


def residual(data, polynomials):
    emphasised = numpy.append(data[0], data[1:] - EMPHASIS * data[:-1])
    out = numpy.zeros(len(data))
    for index in range(len(data)):
        frame = min(index // HOP, len(polynomials) - 1)
        a = polynomials[frame]
        history = emphasised[max(index - ORDER, 0):index + 1][::-1]
        out[index] = numpy.dot(a[:len(history)], history)
    return out


def synthesise(excitation, polynomials, mapping):
    """All-pole filtering of excitation with the mapped frame's polynomial."""
    out = numpy.zeros(len(excitation))
    for index in range(len(excitation)):
        frame = min(index // HOP, len(mapping) - 1)
        a = polynomials[mapping[frame]]
        past = out[max(index - ORDER, 0):index][::-1]
        out[index] = excitation[index] - numpy.dot(a[1:1 + len(past)], past)
    # Undo the pre-emphasis.
    restored = numpy.zeros(len(out))
    for index in range(len(out)):
        restored[index] = out[index] + (EMPHASIS * restored[index - 1] if index > 0 else 0.0)
    return restored


def normalised(data):
    return 0.5 * data / (numpy.abs(data).max() + 1e-12)


def main():
    out, name = sys.argv[1], sys.argv[2]
    speaker = os.environ.get("SPEAKER", "rms")
    human, rate = ca.read_wave(os.path.join(ca.ROOT, "ref", "arctic_" + speaker, "wav", name + ".wav"))
    human = ca.resample(human, rate, RATE)
    wave48 = os.path.join("/tmp", name + "_cross48.wav")
    subprocess.run([os.path.join(ca.ROOT, "build", "say"), ca.transcripts()[name], wave48] + sys.argv[3:],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    ours, our_rate = ca.read_wave(wave48)
    ours = ca.resample(ours, our_rate, RATE)
    # Trim both to speech (35 dB under the peak, 50 ms margin).
    trimmed = []
    for data in (human, ours):
        level = numpy.array([numpy.mean(data[i:i + HOP] ** 2) for i in range(0, len(data) - HOP, HOP)])
        keep = numpy.where(10 * numpy.log10(level + 1e-20) > 10 * numpy.log10(level.max()) - 35)[0]
        trimmed.append(data[max(keep[0] - 10, 0) * HOP:(keep[-1] + 10) * HOP])
    human, ours = trimmed
    human_poly = lpc_frames(human)
    our_poly = lpc_frames(ours)
    human_spec = log_spectra(human)
    our_spec = log_spectra(ours)
    ours_to_human = warp_map(our_spec, human_spec)
    human_to_ours = warp_map(human_spec, our_spec)
    our_source = synthesise(residual(ours, our_poly), human_poly, ours_to_human)
    human_source = synthesise(residual(human, human_poly), our_poly, human_to_ours)
    for suffix, data in (("human", human), ("ours", ours), ("oursource", our_source), ("humansource", human_source)):
        ca.write_wave(out + "_" + suffix + ".wav", normalised(data), RATE)
        print("wrote", out + "_" + suffix + ".wav")


if __name__ == "__main__":
    main()
