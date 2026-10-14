"""Articulatory inversion by codebook lookup.

Every 10 ms frame of speech becomes a spectral envelope (LPC, order 18 at
16 kHz) sampled on the codebook's 1/24-octave grid. Envelope and codebook
responses both lose their best straight line in log-frequency (source tilt
and recording colour are not mouth shape), then each frame takes the
nearest codebook shape by one matrix product over all frames.

usage:
  invert.py CODEBOOK_PREFIX selftest              # recover our own vowels' shapes
  invert.py CODEBOOK_PREFIX arctic OUT.json       # human phone targets from ARCTIC bdl
"""
import json
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
from analyze import lpc

GRID = 100.0 * 2 ** (numpy.arange(0, int(numpy.log2(7600 / 100) * 24) + 1) / 24.0)
LOGF = numpy.log2(GRID)


def detrend(rows):
    design = numpy.stack([numpy.ones_like(LOGF), LOGF - LOGF.mean()], axis=1)
    coefficients, _, _, _ = numpy.linalg.lstsq(design, rows.T, rcond=None)
    return rows - (design @ coefficients).T


def envelope(segment, rate):
    windowed = numpy.append(segment[0], segment[1:] - 0.5 * segment[:-1]) * numpy.hamming(len(segment))
    a = lpc(windowed, 18)
    frequencies = numpy.exp(-2j * numpy.pi * numpy.outer(GRID / rate, numpy.arange(len(a))))
    return -20 * numpy.log10(numpy.abs(frequencies @ a) + 1e-9)


def frames(data, rate, start, end):
    out = []
    for t in numpy.arange(start, end, 0.01):
        seg = data[int((t - 0.0125) * rate):int((t + 0.0125) * rate)]
        if len(seg) == int(0.025 * rate) and numpy.sqrt(numpy.mean(seg ** 2)) > 1e-4:
            out.append(envelope(seg, rate))
    return out


def lpc_formants(segment, rate):
    """F1-F3 from LPC roots (bandwidth < 500 Hz), or None."""
    windowed = numpy.append(segment[0], segment[1:] - 0.94 * segment[:-1]) * numpy.hamming(len(segment))
    found = []
    for root in numpy.roots(lpc(windowed, 14)):
        if root.imag <= 0:
            continue
        f = numpy.angle(root) * rate / (2 * numpy.pi)
        b = -numpy.log(abs(root)) * rate / numpy.pi
        if 180 < f < 3900 and b < 500:
            found.append(f)
    found.sort()
    return found[:3] if len(found) >= 3 else None


def formant_frames(data, rate, start, end):
    out = []
    for t in numpy.arange(start, end, 0.01):
        seg = data[int((t - 0.0125) * rate):int((t + 0.0125) * rate)]
        if len(seg) == int(0.025 * rate) and numpy.sqrt(numpy.mean(seg ** 2)) > 1e-4:
            f = lpc_formants(seg, rate)
            if f:
                out.append(f)
    return out


def row_formants(row):
    """First three resonance peaks of a codebook row (dB on GRID), refined by
    parabolic interpolation in log frequency."""
    peaks = []
    for i in range(1, len(row) - 1):
        if row[i] > row[i - 1] and row[i] >= row[i + 1] and GRID[i] > 180:
            denominator = row[i - 1] - 2 * row[i] + row[i + 1]
            shift = 0.5 * (row[i - 1] - row[i + 1]) / denominator if denominator != 0 else 0.0
            peaks.append(2 ** (LOGF[i] + shift / 24.0))
            if len(peaks) == 3:
                break
    return peaks if len(peaks) == 3 else None


class Codebook:
    def __init__(self, prefix):
        if os.path.exists(prefix + "_params.npy"):
            params = [numpy.load(prefix + "_params.npy")]
            rows = [numpy.load(prefix + "_response.npy")]
        else:
            shard = 0
            params, rows = [], []
            while os.path.exists("%s_%d_params.npy" % (prefix, shard)):
                params.append(numpy.load("%s_%d_params.npy" % (prefix, shard)))
                rows.append(numpy.load("%s_%d_response.npy" % (prefix, shard)))
                shard += 1
        params = numpy.concatenate(params)
        raw = numpy.concatenate(rows).astype(float)
        keep, formants = [], []
        for index, row in enumerate(raw):
            f = row_formants(row)
            if f:
                keep.append(index)
                formants.append(f)
        self.params = params[keep]
        self.log_formants = numpy.log(numpy.array(formants))

    def nearest(self, formant_list):
        """Nearest shape by weighted log-formant distance (F1, F2 weigh most)."""
        weights = numpy.array([1.0, 1.0, 0.5])
        target = numpy.log(numpy.array(formant_list))
        out = []
        for t in target:
            distance = numpy.sum(weights * (self.log_formants - t) ** 2, axis=1)
            out.append(self.params[numpy.argmin(distance)])
        return numpy.array(out)


def selftest(book):
    """Render sustained vowels at known postures and invert them."""
    import subprocess
    cases = [(13.0, 0.29, 9.1, 0.0), (4.0, 1.69, 10.0, 0.0), (7.75, 0.58, 0.71, 0.8), (10.75, 1.19, 6.3, 0.0)]
    for place, area, lips, protrusion in cases:
        path = "/tmp/invert_selftest.wav"
        subprocess.run([os.path.join(ca.ROOT, "build", "vowel_probe"), path, str(place), str(area), str(lips), "110", "750",
                        "0.0", "1.0"], stdout=subprocess.DEVNULL, check=True)
        data, rate = ca.read_wave(path)
        data = ca.resample(data, rate, 16000)
        found = book.nearest(formant_frames(data, 16000, 0.4, 0.9))
        median = numpy.median(found, axis=0)
        print("true place %.2f area %.2f lips %.2f prot %.1f  -> inverted %.2f %.2f %.2f %.1f jaw %.1f" %
              (place, area, lips, protrusion, *median))


def arctic(book, output):
    per_phone = {}
    for number in range(1, 61):
        name = "arctic_a%04d" % number
        data, rate = ca.read_wave(os.path.join(ca.REF, "wav", name + ".wav"))
        for phone, start, end in ca.reference_phones(name):
            if end - start < 0.04:
                continue
            middle = formant_frames(data, rate, start + 0.2 * (end - start), end - 0.2 * (end - start))
            if middle:
                per_phone.setdefault(phone, []).extend(book.nearest(middle).tolist())
    result = {}
    for phone, values in sorted(per_phone.items()):
        values = numpy.array(values)
        median = numpy.median(values, axis=0)
        result[phone] = dict(n=len(values), place=float(median[0]), area=float(median[1]), lips=float(median[2]),
                             protrusion=float(median[3]), jaw=float(median[4]))
        print("%-3s n=%4d  place %5.2f  area %5.2f  lips %5.2f  protrusion %.1f  jaw %.1f" %
              (phone, len(values), *median))
    json.dump(result, open(output, "w"), indent=1)


if __name__ == "__main__":
    book = Codebook(sys.argv[1])
    if sys.argv[2] == "selftest":
        selftest(book)
    else:
        arctic(book, sys.argv[3])
