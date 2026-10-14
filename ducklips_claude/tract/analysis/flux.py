"""How fast the spectrum moves: 20 log-spaced bands (100 Hz-7 kHz), 10 ms
hop, log band energies. Prints, over frames within 40 dB of each file's
peak: the median and 90th-percentile spectral flux (RMS change across bands,
dB per 10 ms), the share of the bands' log-envelope modulation energy in
2-4, 4-8, 8-16 and 16-32 Hz, and the median band-to-band roughness (RMS
difference between neighbouring bands in a frame, dB).
usage: flux.py LABEL GLOB [LABEL GLOB ...]"""
import sys, glob, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly

RATE = 16000
HOP = 160
EDGES = numpy.geomspace(100, 7000, 21)
MODS = [(2, 4), (4, 8), (8, 16), (16, 32)]


def bands(path):
    r, x = wavfile.read(path)
    x = x.astype(float)
    if x.ndim > 1:
        x = x[:, 0]
    if r != RATE:
        x = resample_poly(x, RATE, r)
    n = 400
    window = numpy.hanning(n)
    freqs = numpy.fft.rfftfreq(n, 1 / RATE)
    rows = []
    for i in range(0, len(x) - n, HOP):
        s = numpy.abs(numpy.fft.rfft(x[i:i + n] * window)) ** 2
        rows.append([s[(freqs >= a) & (freqs < b)].sum() for a, b in zip(EDGES[:-1], EDGES[1:])])
    e = 10 * numpy.log10(numpy.array(rows) + 1e-9)
    return e


def measure(paths):
    flux, rough, mods = [], [], numpy.zeros(len(MODS))
    for p in paths:
        e = bands(p)
        total = 10 * numpy.log10(numpy.sum(10 ** (e / 10), axis=1))
        live = total > total.max() - 40
        e = numpy.maximum(e, total.max() - 60)
        d = numpy.sqrt(numpy.mean(numpy.diff(e, axis=0) ** 2, axis=1))
        flux.extend(d[live[1:] & live[:-1]])
        rough.extend(numpy.sqrt(numpy.mean(numpy.diff(e[live], axis=1) ** 2, axis=1)))
        env = e - e.mean(axis=0)
        spec = numpy.abs(numpy.fft.rfft(env * numpy.hanning(len(env))[:, None], axis=0)) ** 2
        mf = numpy.fft.rfftfreq(len(env), HOP / RATE)
        whole = spec[(mf >= 1) & (mf < 32)].sum()
        mods += numpy.array([spec[(mf >= a) & (mf < b)].sum() / whole for a, b in MODS])
    return numpy.median(flux), numpy.percentile(flux, 90), mods / len(paths), numpy.median(rough)


for label, pattern in zip(sys.argv[1::2], sys.argv[2::2]):
    paths = sorted(glob.glob(pattern))
    med, p90, mods, rough = measure(paths)
    print("%-8s files %3d  flux med %.2f p90 %.2f  mod %s  rough %.2f" % (
        label, len(paths), med, p90, " ".join("%d-%d:%.2f" % (a, b, m) for (a, b), m in zip(MODS, mods)), rough))
