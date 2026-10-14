"""Cycle-level transforms of voiced speech, built to see what the ear hears
rather than long-term averages.

For each file (resampled to 16 kHz): voiced regions from a normalized
autocorrelation F0 track (10 ms hop); glottal cycles from the LP residual
(order 18), one strongest negative residual peak per period. Then, per
auditory band (2nd-order Butterworth band-passes around the centres):

  periodicity  normalized autocorrelation of the band's Hilbert envelope at
               the local period (how strongly the band pulses with the voice)
  repeat       correlation of each cycle's band waveform with the next
               cycle's (both resampled to the same length)
  crest        peak / RMS of the 3-7 kHz band within each cycle

and the cycle-domain spectrum: every cycle warped to 64 samples and the
cycles of each voiced run concatenated; the FFT then holds harmonics at
exact bins. For harmonic k: energy within +-0.02 harmonic of k (line),
within 0.05-0.45 of k (sidebands: modulation slower than the pitch, a few
to tens of Hz), reported as sideband/line in dB for harmonics 1-3, 4-8 and
9-20.

usage: hear.py LABEL WAV[:start:end] [LABEL WAV[:start:end] ...]"""
import sys, glob, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, butter, sosfiltfilt, lfilter, hilbert
from scipy.linalg import solve_toeplitz

RATE = 16000
CENTRES = [300, 700, 1400, 2400, 3600, 5200]


def load(spec):
    parts = spec.split(":")
    paths = sorted(glob.glob(parts[0]))
    out = []
    for p in paths:
        r, x = wavfile.read(p); x = x.astype(float)
        if x.ndim > 1: x = x[:, 0]
        x = resample_poly(x, RATE, r)
        if len(parts) == 3:
            x = x[int(float(parts[1]) * RATE):int(float(parts[2]) * RATE)]
        out.append(x / (numpy.abs(x).max() + 1e-9))
    return out


def f0_track(x):
    hop, win = 160, 640
    f0 = []
    for i in range(0, len(x) - win - 300, hop):
        a = x[i:i + win] - x[i:i + win].mean()
        c = numpy.correlate(a, a, "full")[win - 1:]
        if c[0] <= 0:
            f0.append(0.0); continue
        lo, hi = 53, 230
        k = lo + int(numpy.argmax(c[lo:hi]))
        f0.append(RATE / k if c[k] > 0.55 * c[0] else 0.0)
    return numpy.array(f0)


def residual(x):
    out = numpy.zeros_like(x)
    win, hop = 400, 160
    for i in range(0, len(x) - win, hop):
        seg = x[i:i + win] * numpy.hamming(win)
        c = numpy.correlate(seg, seg, "full")[win - 1:win + 19]
        if c[0] <= 1e-9:
            continue
        a = solve_toeplitz(c[:18], -c[1:19])
        r = lfilter(numpy.r_[1, a], [1], x[i:i + win])
        out[i + 120:i + 280] = r[120:280]
    return out


def cycles(x, f0):
    """Glottal cycle starts (samples) per voiced run: lists of indices."""
    res = residual(x)
    runs, i = [], 0
    while i < len(f0):
        if f0[i] <= 0:
            i += 1; continue
        j = i
        while j < len(f0) and f0[j] > 0: j += 1
        if j - i >= 8:  # 80 ms
            start, end = i * 160 + 320, j * 160 + 320
            marks = []
            t = start
            while t < end:
                period = int(RATE / f0[min(len(f0) - 1, max(0, (t - 320) // 160))])
                if not marks:
                    w = res[t:t + period]
                    if len(w) == 0: break
                    marks.append(t + int(numpy.argmin(w)))
                else:
                    lo, hi = marks[-1] + int(0.7 * period), marks[-1] + int(1.3 * period)
                    if hi >= end: break
                    marks.append(lo + int(numpy.argmin(res[lo:hi])))
                t = marks[-1] + int(0.5 * period) if len(marks) > 1 else marks[-1] + 1
                if len(marks) > 1 and marks[-1] - marks[-2] < 0.6 * period:
                    marks.pop(); break
            if len(marks) >= 8:
                runs.append(marks)
        i = j
    return runs


def band(x, centre):
    lo, hi = centre / 1.35, min(centre * 1.35, 7600)
    return sosfiltfilt(butter(2, [lo, hi], btype="band", fs=RATE, output="sos"), x)


def warp(seg, n):
    return numpy.interp(numpy.linspace(0, len(seg) - 1, n), numpy.arange(len(seg)), seg)


def analyse(xs):
    per = {c: [] for c in CENTRES}; rep = {c: [] for c in CENTRES}; crest = []
    side = {"1-3": [], "4-8": [], "9-20": []}
    for x in xs:
        f0 = f0_track(x)
        runs = cycles(x, f0)
        if not runs:
            continue
        bands = {c: band(x, c) for c in CENTRES}
        envs = {c: numpy.abs(hilbert(bands[c])) for c in CENTRES}
        hiband = sosfiltfilt(butter(4, [3000, 7000], btype="band", fs=RATE, output="sos"), x)
        for marks in runs:
            marks = numpy.array(marks)
            periods = numpy.diff(marks)
            for c in CENTRES:
                e = envs[c]; w = bands[c]
                for a, b, p in zip(marks[:-2], marks[1:-1], periods[1:]):
                    s1, s2 = e[a:b], e[b:b + (b - a)]
                    if len(s2) == len(s1) and s1.std() > 0 and s2.std() > 0:
                        per[c].append(numpy.corrcoef(s1, s2)[0, 1])
                    w1, w2 = warp(w[a:b], 64), warp(w[b:b + p], 64)
                    if w1.std() > 0 and w2.std() > 0:
                        rep[c].append(numpy.corrcoef(w1, w2)[0, 1])
            for a, b in zip(marks[:-1], marks[1:]):
                h = hiband[a:b]
                if len(h) > 8:
                    crest.append(20 * numpy.log10(numpy.abs(h).max() / (h.std() + 1e-12)))
            # cycle-domain spectrum
            seq = numpy.concatenate([warp(x[a:b], 64) for a, b in zip(marks[:-1], marks[1:])])
            n = len(seq)
            spec = numpy.abs(numpy.fft.rfft((seq - seq.mean()) * numpy.hanning(n))) ** 2
            k = numpy.arange(len(spec)) * 64.0 / n  # in harmonics
            for name, (h0, h1) in (("1-3", (1, 3)), ("4-8", (4, 8)), ("9-20", (9, 20))):
                line = sb = 0.0
                for h in range(h0, h1 + 1):
                    d = numpy.abs(k - h)
                    line += spec[d <= 0.02].sum(); sb += spec[(d >= 0.05) & (d <= 0.45)].sum()
                if line > 0:
                    side[name].append(10 * numpy.log10(sb / line + 1e-12))
    return per, rep, crest, side


for label, spec in zip(sys.argv[1::2], sys.argv[2::2]):
    per, rep, crest, side = analyse(load(spec))
    print("%-10s periodicity %s | repeat %s | crest3-7k %.1f dB | sideband/line h1-3 %.1f h4-8 %.1f h9-20 %.1f dB" % (
        label, " ".join("%.2f" % numpy.median(per[c]) if per[c] else "nan" for c in CENTRES),
        " ".join("%.2f" % numpy.median(rep[c]) if rep[c] else "nan" for c in CENTRES),
        numpy.median(crest) if crest else numpy.nan,
        numpy.median(side["1-3"]) if side["1-3"] else numpy.nan, numpy.median(side["4-8"]) if side["4-8"] else numpy.nan,
        numpy.median(side["9-20"]) if side["9-20"] else numpy.nan))
print("bands (Hz):", CENTRES)
