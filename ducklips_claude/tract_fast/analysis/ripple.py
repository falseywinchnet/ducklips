"""Fixed colorations: the fine structure of the long-term spectrum of voiced
loud frames, after removing its 1/3-octave trend. A resonance or notch that
stays put whatever the phone survives the averaging as a ripple; moving
formants average out. Prints the largest ripples (Hz, dB) and the ripple's
RMS over 300 Hz-8 kHz for each group.
usage: ripple.py LABEL GLOB [LABEL GLOB ...]"""
import sys, glob, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly

RATE = 32000


def ltas(paths):
    n = int(0.064 * RATE)
    acc = numpy.zeros(n // 2 + 1)
    for p in paths:
        r, x = wavfile.read(p)
        x = x.astype(float)
        if x.ndim > 1:
            x = x[:, 0]
        if r != RATE:
            x = resample_poly(x, RATE, r)
        frames = [x[i:i + n] for i in range(0, len(x) - n, n // 4)]
        if not frames:
            continue
        level = numpy.array([numpy.sqrt(numpy.mean(f ** 2)) for f in frames])
        top = level.max()
        for f, l in zip(frames, level):
            if l < top * 10 ** (-25 / 20):
                continue
            s = f - f.mean()
            c = numpy.correlate(s[: n // 2], s[: n // 2], "full")[n // 2 - 1:]
            lo, hi = int(RATE / 300), int(RATE / 60)
            if c[0] <= 0 or c[lo + numpy.argmax(c[lo:hi])] < 0.3 * c[0]:
                continue
            spec = numpy.abs(numpy.fft.rfft(s * numpy.hanning(n))) ** 2
            acc += spec / (spec.sum() + 1e-30)
    return numpy.fft.rfftfreq(n, 1 / RATE), 10 * numpy.log10(acc + 1e-30)


def trend(freqs, db):
    out = numpy.empty_like(db)
    for i, f in enumerate(freqs):
        lo, hi = max(f, 40) * 2 ** (-1 / 6), max(f, 40) * 2 ** (1 / 6)
        sel = (freqs >= lo) & (freqs <= hi)
        out[i] = 10 * numpy.log10(numpy.mean(10 ** (db[sel] / 10)))
    return out


for label, pattern in zip(sys.argv[1::2], sys.argv[2::2]):
    paths = sorted(glob.glob(pattern))
    freqs, db = ltas(paths)
    ripple = db - trend(freqs, db)
    band = (freqs >= 300) & (freqs <= 7000)
    idx = numpy.where(band)[0]
    order = idx[numpy.argsort(-numpy.abs(ripple[idx]))]
    picked = []
    for i in order:
        if all(abs(freqs[i] - freqs[j]) > 150 for j in picked):
            picked.append(i)
        if len(picked) == 8:
            break
    print("%-8s files %3d  ripple rms %.2f dB  top: %s" % (
        label, len(paths), numpy.sqrt(numpy.mean(ripple[band] ** 2)),
        " ".join("%d:%+.1f" % (freqs[i], ripple[i]) for i in sorted(picked))))
