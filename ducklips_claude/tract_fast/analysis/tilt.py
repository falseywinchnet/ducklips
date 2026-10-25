"""Two things the owner's recording shows against ours, per file:
  vowel tilt   in periodic frames (autocorrelation >= 0.6, within 20 dB of
               the loudest), energy at 2.5-5 kHz re 0-1.2 kHz (dB), median;
  bursts       the sharpest onsets: frames whose 3-8 kHz energy rises by
               over 20 dB within 5 ms; their 2 ms peak level re the median
               voiced frame level (dB), median over bursts, and count per
               second of speech.
usage: tilt.py LABEL WAV[:start:end] [...]"""
import sys, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, butter, sosfiltfilt

RATE = 16000


def load(spec):
    parts = spec.split(":")
    r, x = wavfile.read(parts[0]); x = x.astype(float)
    if x.ndim > 1: x = x[:, 0]
    x = resample_poly(x, RATE, r)
    if len(parts) == 3:
        x = x[int(float(parts[1]) * RATE):int(float(parts[2]) * RATE)]
    return x


for label, spec in zip(sys.argv[1::2], sys.argv[2::2]):
    x = load(spec)
    n = 480; hop = 160
    frames = [x[i:i + n] for i in range(0, len(x) - n, hop)]
    lev = numpy.array([numpy.sqrt(numpy.mean(f ** 2)) + 1e-9 for f in frames])
    top = lev.max()
    tilts, voiced_lev = [], []
    for f, l in zip(frames, lev):
        if l < top * 0.1:
            continue
        a = f - f.mean(); c = numpy.correlate(a, a, "full")[n - 1:]
        k = 53 + int(numpy.argmax(c[53:230]))
        if c[0] <= 0 or c[k] < 0.6 * c[0]:
            continue
        s = numpy.abs(numpy.fft.rfft(a * numpy.hanning(n))) ** 2; fr = numpy.fft.rfftfreq(n, 1 / RATE)
        tilts.append(10 * numpy.log10(s[(fr >= 2500) & (fr < 5000)].sum() / s[fr < 1200].sum()))
        voiced_lev.append(20 * numpy.log10(l))
    hb = sosfiltfilt(butter(4, [3000, 7900], btype="band", fs=RATE, output="sos"), x)
    e = 10 * numpy.log10(numpy.convolve(hb ** 2, numpy.ones(32) / 32, "same") + 1e-12)  # 2 ms
    bursts = []
    i = 80
    while i < len(e) - 80:
        if e[i] - e[i - 80] > 20 and e[i] > e.max() - 40:
            peak = e[i:i + 80].max()
            bursts.append(peak)
            i += 800
        else:
            i += 8
    vmed = numpy.median(voiced_lev) if voiced_lev else 0
    speech = numpy.sum(lev > top * 0.03) * hop / RATE
    print("%-10s vowel tilt %+.1f dB (n %d)  bursts %d (%.1f/s), peak 2ms re voiced frame median %+.1f dB" % (
        label, numpy.median(tilts), len(tilts), len(bursts), len(bursts) / max(speech, 1e-9),
        numpy.median(bursts) - vmed if bursts else numpy.nan))
