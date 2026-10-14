"""Fast modulation ("a vibrator on the throat"): in loud voiced stretches,
the amplitude envelope (cycle-rate smoothing removed: RMS over 2 ms, then
the part between 10 and 150 Hz) and its spectrum, as dB re the envelope's
mean, in bands; also the cycle-to-cycle peak-amplitude variation (shimmer,
dB) and period variation (jitter, %) from glottal pulses picked as the
waveform's main negative peaks.
usage: vibrator.py LABEL WAV[:start:end] [LABEL WAV[:start:end] ...]"""
import sys, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, butter, sosfiltfilt

RATE = 16000
BANDS = [(10, 20), (20, 40), (40, 70), (70, 110), (110, 150)]


def load(spec):
    parts = spec.split(":")
    r, x = wavfile.read(parts[0]); x = x.astype(float)
    if x.ndim > 1: x = x[:, 0]
    x = resample_poly(x, RATE, r)
    if len(parts) == 3:
        x = x[int(float(parts[1]) * RATE):int(float(parts[2]) * RATE)]
    return x


def stretches(x):
    n = 320
    lev = numpy.array([numpy.sqrt(numpy.mean(x[i:i + n] ** 2)) for i in range(0, len(x) - n, n)])
    loud = lev > lev.max() * 10 ** (-12 / 20)
    i = 0
    while i < len(loud):
        if not loud[i]:
            i += 1; continue
        j = i
        while j < len(loud) and loud[j]: j += 1
        if j - i >= 8:  # 160 ms
            yield x[i * n:j * n]
        i = j


for label, spec in zip(sys.argv[1::2], sys.argv[2::2]):
    x = load(spec)
    power = numpy.zeros(len(BANDS)); count = 0
    shim, jit = [], []
    for seg in stretches(x):
        env = numpy.sqrt(numpy.convolve(seg ** 2, numpy.ones(32) / 32, "same"))
        env = resample_poly(env, 1, 16)  # 1 kHz
        sos = butter(4, [10, 150], btype="band", fs=1000, output="sos")
        mod = sosfiltfilt(sos, env) / (env.mean() + 1e-9)
        s = numpy.abs(numpy.fft.rfft(mod * numpy.hanning(len(mod)))) ** 2 / len(mod)
        f = numpy.fft.rfftfreq(len(mod), 1 / 1000)
        power += numpy.array([s[(f >= a) & (f < b)].sum() for a, b in BANDS]); count += 1
        # pulses: local minima below 60% of the stretch's minimum, at least 3 ms apart
        thr = 0.5 * seg.min()
        idx = [k for k in range(1, len(seg) - 1) if seg[k] < thr and seg[k] <= seg[k - 1] and seg[k] <= seg[k + 1]]
        picked = []
        for k in idx:
            if not picked or k - picked[-1] > 48:
                picked.append(k)
            elif seg[k] < seg[picked[-1]]:
                picked[-1] = k
        amps = numpy.abs(seg[picked]); per = numpy.diff(picked)
        good = (per > 40) & (per < 220)
        if good.sum() > 5:
            jit.append(numpy.mean(numpy.abs(numpy.diff(per[good]))) / numpy.mean(per[good]) * 100)
            shim.append(numpy.mean(numpy.abs(20 * numpy.log10(amps[1:] / amps[:-1]))))
    power /= max(count, 1)
    print("%-10s stretches %2d  envelope modulation dB: %s   shimmer %.2f dB  jitter %.2f%%" % (
        label, count, " ".join("%d-%d:%.1f" % (a, b, 10 * numpy.log10(p + 1e-12)) for (a, b), p in zip(BANDS, power)),
        numpy.median(shim) if shim else 0, numpy.median(jit) if jit else 0))
