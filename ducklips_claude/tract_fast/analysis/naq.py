"""Glottal flow from the audio by iterative adaptive inverse filtering
(IAIF, Alku 1992; at 16 kHz), then per cycle the normalized amplitude
quotient NAQ = (flow peak-to-peak / largest closing rate) / period, and the
closed share of the cycle (flow within 10% of its minimum): pressed voices
have low NAQ, breathy high (male modal about 0.1). The same for any voice.
usage: naq.py ours WAV:START:END ...  |  naq.py arctic LABDIR WAVDIR VOWELS MAXFILES"""
import os, sys, numpy, soundfile, librosa
from scipy.signal import lfilter, butter, sosfiltfilt

RATE = 16000

def lpc(x, order):
    return librosa.lpc(x * numpy.hanning(len(x)), order=order)

def iaif(x):
    x = sosfiltfilt(butter(4, 40, "highpass", fs=RATE, output="sos"), x)
    p = RATE // 1000 + 2
    g1 = lpc(x, 1); y = lfilter(g1, [1.0], x)
    v1 = lpc(y, p); d = lfilter(v1, [1.0], x)
    gi = lfilter([1.0], [1.0, -0.99], d)
    g2 = lpc(gi, 4); y2 = lfilter(g2, [1.0], x)
    v2 = lpc(y2, p); d2 = lfilter(v2, [1.0], x)
    flow = lfilter([1.0], [1.0, -0.99], d2)
    return flow, d2

def measure(x, r):
    x = librosa.resample(x, orig_sr=r, target_sr=RATE)
    f0, v, _ = librosa.pyin(x, fmin=70, fmax=260, sr=RATE, frame_length=1024, hop_length=128)
    if numpy.isfinite(f0).sum() < 3: return None
    period = RATE / numpy.nanmedian(f0)
    flow, deriv = iaif(x)
    naqs, closed = [], []
    i = int(numpy.argmin(deriv[:int(period)]))
    while i + 1.4 * period < len(deriv):
        lo, hi = int(i + 0.6 * period), int(i + 1.4 * period)
        j = lo + int(numpy.argmin(deriv[lo:hi]))
        cyc = flow[i:j]; dc = deriv[i:j]
        if len(cyc) > 10:
            ac = cyc.max() - cyc.min(); dpeak = -dc.min()
            if dpeak > 0:
                naqs.append(ac / dpeak / len(cyc))
                closed.append(numpy.mean(cyc < cyc.min() + 0.1 * ac))
        i = j
    if len(naqs) < 4: return None
    return float(numpy.median(naqs)), float(numpy.median(closed))

def main():
    if sys.argv[1] == "ours":
        for spec in sys.argv[2:]:
            path, a, b = spec.split(":")
            x, r = soundfile.read(path)
            res = measure(x[int(float(a) * r):int(float(b) * r)], r)
            print("%-14s NAQ %.3f  closed share %.2f" % (path.split("/")[-1], *res))
        return
    labdir, wavdir, vowels, maxfiles = sys.argv[2], sys.argv[3], sys.argv[4].split(","), int(sys.argv[5])
    rows = {v: [] for v in vowels}
    for name in sorted(os.listdir(labdir))[:maxfiles]:
        lines = [l.split() for l in open(os.path.join(labdir, name)) if len(l.split()) == 3]
        wav = os.path.join(wavdir, name.replace(".lab", ".wav"))
        if not os.path.exists(wav): continue
        x, r = soundfile.read(wav); prev = 0.0
        for end, _, ph in lines:
            end = float(end)
            if ph in rows and end - prev >= 0.12:
                mid = (prev + end) / 2; half = (end - prev) * 0.35
                res = measure(x[int((mid - half) * r):int((mid + half) * r)], r)
                if res: rows[ph].append(res)
            prev = end
    for v in vowels:
        a = numpy.array(rows[v])
        print("bdl %-3s n=%-3d NAQ %.3f (quartiles %.3f-%.3f)  closed share %.2f" % (v, len(a), numpy.median(a[:, 0]), numpy.percentile(a[:, 0], 25), numpy.percentile(a[:, 0], 75), numpy.median(a[:, 1])))

if __name__ == "__main__":
    main()
