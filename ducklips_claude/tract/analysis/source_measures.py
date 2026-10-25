"""Voice-source measures on vowels (tract shape largely factored out by
comparing like vowels):
  tilt   - slope of the harmonic peaks' levels (dB per octave), 0.5-4 kHz;
  valley - formant contrast: 90th minus 10th percentile of the detrended
           harmonic-peak envelope in 0.8-4.5 kHz (dB);
  h1h2   - first harmonic minus second (dB): high = soft, open closure;
  hnr    - harmonic peaks over the troughs between them, 1.5-4 kHz (dB);
  jitter - mean absolute period change / mean period (%), from pyin f0;
  flutter- f0 standard deviation over the segment after removing its linear
           trend (Hz).
usage: source_measures.py ours WAV:START:END:VOWEL ...
       source_measures.py arctic LABDIR WAVDIR VOWEL[,VOWEL] MAXFILES"""
import os, sys, numpy, soundfile, librosa

def measure(x, r):
    f0, voiced, _ = librosa.pyin(x, fmin=60, fmax=300, sr=r, frame_length=1024, hop_length=128)
    f0v = f0[numpy.isfinite(f0)]
    if len(f0v) < 5:
        return None
    pitch = float(numpy.median(f0v))
    period = 1.0 / f0v
    jitter = 100.0 * numpy.mean(numpy.abs(numpy.diff(period))) / numpy.mean(period)
    trend = numpy.polyval(numpy.polyfit(numpy.arange(len(f0v)), f0v, 1), numpy.arange(len(f0v)))
    flutter = float(numpy.std(f0v - trend))
    window = numpy.hanning(len(x))
    n = 1 << int(numpy.ceil(numpy.log2(len(x) * 4)))
    spectrum = 20 * numpy.log10(numpy.abs(numpy.fft.rfft(x * window, n)) + 1e-9)
    freqs = numpy.fft.rfftfreq(n, 1.0 / r)
    peaks, troughs, peak_f = [], [], []
    k = 1
    while k * pitch < 5000:
        lo, hi = (k - 0.3) * pitch, (k + 0.3) * pitch
        m = (freqs >= lo) & (freqs <= hi)
        peaks.append(spectrum[m].max()); peak_f.append(k * pitch)
        m2 = (freqs >= (k + 0.35) * pitch) & (freqs <= (k + 0.65) * pitch)
        troughs.append(spectrum[m2].min() if m2.any() else numpy.nan)
        k += 1
    peaks, troughs, peak_f = numpy.array(peaks), numpy.array(troughs), numpy.array(peak_f)
    sel = (peak_f >= 500) & (peak_f <= 4000)
    tilt = numpy.polyfit(numpy.log2(peak_f[sel]), peaks[sel], 1)[0]
    sel2 = (peak_f >= 800) & (peak_f <= 4500)
    env = peaks[sel2] - numpy.polyval(numpy.polyfit(numpy.log2(peak_f[sel2]), peaks[sel2], 1), numpy.log2(peak_f[sel2]))
    valley = numpy.percentile(env, 90) - numpy.percentile(env, 10)
    h1h2 = peaks[0] - peaks[1]
    sel3 = (peak_f >= 1500) & (peak_f <= 4000)
    hnr = numpy.nanmean(peaks[sel3] - troughs[sel3])
    return dict(pitch=pitch, tilt=tilt, valley=valley, h1h2=h1h2, hnr=hnr, jitter=jitter, flutter=flutter)

def show(name, rows):
    keys = ("pitch", "tilt", "valley", "h1h2", "hnr", "jitter", "flutter")
    med = {k: numpy.median([row[k] for row in rows]) for k in keys}
    print("%-14s n=%-3d " % (name, len(rows)) + "  ".join("%s %.2f" % (k, med[k]) for k in keys))

def main():
    if sys.argv[1] == "ours":
        for spec in sys.argv[2:]:
            path, a, b, vowel = spec.split(":")
            x, r = soundfile.read(path)
            x = x[int(float(a) * r):int(float(b) * r)]
            # 0.25 s windows through the hold
            rows = [measure(x[i:i + int(0.25 * r)], r) for i in range(0, len(x) - int(0.25 * r), int(0.25 * r))]
            show("ours " + vowel, [row for row in rows if row])
    else:
        labdir, wavdir, vowels, maxfiles = sys.argv[2], sys.argv[3], sys.argv[4].split(","), int(sys.argv[5])
        by_vowel = {v: [] for v in vowels}
        for name in sorted(os.listdir(labdir))[:maxfiles]:
            lines = [l.split() for l in open(os.path.join(labdir, name)) if len(l.split()) == 3]
            wav = os.path.join(wavdir, name.replace(".lab", ".wav"))
            if not os.path.exists(wav):
                continue
            x, r = soundfile.read(wav)
            prev = 0.0
            for end, _, phone in lines:
                end = float(end)
                if phone in by_vowel and end - prev >= 0.12:
                    mid = (prev + end) / 2; half = (end - prev) * 0.3
                    seg = x[int((mid - half) * r):int((mid + half) * r)]
                    row = measure(seg, r)
                    if row: by_vowel[phone].append(row)
                prev = end
        for v in vowels:
            if by_vowel[v]: show("%s %s" % (os.path.basename(os.path.dirname(labdir)), v), by_vowel[v])


if __name__ == "__main__":
    main()
