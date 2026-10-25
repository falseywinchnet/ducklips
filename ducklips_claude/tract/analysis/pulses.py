"""Cycle-to-cycle jitter and shimmer from the audio alone, the same way for
any voice: the signal low-passed (1.2 kHz); within each expected period
(pyin f0) the most negative peak marks a cycle; jitter = mean |period
change| / mean period, shimmer = mean |peak change| / mean peak (local, %).
usage: pulses.py ours WAV:START:END ...  |  pulses.py arctic LABDIR WAVDIR VOWELS MAXFILES"""
import os, sys, numpy, soundfile, librosa
from scipy.signal import butter, sosfiltfilt

def cycles(x, r):
    f0, v, _ = librosa.pyin(x, fmin=70, fmax=260, sr=r, frame_length=2048, hop_length=256)
    if numpy.isfinite(f0).sum() < 3:
        return None
    period = r / numpy.nanmedian(f0)
    y = sosfiltfilt(butter(4, 1200, "lowpass", fs=r, output="sos"), x)
    marks, peaks = [], []
    i = int(numpy.argmin(y[:int(period)]))
    while i + 1.4 * period < len(y):
        lo, hi = int(i + 0.6 * period), int(i + 1.4 * period)
        j = lo + int(numpy.argmin(y[lo:hi]))
        marks.append(j); peaks.append(-y[j]); i = j
    if len(marks) < 6:
        return None
    per = numpy.diff(marks).astype(float); amp = numpy.array(peaks)
    return (100 * numpy.mean(numpy.abs(numpy.diff(per))) / per.mean(), 100 * numpy.mean(numpy.abs(numpy.diff(amp))) / amp.mean())

def main():
    if sys.argv[1] == "ours":
        for spec in sys.argv[2:]:
            path, a, b = spec.split(":")
            x, r = soundfile.read(path)
            res = cycles(x[int(float(a) * r):int(float(b) * r)], r)
            print("%-12s jitter %.2f%%  shimmer %.1f%%" % (path.split("/")[-1], *res))
    else:
        labdir, wavdir, vowels, maxfiles = sys.argv[2], sys.argv[3], sys.argv[4].split(","), int(sys.argv[5])
        rows = {v: [] for v in vowels}
        for name in sorted(os.listdir(labdir))[:maxfiles]:
            lines = [l.split() for l in open(os.path.join(labdir, name)) if len(l.split()) == 3]
            wav = os.path.join(wavdir, name.replace(".lab", ".wav"))
            if not os.path.exists(wav): continue
            x, r = soundfile.read(wav)
            prev = 0.0
            for end, _, ph in lines:
                end = float(end)
                if ph in rows and end - prev >= 0.12:
                    mid = (prev + end) / 2; half = (end - prev) * 0.35
                    res = cycles(x[int((mid - half) * r):int((mid + half) * r)], r)
                    if res: rows[ph].append(res)
                prev = end
        for v in vowels:
            a = numpy.array(rows[v])
            print("bdl %-3s n=%-3d jitter %.2f%%  shimmer %.1f%%" % (v, len(a), numpy.median(a[:, 0]), numpy.median(a[:, 1])))


if __name__ == "__main__":
    main()
