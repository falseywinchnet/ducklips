"""Pitch against the plan pooled over many renders: for every planned vowel
in each SCOREDIR/*.score with WAVDIR/<same>.wav, F0 every 15 ms (40 ms
window, autocorrelation >= 0.4); prints the vowel frames voiced, median
|error| (cents), share within 100 cents, median within-vowel swing, and the
median onset error (first voiced frame of a vowel after silence).
usage: pe29.py SCOREDIR WAVDIR [WAVDIR...]"""
import sys, os, re, glob, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly


def load(path):
    r, x = wavfile.read(path)
    return resample_poly(x.astype(float), 16000, r)


def f0(x, t):
    seg = x[int(t * 16000):int((t + 0.04) * 16000)]
    seg = seg - seg.mean()
    c = numpy.correlate(seg, seg, "full")[len(seg) - 1:]
    lo, hi = 53, 267
    k = lo + numpy.argmax(c[lo:hi])
    return 16000 / k if c[0] > 0 and c[k] > 0.4 * c[0] else 0.0


scores = sorted(glob.glob(os.path.join(sys.argv[1], "*.score")))
for wavdir in sys.argv[2:]:
    errors, swings, frames, voiced, onsets = [], [], 0, 0, []
    for score in scores:
        wav = os.path.join(wavdir, os.path.basename(score)[:-6] + ".wav")
        if not os.path.exists(wav):
            continue
        x = load(wav)
        rows = [m for m in (re.match(r"# ([A-Z]+)-?(\d) ([\d.]+) ([\d.]+) ([\d.]+)", l) for l in open(score)) if m]
        prev_voiced = False
        for m in rows:
            s, d, plan = float(m.group(3)), float(m.group(4)), float(m.group(5))
            if plan <= 0:
                prev_voiced = m.group(1) in ("M", "N", "NG", "L", "R", "W", "Y")
                continue
            track = []
            for t in numpy.arange(s, s + d - 0.03, 0.015):
                frames += 1
                v = f0(x, t)
                if v > 0:
                    voiced += 1
                    track.append(1200 * numpy.log2(v / plan))
            errors += [abs(v) for v in track]
            if len(track) >= 3:
                swings.append(max(track) - min(track))
            if track and not prev_voiced:
                onsets.append(track[0])
            prev_voiced = True
    e = numpy.array(errors)
    print("%-22s voiced %4.1f%%  |error| %3.0f c  within100 %3.0f%%  swing %3.0f c  onset %+4.0f c (|%3.0f|)" % (
        wavdir.rstrip("/").split("/")[-1], 100 * voiced / max(frames, 1), numpy.median(e), 100 * numpy.mean(e < 100),
        numpy.median(swings), numpy.median(onsets), numpy.median(numpy.abs(onsets))))
