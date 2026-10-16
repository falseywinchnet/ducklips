"""Cycle-to-cycle coherence by band: cycles from waveform matching
(wmatch.cycles); for consecutive cycle pairs (broadband corr > 0.9) in loud
stretches, the band-passed signal's normalized correlation between one
period and the next; HNR = 10 log10(r / (1 - r)), median per band.
usage: bandhnr.py FILE..."""
import sys, numpy, librosa
sys.path.insert(0, "/tmp")
from wmatch import cycles
from scipy.signal import butter, sosfiltfilt
BANDS = [(80, 1000), (1000, 2000), (2000, 4000), (4000, 7000)]
for p in sys.argv[1:]:
    y, sr = librosa.load(p, sr=16000)
    rs = cycles(y, sr)
    lev = numpy.sqrt(numpy.convolve(y ** 2, numpy.ones(320) / 320, "same")); loud = lev.max() * 10 ** (-15 / 20)
    out = []
    for lo, hi in BANDS:
        yb = sosfiltfilt(butter(4, [lo, min(hi, 7900)], btype="band", fs=sr, output="sos"), y); h = []
        for r in rs:
            for (t0, per, c), (t1, _, _) in zip(r[:-1], r[1:]):
                if c < 0.9: continue
                a = int(t0 * sr); b = int(t1 * sr); P = b - a
                if lev[a] < loud or b + P > len(y): continue
                u, v = yb[a:b], yb[b:b + P]
                rr = numpy.dot(u, v) / numpy.sqrt(numpy.dot(u, u) * numpy.dot(v, v) + 1e-20)
                rr = min(max(rr, 1e-3), 0.9999); h.append(10 * numpy.log10(rr / (1 - rr)))
        out.append(numpy.median(h) if h else float("nan"))
    print("%-30s HNR dB by band 0.08-1 / 1-2 / 2-4 / 4-7 kHz: %s" % (p.split("/")[-1][:30], "  ".join("%5.1f" % v for v in out)))
