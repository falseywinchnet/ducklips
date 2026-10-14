"""Cycle periods from audio by waveform matching: pyin finds voiced runs and
a coarse period; each next cycle start maximises the normalized correlation
of one period of waveform with the next (0.85-1.15 of the guess, parabolic
refinement). Reports RAP (3-cycle) median / p90 %, jumps (>5% off the
neighbours' mean... i.e. RAP > 3.33%) per 100 cycles, and the median
correlation of consecutive cycles.
usage: wmatch.py FILE [FILE...]   (DUMP=path writes time, period ms)"""
import sys, os, numpy, librosa
from scipy.signal import butter, sosfiltfilt
def cycles(y, sr):
    f0, vf, vp = librosa.pyin(y, fmin=60, fmax=320, sr=sr, frame_length=int(0.064 * sr), hop_length=int(0.005 * sr))
    hop = int(0.005 * sr)
    y = sosfiltfilt(butter(4, [60, 4000], btype="band", fs=sr, output="sos"), y)
    out = []
    i = 0
    while i < len(f0):
        if not numpy.isfinite(f0[i]): i += 1; continue
        j = i
        while j < len(f0) and numpy.isfinite(f0[j]): j += 1
        if j - i >= 12:
            t = i * hop; end = j * hop; run = []
            while True:
                k = min(int(t / hop), len(f0) - 1)
                guess = sr / (f0[k] if numpy.isfinite(f0[k]) else numpy.nanmedian(f0[i:j]))
                P = int(round(guess)); lo, hi = int(0.85 * guess), int(1.15 * guess) + 1
                if t + hi + P >= end: break
                a = y[t:t + P]; ea = numpy.dot(a, a)
                cs = numpy.array([numpy.dot(a, y[t + d:t + d + P]) / numpy.sqrt(ea * numpy.dot(y[t + d:t + d + P], y[t + d:t + d + P]) + 1e-30) for d in range(lo - 1, hi + 2)])
                m = int(numpy.argmax(cs[1:-1])) + 1
                den = cs[m - 1] - 2 * cs[m] + cs[m + 1]
                sh = 0.5 * (cs[m - 1] - cs[m + 1]) / den if den < 0 else 0.0
                tau = lo - 1 + m + sh
                run.append((t / sr, tau / sr * 1000, cs[m]))
                t += int(round(tau))
            if len(run) >= 5: out.append(numpy.array(run))
        i = j
    return out
for p in sys.argv[1:]:
    y, sr = librosa.load(p, sr=None)
    rs = cycles(y, sr); rap, cor = [], []
    dump = []
    for r in rs:
        per = r[:, 1]; m = (per[:-2] + per[1:-1] + per[2:]) / 3
        ok = (r[1:-1, 2] > 0.5)
        rap.extend(list((numpy.abs(per[1:-1] - m) / m * 100)[ok])); cor.extend(list(r[:, 2]))
        dump.extend(list(r))
    rap = numpy.array(rap)
    if os.environ.get("DUMP"): numpy.savetxt(os.environ["DUMP"], numpy.array(dump), fmt="%.5f")
    print("%-32s sr %5d cycles %5d RAP median %.2f%% p90 %.2f%% jumps %.1f/100 corr %.3f" % (
        p.split("/")[-1][:32], sr, len(rap), numpy.median(rap), numpy.percentile(rap, 90), 100 * numpy.mean(rap > 10 / 3), numpy.median(cor)))
