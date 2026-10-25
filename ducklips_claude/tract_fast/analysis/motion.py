"""How speech moves: F0, F1 and F2 trajectories (every 5 ms, voiced runs of
at least 100 ms) and the distribution of their speed in semitones per
second. For each track: the share of time it sits still (|v| < 4 st/s),
moves gently (4-40), or jumps (> 40), and the median |acceleration|
(st/s^2) as the curvature of the path.
Formants from LPC (order 12 at 10 kHz, 25 ms Hamming), continuity by
nearest-candidate tracking, medians of 3 frames.
usage: motion.py LABEL GLOB[:start:end] [LABEL GLOB ...]"""
import sys, glob, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly, lfilter, medfilt
from scipy.linalg import solve_toeplitz

RATE = 10000
HOP = 50


def load(spec):
    parts = spec.split(":")
    out = []
    for p in sorted(glob.glob(parts[0])):
        r, x = wavfile.read(p); x = x.astype(float)
        if x.ndim > 1: x = x[:, 0]
        x = resample_poly(x, RATE, r)
        if len(parts) == 3:
            x = x[int(float(parts[1]) * RATE):int(float(parts[2]) * RATE)]
        out.append(x / (numpy.abs(x).max() + 1e-9))
    return out


def tracks(x):
    win = 250
    pre = lfilter([1, -0.7], [1], x)
    f0, f1, f2 = [], [], []
    for i in range(0, len(x) - win - 200, HOP):
        a = x[i:i + win] - x[i:i + win].mean()
        full = numpy.correlate(x[i:i + win + 150] - x[i:i + win + 150].mean(), a, "valid")[::-1]
        c = numpy.correlate(a, a, "full")[win - 1:]
        lo, hi = 33, 143
        if c[0] <= 1e-9:
            f0.append(0); f1.append(0); f2.append(0); continue
        k = lo + int(numpy.argmax(c[lo:hi]))
        if c[k] < 0.5 * c[0]:
            f0.append(0); f1.append(0); f2.append(0); continue
        # parabolic refinement
        if lo < k < hi - 1:
            y0, y1, y2 = c[k - 1], c[k], c[k + 1]
            den = y0 - 2 * y1 + y2
            k = k + (0.5 * (y0 - y2) / den if den < 0 else 0)
        f0.append(RATE / k)
        seg = pre[i:i + win] * numpy.hamming(win)
        cc = numpy.correlate(seg, seg, "full")[win - 1:win + 13]
        ac = solve_toeplitz(cc[:12], -cc[1:13])
        roots = numpy.roots(numpy.r_[1, ac]); roots = roots[numpy.imag(roots) > 0]
        fr = numpy.angle(roots) * RATE / (2 * numpy.pi); bw = -numpy.log(numpy.abs(roots)) * RATE / numpy.pi
        cand = sorted(f for f, b in zip(fr, bw) if 200 < f < 3500 and b < 500)
        f1.append(cand[0] if len(cand) > 0 else 0); f2.append(cand[1] if len(cand) > 1 else 0)
    return numpy.array(f0), numpy.array(f1), numpy.array(f2)


def stats(values):
    v = numpy.concatenate(values) if values else numpy.array([0.0])
    return (numpy.mean(numpy.abs(v) < 4), numpy.mean((numpy.abs(v) >= 4) & (numpy.abs(v) <= 40)), numpy.mean(numpy.abs(v) > 40))


for label, spec in zip(sys.argv[1::2], sys.argv[2::2]):
    vel = {"F0": [], "F1": [], "F2": []}; acc = {"F0": [], "F1": [], "F2": []}
    for x in load(spec):
        f0, f1, f2 = tracks(x)
        voiced = (f0 > 0) & (f1 > 0) & (f2 > 0)
        i = 0
        while i < len(voiced):
            if not voiced[i]:
                i += 1; continue
            j = i
            while j < len(voiced) and voiced[j]: j += 1
            if j - i >= 20:
                for name, tr in (("F0", f0[i:j]), ("F1", f1[i:j]), ("F2", f2[i:j])):
                    st = 12 * numpy.log2(medfilt(tr, 3) / numpy.median(tr))
                    st = numpy.convolve(st, numpy.ones(3) / 3, "valid")
                    v = numpy.diff(st) / (HOP / RATE)
                    if numpy.max(numpy.abs(numpy.diff(st))) > 6:  # a tracking jump of half an octave: skip the run
                        continue
                    vel[name].append(v); acc[name].append(numpy.diff(v) / (HOP / RATE))
            i = j
    row = []
    for name in ("F0", "F1", "F2"):
        still, gentle, jump = stats(vel[name])
        a = numpy.median(numpy.abs(numpy.concatenate(acc[name]))) if acc[name] else 0
        row.append("%s still %.2f gentle %.2f jump %.2f curv %4.0f" % (name, still, gentle, jump, a))
    print("%-10s %s" % (label, " | ".join(row)))
