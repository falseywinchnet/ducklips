"""Pitch contour shape: pyin F0 (5 ms, 0.02 st resolution), 25 ms median
smoothing; per second of voicing: reversals (turns of at least H st), mean
|slope| st/s, and the jump across unvoiced gaps under 150 ms (median |st|).
usage: contour.py H FILE..."""
import sys, numpy, librosa, scipy.signal
H = float(sys.argv[1])
def track(p):
    y, sr = librosa.load(p, sr=16000)
    f0, vf, vp = librosa.pyin(y, fmin=60, fmax=320, sr=sr, frame_length=1024, hop_length=80, resolution=0.02)
    return 12 * numpy.log2(f0 / 100)
def runs(st):
    v = numpy.isfinite(st); out = []; i = 0
    while i < len(st):
        if not v[i]: i += 1; continue
        j = i
        while j < len(st) and v[j]: j += 1
        out.append((i, j)); i = j
    return out
for p in sys.argv[2:]:
    st = track(p); rs = [r for r in runs(st) if r[1] - r[0] >= 20]
    rev, slope, total, gaps = 0, [], 0, []
    for a, b in rs:
        s = scipy.signal.medfilt(st[a:b], 5); total += (b - a) * 0.005
        slope += list(numpy.abs(numpy.diff(s)) / 0.005)
        ext = s[0]; d = 0
        for x in s[1:]:
            if d >= 0 and x > ext: ext = x; d = 1 if d == 1 or x - s[0] > H else d
            if d <= 0 and x < ext: ext = x; d = -1 if d == -1 or s[0] - x > H else d
            if d == 1 and ext - x >= H: rev += 1; d = -1; ext = x
            elif d == -1 and x - ext >= H: rev += 1; d = 1; ext = x
    for (a, b), (c, e) in zip(rs[:-1], rs[1:]):
        if (c - b) * 0.005 < 0.15: gaps.append(abs(numpy.median(st[c:c + 4]) - numpy.median(st[b - 4:b])))
    print("%-34s voiced %5.2fs reversals/s %.2f  |slope| %.1f st/s  gap jump %.1f st" % (p.split("/")[-1][:34], total, rev / total, numpy.mean(slope), numpy.median(gaps) if gaps else 0))
