"""Cycle statistics from speak's CYCLES log: per voiced run (periods 2.5-15 ms,
at least 8 cycles), local jitter %, shimmer dB (sound peak-to-peak and peak
glottal flow), lag-1 correlation of period differences (-0.5 = white jitter,
toward -1 = cycle alternation), and rms of period and amplitude in the
3-30 Hz modulation band (after removing the trend below 3 Hz)."""
import sys, numpy
def runs(path):
    d = numpy.loadtxt(path)
    out, cur = [], []
    for row in d:
        if 2.5 < row[1] < 15:
            cur.append(row)
        else:
            if len(cur) >= 8: out.append(numpy.array(cur))
            cur = []
    if len(cur) >= 8: out.append(numpy.array(cur))
    return out
def band(t, v, lo=3.0, hi=30.0):
    # resample on 2 ms grid, FFT band rms relative to mean
    g = numpy.arange(t[0], t[-1], 0.002)
    if len(g) < 32: return None
    y = numpy.interp(g, t, v); m = y.mean(); y = y / m - 1
    f = numpy.fft.rfftfreq(len(y), 0.002); Y = numpy.fft.rfft(y * numpy.hanning(len(y)))
    sel = (f >= lo) & (f <= hi)
    return numpy.sqrt(2 * numpy.sum(numpy.abs(Y[sel]) ** 2)) / numpy.sqrt(numpy.sum(numpy.hanning(len(y)) ** 2)) * 100
for path in sys.argv[1:]:
    J, S, F, C, BP, BA, n = [], [], [], [], [], [], 0
    for r in runs(path):
        p, a, fl = r[:, 1], r[:, 3], r[:, 2]
        n += len(r)
        J.append(numpy.mean(numpy.abs(numpy.diff(p))) / p.mean() * 100)
        S.append(numpy.mean(numpy.abs(numpy.diff(20 * numpy.log10(a + 1e-12)))))
        F.append(numpy.mean(numpy.abs(numpy.diff(20 * numpy.log10(fl + 1e-3)))))
        dp = numpy.diff(p)
        if len(dp) > 3 and dp.std() > 0: C.append(numpy.corrcoef(dp[:-1], dp[1:])[0, 1])
        b = band(r[:, 0], p); ba = band(r[:, 0], a)
        if b is not None: BP.append(b); BA.append(ba)
    print("%-14s cycles %4d jitter %.2f%% shimmer %.2f dB flowshim %.2f dB dcorr %+.2f  mod3-30 period %.2f%% amp %.1f%%" % (
        path.split("/")[-1], n, numpy.median(J), numpy.median(S), numpy.median(F), numpy.median(C), numpy.median(BP), numpy.median(BA)))
