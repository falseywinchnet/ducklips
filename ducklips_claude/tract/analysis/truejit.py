"""True cycle irregularity of the sheet's folds from a per-sample trace (glottal
area column 6): closures as area falling under 0.002 cm^2; per voiced run
(periods 3-15 ms), jitter (mean |dT| / T, %), the share of cycles alternating
long-short (period doubling: successive differences changing sign with
|dT| > 3%), and shimmer of peak flow (dB).
usage: truejit.py TRACE [rate]"""
import sys, numpy
rate = float(sys.argv[2]) if len(sys.argv) > 2 else 96000.0
A = []; U = []
for l in open(sys.argv[1]):
    f = l.split(); U.append(float(f[4])); A.append(float(f[5]))
A = numpy.array(A); U = numpy.array(U)
closed = A < 0.002
starts = numpy.where(closed[1:] & ~closed[:-1])[0]
T = numpy.diff(starts) / rate
peaks = numpy.array([U[a:b].max() for a, b in zip(starts[:-1], starts[1:])])
ok = (T > 0.003) & (T < 0.015)
jit, alt, shim = [], [], []
i = 0
while i < len(T):
    if not ok[i]:
        i += 1; continue
    j = i
    while j < len(T) and ok[j]: j += 1
    t = T[i:j]; p = peaks[i:j]
    if len(t) >= 6:
        d = numpy.diff(t)
        jit.append(numpy.mean(numpy.abs(d)) / numpy.mean(t) * 100)
        big = numpy.abs(d) > 0.03 * numpy.mean(t)
        alt.append(numpy.mean((numpy.sign(d[1:]) != numpy.sign(d[:-1])) & big[1:] & big[:-1]))
        shim.append(numpy.mean(numpy.abs(20 * numpy.log10((p[1:] + 1) / (p[:-1] + 1)))))
    i = j
print("runs %d  jitter %.2f%%  alternating %.2f  shimmer %.2f dB" % (len(jit), numpy.median(jit), numpy.median(alt), numpy.median(shim)))
