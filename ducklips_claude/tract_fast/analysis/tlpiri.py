"""Ideal transmission-line effect of the piriform recesses on F1-F3: the
model's area profile for a posture, with and without a closed side branch
(2 cm deep, 1.6 -> 0.2 cm^2 times SCALE, through a one-node port of PORT
cm^2) at 3.0 cm above the glottis. usage: tlpiri.py POSTURES.txt vowel..."""
import sys, subprocess, os, numpy
from scipy.signal import find_peaks
sys.argv_ = sys.argv
exec(open("/tmp/tlcheck.py").read().split("def engine")[0].replace('story = {}', 'story = {}').split("def tl_formants")[0])
exec("def profile" + open("/tmp/tlcheck.py").read().split("def profile")[1].split("def engine")[0])
def chain(lengths, areas, k):
    M = numpy.eye(2, dtype=complex)
    for L, A in zip(lengths, areas):
        Z = RHO * C / A; kl = k * L
        M = M @ numpy.array([[numpy.cos(kl), 1j * Z * numpy.sin(kl)], [1j * numpy.sin(kl) / Z, numpy.cos(kl)]])
    return M
def formants(lengths, areas, scale, port, at=3.0):
    edges = numpy.cumsum(lengths); fr = numpy.arange(100, 5000, 5.0); H = []
    cut = int(numpy.searchsorted(edges, at))
    bl = numpy.full(6, 2.0 / 6); ba = scale * (1.6 * (1 - (numpy.arange(6) + 0.5) / 6) + 0.2)
    for f in fr:
        w = 2 * numpy.pi * f; k = w / C - 1j * 0.0004 * numpy.sqrt(f)
        a = areas[-1]; rad = numpy.sqrt(a / numpy.pi)
        Zr = RHO * C / a * ((k.real * rad) ** 2 / 2 + 1j * 8 * k.real * rad / (3 * numpy.pi))
        low = chain(lengths[:cut], areas[:cut], k); high = chain(lengths[cut:], areas[cut:], k)
        if scale > 0:
            # branch input impedance: port duct then the closed recess
            Mb = chain(numpy.r_[0.365, bl], numpy.r_[port, ba], k)
            Zb = Mb[0, 0] / Mb[1, 0]          # closed far end: U=0 -> Z = A/C
            S = numpy.array([[1, 0], [1 / Zb, 1]], dtype=complex)
            M = low @ S @ high
        else:
            M = low @ high
        H.append(1 / abs(M[1, 0] * Zr + M[1, 1]))
    H = 20 * numpy.log10(numpy.array(H)); pk, _ = find_peaks(H)
    return [fr[i] for i in pk[:3]]
post = {}
for l in open(sys.argv_[1]):
    f = l.split("#")[0].split()
    if f and not l.startswith("#") and "=" in l and not f[0].startswith(("rule", "anat")):
        post[f[0].lower()] = " ".join(f[1:])
print("vowel   none              x1                x0.5              x0.25")
for v in sys.argv_[2:]:
    L, A = profile(post[v])
    print("%-4s " % v + "  ".join("%4.0f %4.0f %4.0f" % tuple(formants(L, A, s, 0.6)) for s in (0, 1, 0.5, 0.25)))
