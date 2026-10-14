"""Mean spectrum of the middle half of each instance of the given phones,
ours (out/arctic, our scores) against bdl (ref/arctic_bdl labels), over
ARCTIC a0001-a0060; band levels in dB re the utterance's overall level.
usage: phone_spectrum.py PHONE..."""
import os, sys, numpy
sys.path.insert(0, os.path.expanduser("~/ducklips/voice/tools"))
import compare_arctic as ca
BANDS = [(80, 200), (200, 400), (400, 700), (700, 1000), (1000, 1500), (1500, 2200), (2200, 3000), (3000, 4500), (4500, 7000)]
def bands(x, r):
    n = len(x); w = numpy.hanning(n); S = numpy.abs(numpy.fft.rfft(x * w)) ** 2 / (w ** 2).sum(); f = numpy.fft.rfftfreq(n, 1 / r)
    return numpy.array([S[(f >= a) & (f < b)].sum() for a, b in BANDS])
for phone in sys.argv[1:]:
    rows = {"ours": [], "bdl": []}
    for number in range(1, 61):
        name = "arctic_a%04d" % number
        for who in ("ours", "bdl"):
            if who == "ours":
                path = os.path.join(ca.OUT, name + "_ours.wav"); phones = ca.our_phones(os.path.join(ca.OUT, name + ".score"))
            else:
                path = os.path.join(ca.REF, "wav", name + ".wav"); phones = ca.reference_phones(name)
            if not os.path.exists(path): continue
            x, r = ca.read_wave(path); x = numpy.asarray(x, float)
            level = numpy.mean(x ** 2) + 1e-20
            for p, s, e in phones:
                if p != phone or e - s < 0.03: continue
                a, b = int((s + 0.25 * (e - s)) * r), int((e - 0.25 * (e - s)) * r)
                rows[who].append(bands(x[a:b], r) / (b - a) / level)
    for who in ("bdl", "ours"):
        m = numpy.mean(rows[who], axis=0) if rows[who] else numpy.zeros(len(BANDS)) + 1e-20
        print("%-3s %-4s n=%3d " % (phone, who, len(rows[who])) + " ".join("%6.1f" % (10 * numpy.log10(v + 1e-20)) for v in m))
print("bands: " + " ".join("%d-%d" % b for b in BANDS))
