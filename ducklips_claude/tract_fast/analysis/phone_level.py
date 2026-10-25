"""Loudness of each phone against its sentence: RMS over the phone's middle
60% (dB re the sentence's loudest 20 ms), median per phone, for bdl and our
renders of the same sentences; prints phones with the difference.
usage: phone_level.py LABDIR BDLWAVDIR SCOREDIR OURWAVDIR [N]"""
import sys, os, numpy
sys.path.insert(0, "/tmp")
from vw import load
from dips import phones_bdl, phones_ours

ALIAS = {"ax": "ah", "axr": "er", "ix": "ih", "dx": "d", "el": "l", "em": "m", "en": "n", "nx": "n"}


def levels(source):
    out = {}
    for wav, phones in source:
        if not os.path.exists(wav):
            continue
        x = load(wav)
        n = 320
        frames = numpy.array([numpy.sqrt(numpy.mean(x[i:i + n] ** 2)) for i in range(0, len(x) - n, 80)])
        top = 20 * numpy.log10(frames.max() + 1e-9)
        for name, s, e in phones:
            name = ALIAS.get(name, name)
            if name in ("pau", "sil", "h#", "ssil") or e - s < 0.03:
                continue
            a, b = s + 0.2 * (e - s), e - 0.2 * (e - s)
            seg = x[int(a * 16000):int(b * 16000)]
            if len(seg) < 80:
                continue
            out.setdefault(name, []).append(20 * numpy.log10(numpy.sqrt(numpy.mean(seg ** 2)) + 1e-9) - top)
    return out


count = int(sys.argv[5]) if len(sys.argv) > 5 else 29
ref = levels(phones_bdl(sys.argv[1], sys.argv[2], count))
ours = levels(phones_ours(sys.argv[3], sys.argv[4]))
rows = []
for name in sorted(ref):
    if name in ours and len(ref[name]) >= 3 and len(ours[name]) >= 3:
        rows.append((numpy.median(ours[name]) - numpy.median(ref[name]), name, numpy.median(ref[name]),
                     numpy.median(ours[name]), len(ref[name]), len(ours[name])))
for d, name, r, o, nr, no in sorted(rows):
    print("%-3s bdl %6.1f  ours %6.1f  diff %+5.1f dB  (n %d/%d)" % (name, r, o, d, nr, no))
