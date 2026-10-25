"""Pitch at voiced consonants between vowels (microprosody): for each vowel-
consonant-vowel with the consonant in CLASSES, F0 at the consonant's middle
against the mean of F0 at the two vowels' middles (cents). Prints per class
the count, median signed dip and median |dip|, and the 90th percentile |dip|.
usage: dips.py bdl LABDIR WAVDIR [N] | dips.py ours SCOREDIR WAVDIR"""
import sys, os, re, glob, numpy
sys.path.insert(0, "/tmp")
from vw import load, pitch_at, VOWELS

CLASSES = {"nasal": {"m", "n", "ng"}, "liquid": {"l", "r", "w", "y"}, "voiced": {"b", "d", "g", "v", "dh", "z", "zh", "jh"}}


def f0_mid(x, s, e):
    vals = []
    for t in numpy.arange(0.5 * (s + e) - 0.015, 0.5 * (s + e) + 0.016, 0.005):
        f, _ = pitch_at(x, int((t - 0.015) * 16000))
        if f > 0:
            vals.append(f)
    return numpy.median(vals) if len(vals) >= 3 else 0.0


def phones_bdl(labdir, wavdir, count):
    for lab in sorted(glob.glob(os.path.join(labdir, "*.lab")))[:count]:
        rows = [l.split() for l in open(lab) if len(l.split()) == 3]
        prev, out = 0.0, []
        for end, _, phone in rows:
            out.append((phone, prev, float(end)))
            prev = float(end)
        yield os.path.join(wavdir, os.path.basename(lab)[:-4] + ".wav"), out


def phones_ours(scoredir, wavdir):
    for score in sorted(glob.glob(os.path.join(scoredir, "*.score"))):
        out = []
        for l in open(score):
            m = re.match(r"# ([A-Z]+)-?\d? ([\d.]+) ([\d.]+)", l)
            if m and m.group(1) != "phones":
                s = float(m.group(2))
                out.append((m.group(1).lower(), s, s + float(m.group(3))))
        yield os.path.join(wavdir, os.path.basename(score)[:-6] + ".wav"), out


def main():
    kind = sys.argv[1]
    source = phones_bdl(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 29) if kind == "bdl" \
        else phones_ours(sys.argv[2], sys.argv[3])
    dips = {k: [] for k in CLASSES}
    for wav, phones in source:
        if not os.path.exists(wav):
            continue
        x = load(wav)
        for a, c, b in zip(phones, phones[1:], phones[2:]):
            if a[0] not in VOWELS or b[0] not in VOWELS:
                continue
            for name, members in CLASSES.items():
                if c[0] in members:
                    fa, fc, fb = f0_mid(x, *a[1:]), f0_mid(x, *c[1:]), f0_mid(x, *b[1:])
                    if fa > 0 and fc > 0 and fb > 0:
                        dips[name].append(1200 * numpy.log2(fc / numpy.sqrt(fa * fb)))
    row = []
    for name, v in dips.items():
        v = numpy.array(v)
        if len(v):
            row.append("%s n%d med %+.0f |%.0f| p90 %.0f" % (name, len(v), numpy.median(v), numpy.median(numpy.abs(v)),
                                                            numpy.percentile(numpy.abs(v), 90)))
    print("%-10s %s" % (sys.argv[3].rstrip("/").split("/")[-1], "   ".join(row)))


if __name__ == "__main__":
    main()
