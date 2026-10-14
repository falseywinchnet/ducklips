"""Wobble inside long vowels: F0 and level every 5 ms (30 ms window,
normalized autocorrelation with parabolic refinement) through each vowel of
at least 120 ms, 25 ms trimmed at each end; a quadratic fit removes the
intonation glide, and the residual's RMS is the wobble (cents; dB). Medians
over vowels.
usage: vowel_wobble.py bdl LABDIR WAVDIR [N] | vowel_wobble.py ours SCOREDIR WAVDIR"""
import sys, os, re, glob, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly

RATE = 16000
VOWELS = {"aa", "ae", "ah", "ao", "aw", "ay", "eh", "er", "ey", "ih", "iy", "ow", "oy", "uh", "uw"}


def load(path):
    r, x = wavfile.read(path)
    x = x.astype(float)
    if x.ndim > 1:
        x = x[:, 0]
    return resample_poly(x, RATE, r) if r != RATE else x


def pitch_at(x, i, win=480):
    lo, hi = int(RATE / 300), int(RATE / 60)
    a = x[i:i + win] - x[i:i + win].mean()
    ea = numpy.dot(a, a)
    if ea <= 0 or i + hi + 1 + win > len(x):
        return 0.0, 0.0
    cs = numpy.empty(hi + 2)
    for k in range(lo - 1, hi + 2):
        b = x[i + k:i + k + win] - x[i + k:i + k + win].mean()
        cs[k] = numpy.dot(a, b) / numpy.sqrt(ea * numpy.dot(b, b) + 1e-30)
    k = lo + int(numpy.argmax(cs[lo:hi]))
    if cs[k] < 0.6:
        return 0.0, 10 * numpy.log10(ea / win + 1e-9)
    den = cs[k - 1] - 2 * cs[k] + cs[k + 1]
    shift = 0.5 * (cs[k - 1] - cs[k + 1]) / den if den < 0 else 0.0
    return RATE / (k + shift), 10 * numpy.log10(ea / win + 1e-9)


def wobble(x, start, end):
    f0, lev = [], []
    for t in numpy.arange(start + 0.025, end - 0.025 - 0.03, 0.005):
        f, l = pitch_at(x, int(t * RATE))
        f0.append(f)
        lev.append(l)
    f0, lev = numpy.array(f0), numpy.array(lev)
    if len(f0) < 12 or numpy.mean(f0 > 0) < 0.9:
        return None
    keep = f0 > 0
    t = numpy.arange(len(f0))[keep]
    c = 1200 * numpy.log2(f0[keep] / numpy.median(f0[keep]))
    if numpy.max(numpy.abs(numpy.diff(c))) > 250:
        return None
    rc = c - numpy.polyval(numpy.polyfit(t, c, 2), t)
    rl = lev[keep] - numpy.polyval(numpy.polyfit(t, lev[keep], 2), t)
    return numpy.sqrt(numpy.mean(rc ** 2)), numpy.sqrt(numpy.mean(rl ** 2))


def vowels_bdl(labdir, wavdir, count):
    for lab in sorted(glob.glob(os.path.join(labdir, "*.lab")))[:count]:
        rows = [l.split() for l in open(lab) if len(l.split()) == 3]
        wav = os.path.join(wavdir, os.path.basename(lab)[:-4] + ".wav")
        if not os.path.exists(wav):
            continue
        prev = 0.0
        out = []
        for end, _, phone in rows:
            end = float(end)
            if phone in VOWELS and end - prev >= 0.12:
                out.append((prev, end))
            prev = end
        yield wav, out


def vowels_ours(scoredir, wavdir):
    for score in sorted(glob.glob(os.path.join(scoredir, "*.score"))):
        wav = os.path.join(wavdir, os.path.basename(score)[:-6] + ".wav")
        if not os.path.exists(wav):
            continue
        out = []
        for l in open(score):
            m = re.match(r"# ([A-Z]+)\d ([\d.]+) ([\d.]+) ([\d.]+)", l)
            if m and m.group(1).lower() in VOWELS and float(m.group(3)) >= 0.12:
                s = float(m.group(2))
                out.append((s, s + float(m.group(3))))
        yield wav, out


def main():
    kind = sys.argv[1]
    source = vowels_bdl(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 29) if kind == "bdl" \
        else vowels_ours(sys.argv[2], sys.argv[3])
    pw, lw = [], []
    for wav, spans in source:
        x = load(wav)
        for s, e in spans:
            w = wobble(x, s, e)
            if w:
                pw.append(w[0])
                lw.append(w[1])
    print("%-30s vowels %3d  pitch wobble %5.1f c  level wobble %.2f dB" % (
        sys.argv[3].rstrip("/").split("/")[-1], len(pw), numpy.median(pw), numpy.median(lw)))


if __name__ == "__main__":
    main()
