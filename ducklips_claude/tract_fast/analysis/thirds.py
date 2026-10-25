"""Third-octave spectra (dB re the 500 Hz band) of our held vowels against
bdl's stressed tokens of the same vowel (ARCTIC labels, middle 60%).
usage: thirds.py LABDIR WAVDIR MAXFILES vowel:WAV:START:END ..."""
import os, sys, numpy, soundfile
CENTRES = [250 * 2 ** (k / 3) for k in range(0, 15)]  # 250 Hz .. 6.3 kHz
def bands(x, r):
    n = 1 << int(numpy.ceil(numpy.log2(len(x))))
    S = numpy.abs(numpy.fft.rfft(x * numpy.hanning(len(x)), n)) ** 2
    f = numpy.fft.rfftfreq(n, 1 / r)
    out = numpy.array([S[(f >= c / 2 ** (1 / 6)) & (f < c * 2 ** (1 / 6))].sum() for c in CENTRES])
    out = 10 * numpy.log10(out + 1e-20)
    return out - out[CENTRES.index(min(CENTRES, key=lambda c: abs(c - 500)))]
labdir, wavdir, maxfiles = sys.argv[1], sys.argv[2], int(sys.argv[3])
specs = [s.split(":") for s in sys.argv[4:]]
want = {v for v, _, _, _ in specs}
tokens = {v: [] for v in want}
for name in sorted(os.listdir(labdir))[:maxfiles]:
    lines = [l.split() for l in open(os.path.join(labdir, name)) if len(l.split()) == 3]
    wav = os.path.join(wavdir, name.replace(".lab", ".wav"))
    if not os.path.exists(wav): continue
    x, r = soundfile.read(wav); prev = 0.0
    for end, _, ph in lines:
        end = float(end)
        if ph in want and end - prev >= 0.12:
            mid = (prev + end) / 2; half = (end - prev) * 0.3
            tokens[ph].append(bands(x[int((mid - half) * r):int((mid + half) * r)], r))
        prev = end
print("band Hz      " + " ".join("%5d" % c for c in CENTRES))
for v, path, a, b in specs:
    x, r = soundfile.read(path)
    ours = bands(x[int(float(a) * r):int(float(b) * r)], r)
    ref = numpy.median(numpy.array(tokens[v]), axis=0)
    print("bdl %-3s n=%-3d" % (v, len(tokens[v])) + " ".join("%5.0f" % q for q in ref))
    print("%-12s" % path.split("/")[-1][:12] + " ".join("%5.0f" % q for q in ours))
    print("difference  " + " ".join("%+5.0f" % q for q in ours - ref))
