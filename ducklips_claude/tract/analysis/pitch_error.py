"""Pitch against the plan: every 15 ms in each planned (voiced) vowel, the
heard F0 against the score's pitch; prints the median |error| in cents,
the share within 100 cents, and the median within-vowel swing (cents).
usage: pitch_error.py SCORE WAV..."""
import sys, re, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly
score = sys.argv[1]
vowels = [(float(m.group(2)), float(m.group(3)), float(m.group(4))) for m in
          (re.match(r"# ([A-Z]+\d) ([\d.]+) ([\d.]+) ([\d.]+)", l) for l in open(score)) if m and float(m.group(4)) > 0]
for wav in sys.argv[2:]:
    r, x = wavfile.read(wav); x = resample_poly(x.astype(float), 16000, r); r = 16000
    errors, swings = [], []
    for s, d, plan in vowels:
        track = []
        for t in numpy.arange(s, s + d - 0.03, 0.015):
            seg = x[int(t*r):int((t+0.04)*r)]; seg = seg - seg.mean()
            c = numpy.correlate(seg, seg, "full")[len(seg)-1:]; lo, hi = int(r/300), int(r/60); k = lo + numpy.argmax(c[lo:hi])
            if c[0] > 0 and c[k] > 0.4*c[0]: track.append(1200*numpy.log2(r/k/plan))
        errors += [abs(v) for v in track]
        if len(track) >= 3: swings.append(max(track) - min(track))
    print("%-26s median |error| %4.0f cents, within 100 cents %3.0f%%, median swing %4.0f cents" % (
        wav.split("/")[-1], numpy.median(errors), 100*numpy.mean(numpy.array(errors) < 100), numpy.median(swings)))
