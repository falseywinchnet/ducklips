"""Average narrowband spectrum of the 5 loudest 60 ms frames of each file,
dB re the peak below 1.5 kHz; prints levels of the peaks and floors in
bands and plots the slices.
usage: slice.py OUT.png LABEL WAV [LABEL WAV ...]"""
import sys, numpy, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import resample_poly


def load(p):
    r, x = wavfile.read(p); x = x.astype(float)
    if x.ndim > 1: x = x[:, 0]
    return resample_poly(x, 16000, r)


out = sys.argv[1]; pairs = list(zip(sys.argv[2::2], sys.argv[3::2]))
fig, ax = plt.subplots(len(pairs), 1, figsize=(13, 2.2 * len(pairs)), sharex=True, squeeze=False)
for k, (lab, p) in enumerate(pairs):
    x = load(p); n = 960
    e = numpy.array([numpy.sqrt(numpy.mean(x[i:i + n] ** 2)) for i in range(0, len(x) - n, 160)])
    acc = 0
    for o in numpy.argsort(e)[::-1][:5]:
        acc = acc + numpy.abs(numpy.fft.rfft(x[o * 160:o * 160 + n] * numpy.hanning(n), 4096)) ** 2
    f = numpy.fft.rfftfreq(4096, 1 / 16000); db = 10 * numpy.log10(acc + 1e-9); db -= db[f < 1500].max()
    bands = [(1500, 2500), (2500, 4000), (4000, 6000), (6000, 8000)]
    print("%-14s %s" % (lab, "  ".join("%d-%d: peak %5.1f floor %5.1f" % (a, b, db[(f >= a) & (f < b)].max(), numpy.percentile(db[(f >= a) & (f < b)], 10)) for a, b in bands)))
    ax[k, 0].plot(f, db, lw=0.7); ax[k, 0].set_ylim(-80, 5); ax[k, 0].set_title(lab, fontsize=9); ax[k, 0].grid(alpha=0.3)
plt.tight_layout(); plt.savefig(out, dpi=70)
