"""1 kHz-locked content: the signal (48 kHz, high-passed at 2.5 kHz) folded
on 48-sample frames and averaged over a span; the average's rms re the
high-passed signal's rms, dB. Anything not locked to the 1 ms tick averages
away (about -10 log10(frames) for noise).
usage: fold.py WAV t0 t1 [...]"""
import sys, numpy
from scipy.io import wavfile
from scipy.signal import butter, sosfiltfilt
a = sys.argv[1:]
for i in range(0, len(a), 3):
    r, x = wavfile.read(a[i]); x = x.astype(float)
    x = sosfiltfilt(butter(6, 2500, btype="high", fs=r, output="sos"), x)
    s = x[int(float(a[i + 1]) * r):int(float(a[i + 2]) * r)]; n = len(s) // 48 * 48
    f = s[:n].reshape(-1, 48)
    avg = f.mean(axis=0)
    print("%-12s %s-%s locked %.1f dB (noise floor %.1f)" % (a[i].split("/")[-1], a[i + 1], a[i + 2],
          20 * numpy.log10(numpy.sqrt(numpy.mean(avg ** 2)) / numpy.sqrt(numpy.mean(s ** 2))), -10 * numpy.log10(len(f))))
