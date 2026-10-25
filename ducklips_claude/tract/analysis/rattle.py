"""Rattle in a held vowel: flutter of the voice's loudness (whole band and
2-5 kHz) at 10-60 Hz, from 1 ms RMS envelopes over START-END (dB rms of the
band-passed log envelope), with pulse jitter and shimmer (pulses.py).
usage: rattle.py START END WAV..."""
import sys, numpy, soundfile, librosa
from scipy.signal import butter, sosfiltfilt
sys.path.insert(0, "/tmp/es")
import pulses
a, b = float(sys.argv[1]), float(sys.argv[2])
for path in sys.argv[3:]:
    x, r = soundfile.read(path)
    seg = x[int(a * r):int(b * r)]
    hop = int(0.001 * r)
    def flutter(y):
        env = 20 * numpy.log10(librosa.feature.rms(y=y, frame_length=4 * hop, hop_length=hop)[0] + 1e-9)
        band = sosfiltfilt(butter(2, [10, 60], "bandpass", fs=1000, output="sos"), env)[50:-50]
        return band.std()
    high = sosfiltfilt(butter(4, [2000, 5000], "bandpass", fs=r, output="sos"), seg)
    jit, shim = pulses.cycles(seg, r)
    print("%-16s flutter 10-60 Hz: whole %.2f dB, 2-5 kHz %.2f dB | jitter %.2f%% shimmer %.1f%%" % (path.split("/")[-1], flutter(seg), flutter(high), jit, shim))
