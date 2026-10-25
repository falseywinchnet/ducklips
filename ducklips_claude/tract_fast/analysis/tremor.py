"""Tremor in a held vowel: the pitch and loudness contours (every 2.5 ms)
over START-END, slow drift (< 2.5 Hz) removed; their rms and the frequency
of their strongest component in 2.5-20 Hz; plus cycle-to-cycle jitter and
shimmer from the glottal cycles (CYCLES file: period ms in column 2,
amplitude in column 4). usage: tremor.py START END WAV[:CYCLES] ..."""
import sys, numpy, soundfile, librosa
from scipy.signal import butter, sosfiltfilt
a, b = float(sys.argv[1]), float(sys.argv[2])
for spec in sys.argv[3:]:
    path, _, cycles = spec.partition(":")
    x, r = soundfile.read(path)
    seg = x[int(a * r):int(b * r)]
    hop = int(0.0025 * r)
    f0, v, _ = librosa.pyin(seg, fmin=70, fmax=260, sr=r, frame_length=2048, hop_length=hop)
    f0 = numpy.nan_to_num(f0, nan=numpy.nanmedian(f0))
    env = librosa.feature.rms(y=seg, frame_length=1024, hop_length=hop)[0]
    rate = 1.0 / 0.0025
    sos = butter(2, 2.5, "highpass", fs=rate, output="sos")
    out = []
    for name, track, unit in (("pitch", 100 * (f0 / numpy.median(f0) - 1), "%"), ("loudness", 20 * numpy.log10(env / numpy.median(env) + 1e-9), "dB")):
        hp = sosfiltfilt(sos, track)[40:-40]
        spec_ = numpy.abs(numpy.fft.rfft(hp * numpy.hanning(len(hp)), 8192))
        fr = numpy.fft.rfftfreq(8192, 1 / rate)
        band = (fr > 2.5) & (fr < 20)
        out.append("%s tremor %.2f%s rms, strongest %.1f Hz" % (name, hp.std(), unit, fr[band][numpy.argmax(spec_[band])]))
    if cycles:
        d = numpy.loadtxt(cycles)
        d = d[(d[:, 0] > a) & (d[:, 0] < b)]
        per, amp = d[:, 1], d[:, 3]
        out.append("jitter %.2f%% shimmer %.1f%%" % (100 * numpy.mean(numpy.abs(numpy.diff(per))) / per.mean(), 100 * numpy.mean(numpy.abs(numpy.diff(amp))) / amp.mean()))
    print("%-10s" % path.split("/")[-1], " | ".join(out))
