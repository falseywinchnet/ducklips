"""Wideband spectrogram of ours between t0 and t1 with score phone labels; optional second file panel.
usage: labspec.py OUT WAV SCORE t0 t1 [WAV2 t0 t1]"""
import sys, re, numpy, librosa, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
out, wav, sc, t0, t1 = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5])
ph = [(m.group(1), float(m.group(2))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(sc)) if m]
panels = [(wav, t0, t1, ph)]
if len(sys.argv) > 6: panels.append((sys.argv[6], float(sys.argv[7]), float(sys.argv[8]), []))
fig, ax = plt.subplots(1, len(panels), figsize=(9 * len(panels), 5), squeeze=False)
for a, (w, a0, a1, labs) in zip(ax[0], panels):
    y, sr = librosa.load(w, sr=16000); seg = y[int(a0 * sr):int(a1 * sr)]
    S = librosa.amplitude_to_db(numpy.abs(librosa.stft(seg, n_fft=512, hop_length=32, win_length=96)), ref=numpy.max)
    a.imshow(S, origin="lower", aspect="auto", extent=[a0, a1, 0, 8000], vmin=-75, cmap="magma")
    for n, s in labs:
        if a0 < s < a1: a.axvline(s, color="c", lw=0.5); a.text(s + 0.002, 7600, n.rstrip("-0123456789"), color="c")
    a.set_title(w.split("/")[-1])
plt.tight_layout(); plt.savefig(out, dpi=65)
