"""Wideband spectrograms of one word: bdl (from its lab) beside ours (from the score's word marks).
usage: wordspec.py OUT.png WORD BDL_WAV BDL_LAB OURS_WAV OURS_SCORE"""
import sys, re, numpy, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import librosa
out, word, bw, bl, ow, osc = sys.argv[1:7]
# ours: word marks
marks = [(float(m.group(1)), m.group(2)) for m in (re.match(r"([\d.]+) mark word [\d.]+ (\S+)", l) for l in open(osc)) if m]
ends = [t for t, _ in marks[1:]] + [None]
ph = [(m.group(1), float(m.group(2)), float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(osc)) if m]
for (t, w), e in zip(marks, ends):
    if w.lower() == word.lower(): o0, o1 = t - 0.05, (e if e else t + 0.6) + 0.05
rows = [l.split() for l in open(bl) if len(l.split()) == 3]
print("bdl lab:", " ".join(r[2] + "@" + r[0] for r in rows))
b0, b1 = float(sys.argv[7]), float(sys.argv[8])
fig, ax = plt.subplots(1, 2, figsize=(16, 5))
for a, (path, t0, t1, labs) in zip(ax, [(bw, b0, b1, [(r[2], float(r[0])) for r in rows]), (ow, o0, o1, [(n, s) for n, s, d in ph])]):
    y, sr = librosa.load(path, sr=16000); seg = y[int(t0 * sr):int(t1 * sr)]
    S = librosa.amplitude_to_db(numpy.abs(librosa.stft(seg, n_fft=512, hop_length=32, win_length=96)), ref=numpy.max)
    a.imshow(S, origin="lower", aspect="auto", extent=[t0, t1, 0, 8000], vmin=-70, cmap="magma")
    for n, s in labs:
        if t0 < s < t1: a.axvline(s, color="c", lw=0.5); a.text(s, 7600, n, color="c")
    a.set_title(path.split("/")[-1])
plt.tight_layout(); plt.savefig(out, dpi=70)
