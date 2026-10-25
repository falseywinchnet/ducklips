"""Vocal-fold features of held vowels, for fitting the glottis to human
recordings (owner, 2026-10-18: "finetune ... against our human reference
recordings for pure vocal cord sounds ... taking the angle on the specgram
... align the harmonic behavior").

Per stretch (like vowels compared, so the tract largely factors out):
  tilt      - the harmonic peaks' level against log2 frequency, 0.3-5 kHz
              (dB per octave): the spectrogram's angle across frequency;
  h1h2      - first harmonic minus second (dB);
  hnr_<band> - harmonic peaks over the median trough between them, per band
              (0-1, 1-2, 2-4, 4-6 kHz; dB): how far up the harmonics stand;
  drift, tremor, flutter - the pitch track (pyin, 5 ms) in cents re its mean,
              less its straight-line trend, RMS in 0.5-3, 3-8 and 8-20 Hz:
              the harmonics' waving in time, which grows with harmonic
              number (the spectrogram's angle in time);
  wobble    - the level track (dB, 5 ms) likewise, RMS in 3-20 Hz;
  pdd_<band>, pdm_<band> - the phase (owner: "taking np.angle ie the
              phase"): every 5 ms, over three periods (Hann), each harmonic's
              np.angle at k f0; its phase distortion PD_k = phi_(k+1) -
              phi_k - phi_1 (Degottex), which no longer depends on where the
              frame falls; per harmonic over time, the deviation
              sqrt(-2 ln R) (R the length of the mean of exp(i PD)): how far
              the harmonics' relative phases scatter from cycle to cycle
              (rad); and the mean's angle (rad), the pulse's shape. Medians
              over the harmonics in each band.
usage: cord_features.py owner WAV          the owner's voiced stretches, vowel by formants
       cord_features.py arctic LAB WAV N   bdl's long aa / iy / uw
       cord_features.py ours WAV:START:END:VOWEL ..."""
import glob
import os
import sys

import librosa
import numpy
import soundfile
from scipy.signal import lfilter

FEATURES = ["tilt", "h1h2", "hnr_0_1k", "hnr_1_2k", "hnr_2_4k", "hnr_4_6k", "drift", "tremor", "flutter", "wobble",
            "pdd_0_1k", "pdd_1_2k", "pdd_2_4k", "pdm_0_1k", "pdm_1_2k", "pdm_2_4k"]
BANDS = [(0, 1000), (1000, 2000), (2000, 4000), (4000, 6000)]
VOWEL_OF = {"aa": "ah", "iy": "ee", "uw": "oo"}


def band_rms(track, rate, low, high):
    if len(track) < 8:
        return float("nan")
    spectrum = numpy.fft.rfft(track * numpy.hanning(len(track)))
    hz = numpy.fft.rfftfreq(len(track), 1.0 / rate)
    keep = (hz >= low) & (hz < high)
    # Parseval with the Hann window's power (3/8).
    power = 2.0 * numpy.sum(numpy.abs(spectrum[keep]) ** 2) / (len(track) ** 2 * 0.375)
    return float(numpy.sqrt(power))


def detrended(track):
    t = numpy.arange(len(track))
    return track - numpy.polyval(numpy.polyfit(t, track, 1), t)


def features(x, r):
    if x.ndim > 1:
        x = x.mean(axis=1)
    hop = int(0.005 * r)
    f0, voiced, _ = librosa.pyin(x, fmin=60, fmax=320, sr=r, frame_length=2048, hop_length=hop)
    good = numpy.isfinite(f0)
    if good.sum() < 20:
        return None
    pitch = float(numpy.median(f0[good]))
    cents = 1200.0 * numpy.log2(numpy.interp(numpy.arange(len(f0)), numpy.nonzero(good)[0], f0[good]) / pitch)
    cents = detrended(cents)
    track_rate = r / hop
    level = 20.0 * numpy.log10(librosa.feature.rms(y=x, frame_length=4 * hop, hop_length=hop)[0] + 1e-9)
    level = detrended(level)
    window = numpy.hanning(len(x))
    n = 1 << int(numpy.ceil(numpy.log2(len(x) * 4)))
    spectrum = 20.0 * numpy.log10(numpy.abs(numpy.fft.rfft(x * window, n)) + 1e-9)
    hz = numpy.fft.rfftfreq(n, 1.0 / r)
    peak_f, peak_db, trough_db = [], [], []
    k = 1
    while (k + 0.5) * pitch < min(6000.0, r / 2 - pitch):
        around = (hz > (k - 0.25) * pitch) & (hz < (k + 0.25) * pitch)
        between = (hz > (k + 0.3) * pitch) & (hz < (k + 0.7) * pitch)
        peak_f.append(k * pitch)
        peak_db.append(spectrum[around].max())
        trough_db.append(numpy.median(spectrum[between]))
        k += 1
    peak_f = numpy.array(peak_f)
    peak_db = numpy.array(peak_db)
    trough_db = numpy.array(trough_db)
    span = (peak_f >= 300) & (peak_f <= 5000)
    out = {"pitch": pitch,
           "tilt": float(numpy.polyfit(numpy.log2(peak_f[span]), peak_db[span], 1)[0]),
           "h1h2": float(peak_db[0] - peak_db[1])}
    for low, high in BANDS:
        pick = (peak_f >= max(low, 1.5 * pitch)) & (peak_f < high)
        out["hnr_%d_%dk" % (low // 1000, high // 1000)] = float(numpy.median(peak_db[pick] - trough_db[pick])) if pick.any() else float("nan")
    out["drift"] = band_rms(cents, track_rate, 0.5, 3.0)
    out["tremor"] = band_rms(cents, track_rate, 3.0, 8.0)
    out["flutter"] = band_rms(cents, track_rate, 8.0, 20.0)
    out["wobble"] = band_rms(level, track_rate, 3.0, 20.0)
    out.update(phase_distortion(x, r, f0, hop))
    return out


def phase_distortion(x, r, f0, hop):
    centres = []
    rows = []
    top = 4000.0
    for frame, pitch in enumerate(f0):
        if not numpy.isfinite(pitch):
            continue
        half = int(1.5 * r / pitch)
        centre = frame * hop
        if centre - half < 0 or centre + half >= len(x):
            continue
        piece = x[centre - half:centre + half] * numpy.hanning(2 * half)
        count = int(top / pitch)
        k = numpy.arange(1, count + 2)
        n = numpy.arange(-half, half)
        phases = numpy.angle(numpy.exp(-2j * numpy.pi * numpy.outer(k * pitch / r, n)) @ piece)
        pd = numpy.angle(numpy.exp(1j * (phases[1:] - phases[:-1] - phases[0])))
        rows.append((pitch, pd))
        centres.append(centre)
    out = {}
    if len(rows) < 10:
        return {name: float("nan") for name in ("pdd_0_1k", "pdd_1_2k", "pdd_2_4k", "pdm_0_1k", "pdm_1_2k", "pdm_2_4k")}
    count = min(len(pd) for _, pd in rows)
    stack = numpy.array([pd[:count] for _, pd in rows])
    pitch = numpy.median([p for p, _ in rows])
    mean_vector = numpy.mean(numpy.exp(1j * stack), axis=0)
    deviation = numpy.sqrt(-2.0 * numpy.log(numpy.clip(numpy.abs(mean_vector), 1e-6, 1.0)))
    angle = numpy.angle(mean_vector)
    freq = (numpy.arange(count) + 1) * pitch
    for low, high, name in ((0, 1000, "0_1k"), (1000, 2000, "1_2k"), (2000, 4000, "2_4k")):
        pick = (freq >= low) & (freq < high)
        out["pdd_" + name] = float(numpy.median(deviation[pick])) if pick.any() else float("nan")
        out["pdm_" + name] = float(numpy.angle(numpy.mean(numpy.exp(1j * angle[pick])))) if pick.any() else float("nan")
    return out


def formant_vowel(x, r):
    # LPC formants on a 16 kHz copy: ah has F1 over 600; ee F2 over 1800; oo the rest.
    y = librosa.resample(x, orig_sr=r, target_sr=16000)
    y = lfilter([1.0, -0.97], [1.0], y)
    a = librosa.lpc(y * numpy.hamming(len(y)), order=14)
    roots = [z for z in numpy.roots(a) if numpy.imag(z) > 0]
    freqs = sorted(numpy.angle(roots) * 16000 / (2 * numpy.pi))
    freqs = [f for f in freqs if f > 200]
    f1, f2 = freqs[0], freqs[1]
    vowel = "ah" if f1 > 600 else ("ee" if f2 > 1800 else "oo")
    return vowel, f1, f2


def owner(path):
    x, r = soundfile.read(path)
    if x.ndim > 1:
        x = x.mean(axis=1)
    hop = 512
    f0, voiced, _ = librosa.pyin(x, fmin=60, fmax=320, sr=r, frame_length=2048, hop_length=hop)
    rows = []
    start = None
    for index, on in enumerate(list(voiced) + [False]):
        if on and start is None:
            start = index
        elif not on and start is not None:
            if (index - start) * hop / r > 0.3:
                a = start * hop + int(0.05 * r)
                b = index * hop - int(0.05 * r)
                vowel, f1, f2 = formant_vowel(x[a:b], r)
                found = features(x[a:b], r)
                if found:
                    found.update(vowel=vowel, where="%.2f-%.2f" % (a / r, b / r), f1=f1, f2=f2)
                    rows.append(found)
            start = None
    return rows


def arctic(lab_dir, wav_dir, limit):
    rows = []
    for lab in sorted(glob.glob(os.path.join(lab_dir, "*.lab")))[:limit]:
        wav = os.path.join(wav_dir, os.path.basename(lab).replace(".lab", ".wav"))
        if not os.path.exists(wav):
            continue
        x, r = soundfile.read(wav)
        previous = 0.0
        for line in open(lab):
            parts = line.split()
            if len(parts) != 3 or parts[0] == "#":
                continue
            end = float(parts[0])
            phone = parts[2].lower()
            if phone in VOWEL_OF and end - previous > 0.15:
                trim = 0.15 * (end - previous)
                found = features(x[int((previous + trim) * r):int((end - trim) * r)], r)
                if found:
                    found.update(vowel=VOWEL_OF[phone], where=os.path.basename(wav))
                    rows.append(found)
            previous = end
    return rows


def show(label, rows):
    for vowel in ("ah", "ee", "oo"):
        pick = [row for row in rows if row["vowel"] == vowel]
        if not pick:
            continue
        line = " ".join("%s %6.2f" % (name, numpy.nanmedian([row[name] for row in pick])) for name in FEATURES)
        print("%-6s %s n %2d F0 %5.1f | %s" % (label, vowel, len(pick), numpy.median([row["pitch"] for row in pick]), line))


def main():
    mode = sys.argv[1]
    if mode == "owner":
        rows = owner(sys.argv[2])
        for row in rows:
            print("  %s %s F1 %.0f F2 %.0f" % (row["where"], row["vowel"], row["f1"], row["f2"]))
        show("owner", rows)
    elif mode == "arctic":
        show("bdl", arctic(sys.argv[2], sys.argv[3], int(sys.argv[4])))
    else:
        rows = []
        for spec in sys.argv[2:]:
            path, a, b, vowel = spec.split(":")
            x, r = soundfile.read(path)
            found = features(x[int(float(a) * r):int(float(b) * r)], r)
            if found:
                found.update(vowel=vowel)
                rows.append(found)
        show("ours", rows)


if __name__ == "__main__":
    main()
