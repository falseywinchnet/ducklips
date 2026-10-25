"""Fidelity of a voice to a reference reading of the same text, in dB (the
owner, 2026-10-18: "set a goal for yourself and seek to work carefully until
you get to 50 db ... use technical measures only ... use profiles").

Both readings are cut into phones (the reference by its labels, or mapped
from ours by dynamic time warping on MFCCs; ours by its score), and each
phone class (every vowel; the consonants by manner and place) is described
by five groups of features, averaged over its tokens:
  spectrum - the mean level in 24 ERB-spaced bands (100-7600 Hz) over the
             middle 60% of each token, dB re the utterance's speech level;
  dynamics - how far each of six octave bands (125-8000 Hz) moves within a
             token: the standard deviation of its 10 ms level (dB), the
             wavelet-scale behaviour of the class;
  phase    - the harmonics' phase scatter from cycle to cycle in 0-1, 1-2,
             2-4 kHz (np.angle phase distortion's deviation, rad,
             cord_features.phase_distortion), pooled over front, central and
             back vowels (per vowel it was noise: bdl against himself 4 dB);
  pitch    - the reading's F0 at its 5th, 25th, 50th, 75th and 95th
             percentiles, semitones re 100 Hz (per vowel the medians barely
             spread and intonation swamped them: bdl against himself -6 dB);
  timing   - the mean duration, natural log of seconds.
Per group: SNR = 10 log10(spread of the reference's features about their
mean / our squared error against them), over the features both readings
have; for phase, an absolute magnitude, the signal is the features' own
sum of squares. The combined score weighs each group alike: 10 log10(number of groups
/ sum of each group's error over its spread). 0 dB: no better than the
reference's mean; 50 dB: an error power 1e-5 of the spread.
usage: fidelity.py arctic REFDIR IDS OURSDIR      (ours: OURSDIR/ID.wav + ID.score)
       fidelity.py arctic2 REFDIR_A IDS_A REFDIR_B IDS_B  (reader against reader, labels both)
       fidelity.py nato OURS.wav OURS.score      (the NATO reader's clips, mapped by DTW)
IDS: first-last, e.g. a0001-a0030."""
import difflib
import os
import sys

import librosa
import numpy
import soundfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cord_features

RATE = 16000
HOP = 160          # 10 ms
FRAME = 400        # 25 ms
VOWELS = {"aa", "ae", "ah", "ao", "aw", "ay", "eh", "er", "ey", "ih", "iy", "ow", "oy", "uh", "uw"}
VOWEL_GROUP = {"iy": "front", "ih": "front", "ey": "front", "eh": "front", "ae": "front", "ah": "central",
               "er": "central", "aw": "central", "ay": "central", "aa": "back", "ao": "back", "ow": "back",
               "oy": "back", "uh": "back", "uw": "back"}
CLASS = {"ax": "ah", "axr": "er", "p": "stop_voiceless", "t": "stop_voiceless", "k": "stop_voiceless",
         "b": "stop_voiced", "d": "stop_voiced", "g": "stop_voiced", "dx": "stop_voiced",
         "f": "fric_labial", "v": "fric_labial_voiced", "th": "fric_dental", "dh": "fric_dental_voiced",
         "s": "fric_alveolar", "z": "fric_alveolar_voiced", "sh": "fric_postalveolar", "zh": "fric_postalveolar",
         "hh": "fric_glottal", "ch": "affricate", "jh": "affricate_voiced", "m": "nasal_labial", "n": "nasal_alveolar",
         "ng": "nasal_velar", "l": "liquid_lateral", "r": "liquid_rhotic", "w": "glide_labial", "y": "glide_palatal"}
GROUPS = ["spectrum", "dynamics", "phase", "pitch", "timing"]
# Phase pooled over all vowels unless FIDELITY_PHASE_GROUPS=1 (front, central,
# back). Pooled, each feature has three times the tokens: on 20 sentences the
# back vowels' scatter swung 0.41-0.93 rad between render seeds of one voice.
PHASE_GROUPS = os.environ.get("FIDELITY_PHASE_GROUPS", "0") == "1"


def phone_class(name):
    base = name.lower().rstrip("012").rstrip("-")
    if base in VOWELS:
        return base
    return CLASS.get(base)


def erb_bands(count=24, low=100.0, high=7600.0):
    def erb(f):
        return 21.4 * numpy.log10(1 + 0.00437 * f)

    def inverse(e):
        return (10 ** (e / 21.4) - 1) / 0.00437
    edges = inverse(numpy.linspace(erb(low), erb(high), count + 2))
    hz = numpy.fft.rfftfreq(512, 1.0 / RATE)
    weights = numpy.zeros((count, len(hz)))
    for b in range(count):
        lo, mid, hi = edges[b], edges[b + 1], edges[b + 2]
        rise = (hz - lo) / (mid - lo)
        fall = (hi - hz) / (hi - mid)
        weights[b] = numpy.clip(numpy.minimum(rise, fall), 0.0, None)
    return weights


ERB = erb_bands()
OCTAVES = [(125, 250), (250, 500), (500, 1000), (1000, 2000), (2000, 4000), (4000, 8000)]


def frame_levels(x):
    """10 ms frames: ERB band levels (dB) and octave band levels (dB)."""
    window = numpy.hanning(FRAME)
    count = 1 + max(len(x) - FRAME, 0) // HOP
    hz = numpy.fft.rfftfreq(512, 1.0 / RATE)
    erb = numpy.zeros((count, ERB.shape[0]))
    octave = numpy.zeros((count, len(OCTAVES)))
    total = numpy.zeros(count)
    for i in range(count):
        piece = x[i * HOP:i * HOP + FRAME]
        if len(piece) < FRAME:
            piece = numpy.pad(piece, (0, FRAME - len(piece)))
        power = numpy.abs(numpy.fft.rfft(piece * window, 512)) ** 2 + 1e-12
        erb[i] = 10 * numpy.log10(ERB @ power + 1e-12)
        for k, (lo, hi) in enumerate(OCTAVES):
            octave[i, k] = 10 * numpy.log10(power[(hz >= lo) & (hz < hi)].sum() + 1e-12)
        total[i] = 10 * numpy.log10(power.sum())
    return erb, octave, total


def describe(utterances):
    """utterances: list of (x at RATE, [(phone, start, end), ...]). Returns {group: {feature: value}}."""
    pools = {}
    pitch_frames = []
    for x, segments in utterances:
        erb, octave, total = frame_levels(x)
        speech = total > total.max() - 40.0
        level = numpy.mean(total[speech]) if speech.any() else 0.0
        f0, voiced, _ = librosa.pyin(x, fmin=60, fmax=400, sr=RATE, frame_length=1024, hop_length=HOP)
        n = min(len(f0), len(total))
        good = numpy.isfinite(f0[:n]) & (total[:n] > total.max() - 30.0)
        pitch_frames.extend(list(f0[:n][good]))
        for phone, start, end in segments:
            c = phone_class(phone)
            if c is None or end - start < 0.03:
                continue
            a = int((start + 0.2 * (end - start)) * RATE / HOP)
            b = max(int((end - 0.2 * (end - start)) * RATE / HOP), a + 1)
            a0 = int(start * RATE / HOP)
            b0 = max(int(end * RATE / HOP), a0 + 1)
            if b > len(erb):
                continue
            pool = pools.setdefault(c, {"spectrum": [], "dynamics": [], "phase": [], "pitch": [], "timing": []})
            pool["spectrum"].append(erb[a:b].mean(axis=0) - level)
            if b0 - a0 >= 3:
                pool["dynamics"].append(octave[a0:b0].std(axis=0))
            pool["timing"].append(numpy.log(end - start))
            if c in VOWELS:
                f = f0[a:b][numpy.isfinite(f0[a:b])] if b <= len(f0) else numpy.array([])
                if len(f) >= 2 and end - start >= 0.06 and b0 <= len(f0):
                    # The token at its own pitch frame by frame (at its median,
                    # any glide within the vowel counted as phase scatter).
                    piece = x[a0 * HOP:b0 * HOP]
                    pd = cord_features.phase_distortion(piece, RATE, f0[a0:b0], HOP)
                    values = [pd["pdd_0_1k"], pd["pdd_1_2k"], pd["pdd_2_4k"]]
                    if all(numpy.isfinite(values)):
                        group_pool = pools.setdefault("vowels_" + (VOWEL_GROUP[c] if PHASE_GROUPS else "all"),
                                                      {g: [] for g in GROUPS})
                        group_pool["phase"].append(numpy.array(values))
    out = {group: {} for group in GROUPS}
    if len(pitch_frames) >= 20:
        semitones = 12 * numpy.log2(numpy.array(pitch_frames) / 100.0)
        for q in (5, 25, 50, 75, 95):
            out["pitch"]["f0_p%02d" % q] = float(numpy.percentile(semitones, q))
    for c, pool in pools.items():
        for group in GROUPS:
            values = pool[group]
            if not values:
                continue
            # Phase by the median over tokens: a few badly tracked tokens swung
            # its mean (and with it the combined score) by 5-8 dB on 10
            # sentences.
            if group == "phase":
                mean = numpy.median(numpy.array(values), axis=0)
            else:
                mean = numpy.mean(numpy.array(values), axis=0)
            mean = numpy.atleast_1d(mean)
            for k, v in enumerate(mean):
                out[group]["%s_%d" % (c, k)] = float(v)
    return out


def score(reference, ours):
    """Per-group SNR (dB) and the combined score."""
    result = {}
    normalized_error = 0.0
    groups = 0
    for group in GROUPS:
        keys = sorted(set(reference[group]) & set(ours[group]))
        if len(keys) < 3:
            continue
        r = numpy.array([reference[group][k] for k in keys])
        o = numpy.array([ours[group][k] for k in keys])
        # Phase scatter is an absolute magnitude (rad), not a level with an
        # arbitrary offset: its signal is its own size. (As spread about its
        # mean, bdl against himself on two 20-sentence sets scored -1.3 dB:
        # nine features at 0.3-0.8 rad have almost no reliable spread.)
        spread = numpy.sum(r ** 2) if group == "phase" else numpy.sum((r - r.mean()) ** 2)
        error = numpy.sum((o - r) ** 2)
        result[group] = 10 * numpy.log10(spread / max(error, 1e-12))
        normalized_error += error / max(spread, 1e-12)
        groups += 1
    result["combined"] = 10 * numpy.log10(groups / max(normalized_error, 1e-12)) if groups else float("nan")
    return result


def read_audio(path):
    x, r = soundfile.read(path)
    x = x if x.ndim == 1 else x.mean(axis=1)
    return librosa.resample(x, orig_sr=r, target_sr=RATE) if r != RATE else x


def arctic_labels(path):
    segments = []
    previous = 0.0
    started = False
    for line in open(path):
        parts = line.split()
        if line.startswith("#"):
            started = True
            continue
        if not started or len(parts) < 3:
            continue
        end = float(parts[0])
        segments.append((parts[2], previous, end))
        previous = end
    return segments


def score_segments(path):
    segments = []
    for line in open(path):
        parts = line.split()
        if len(parts) == 5 and parts[0] == "#" and parts[1] != "phones:":
            start = float(parts[2])
            segments.append((parts[1], start, start + float(parts[3])))
    return segments


def id_range(text):
    """"a0001-a0030", or with a step "a0001-a0593/20" (every 20th: a spread
    across the recording sessions; bdl's F0 on a0201-a0300 sits 0.7-1.6 st
    above a0101-a0120)."""
    span, _, step = text.partition("/")
    first, last = span.split("-")
    return ["%s%04d" % (first[0], n) for n in range(int(first[1:]), int(last[1:]) + 1, int(step or 1))]


def arctic_set(ref_dir, ids):
    utterances = []
    for name in ids:
        wav = os.path.join(ref_dir, "wav", "arctic_%s.wav" % name)
        lab = os.path.join(ref_dir, "lab", "arctic_%s.lab" % name)
        if os.path.exists(wav) and os.path.exists(lab):
            utterances.append((read_audio(wav), arctic_labels(lab)))
    return utterances


SILENT_LABELS = {"pau", "sil", "h#", "ssil"}


def base_name(phone):
    return phone.lower().rstrip("012").rstrip("-")


def aligned_names(ours, reference):
    """Our phone segments renamed to the reference's labels where the two
    sequences align one to one (difflib on the names; a replaced run of equal
    length maps position by position). So a class holds the same words in
    both readings: our planner's weak forms are AX where ARCTIC's labels keep
    the citation vowel (16 of 319 phones in bdl a0101-a0110: "that", "and",
    "would", "the"), which put 14 /ae/ in bdl's class against our 6."""
    spoken = [(p, a, b) for p, a, b in reference if p.lower() not in SILENT_LABELS]
    matcher = difflib.SequenceMatcher(a=[p.lower() for p, _, _ in spoken], b=[base_name(p) for p, _, _ in ours],
                                      autojunk=False)
    renamed = list(ours)
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "replace" and i2 - i1 == j2 - j1:
            for i, j in zip(range(i1, i2), range(j1, j2)):
                renamed[j] = (spoken[i][0].lower(), ours[j][1], ours[j][2])
    return renamed


# Our boundaries: "planned" (the score's) or "mapped" (the reference's labels
# carried onto our audio by DTW on MFCCs, as its labels were placed on its
# own: by the sound). Planned, our gestures led the phone boundaries, so a
# token's inside was already settled: faster moves lowered the within-token
# spread (dynamics 4.5 -> 3.7 at half the move times) where bdl's labels sit
# on his acoustic changes.
BOUNDARIES = os.environ.get("FIDELITY_BOUNDARIES", "planned")


def ours_set(ours_dir, ids, ref_dir=None):
    """Our renders; with `ref_dir`, each one's phones carry the reference
    sentence's labels where they align (aligned_names), or with BOUNDARIES
    "mapped" its labels outright, placed on our audio by DTW."""
    utterances = []
    for name in ids:
        wav = os.path.join(ours_dir, name + ".wav")
        sc = os.path.join(ours_dir, name + ".score")
        if os.path.exists(wav) and os.path.exists(sc):
            segments = score_segments(sc)
            x = read_audio(wav)
            lab = os.path.join(ref_dir, "lab", "arctic_%s.lab" % name) if ref_dir else None
            if lab and os.path.exists(lab):
                if BOUNDARIES == "mapped":
                    reference = read_audio(os.path.join(ref_dir, "wav", "arctic_%s.wav" % name))
                    spoken = [(p, a, b) for p, a, b in arctic_labels(lab) if p.lower() not in SILENT_LABELS]
                    segments = mapped_segments(x, reference, spoken)
                else:
                    segments = aligned_names(segments, arctic_labels(lab))
            utterances.append((x, segments))
    return utterances


def mfcc(x):
    return librosa.feature.mfcc(y=x, sr=RATE, n_mfcc=13, n_fft=512, hop_length=HOP, win_length=FRAME)[1:]


def mapped_segments(reference, ours, ours_segments):
    """Our phone boundaries carried onto the reference by DTW on MFCCs."""
    _, path = librosa.sequence.dtw(X=mfcc(ours), Y=mfcc(reference), metric="euclidean")
    path = path[::-1]
    ours_to_ref = {}
    for i, j in path:
        ours_to_ref.setdefault(i, j)

    def carry(t):
        i = int(round(t * RATE / HOP))
        keys = numpy.array(sorted(ours_to_ref))
        i = keys[numpy.argmin(numpy.abs(keys - i))]
        return ours_to_ref[i] * HOP / RATE
    return [(phone, carry(a), carry(b)) for phone, a, b in ours_segments]


NATO_WORDS = "alfa bravo charlie delta echo foxtrot golf hotel india juliett kilo lima mike november oscar papa quebec romeo sierra tango uniform victor whiskey xray yankee zulu".split()


def nato_sets(ours_wav, ours_score, ref_dir):
    """Our NATO render split into words by pauses; each word against the reader's clip."""
    x = read_audio(ours_wav)
    segments = score_segments(ours_score)
    words = []
    current = []
    end = None
    for phone, a, b in segments:
        if current and a - end > 0.1:
            words.append(current)
            current = []
        current.append((phone, a, b))
        end = b
    words.append(current)
    ours = []
    reference = []
    for word, segs in zip(NATO_WORDS, words):
        a = max(segs[0][1] - 0.03, 0.0)
        b = segs[-1][2] + 0.05
        piece = x[int(a * RATE):int(b * RATE)]
        shifted = [(p, s - a, e - a) for p, s, e in segs]
        ref = read_audio(os.path.join(ref_dir, "nato_%s.wav" % word[0]))
        level = numpy.convolve(ref ** 2, numpy.ones(HOP) / HOP, "same")
        on = numpy.nonzero(level > level.max() * 10 ** (-3.5))[0]
        ref = ref[max(on[0] - 480, 0):on[-1] + 480]
        ours.append((piece, shifted))
        reference.append((ref, mapped_segments(ref, piece, shifted)))
    return reference, ours


def report(result):
    return " ".join("%s %.1f" % (k, result[k]) for k in GROUPS + ["combined"] if k in result)


def main():
    mode = sys.argv[1]
    if mode == "arctic":
        ids = id_range(sys.argv[3])
        reference = describe(arctic_set(sys.argv[2], ids))
        ours = describe(ours_set(sys.argv[4], ids))
        print(report(score(reference, ours)))
    elif mode == "arctic2":
        reference = describe(arctic_set(sys.argv[2], id_range(sys.argv[3])))
        other = describe(arctic_set(sys.argv[4], id_range(sys.argv[5])))
        print(report(score(reference, other)))
    elif mode == "nato":
        reference, ours = nato_sets(sys.argv[2], sys.argv[3], os.path.expanduser("~/ducklips/voice/ref/nato"))
        print(report(score(describe(reference), describe(ours))))


if __name__ == "__main__":
    main()
