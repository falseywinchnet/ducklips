"""Phone-by-phone comparison of our renders with CMU ARCTIC (bdl) recordings in
four time-frequency views.

For each utterance: render the transcript with build/say, resample to the
reference rate, run build/tf_analyze on both, align the phone sequences
(ARCTIC labels against our planner's score), and compare matched phones:

  steady  - mean level-normalized log spectrum over the middle half of a phone
  onset   - the 30 ms after a phone starts (bursts, transitions into vowels)

Distances are RMS differences in dB on a shared log-frequency grid (1/24
octave, 100 Hz - 7.6 kHz), each spectrum normalized to its own mean level so
that loudness differences do not mask shape differences.

usage: compare_arctic.py [FIRST LAST] [--view PNG_FOR_UTTERANCE_INDEX]
"""
import json
import os
import re
import subprocess
import sys
import wave
from concurrent.futures import ThreadPoolExecutor

import numpy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "ref", "arctic_bdl")
OUT = os.path.join(ROOT, "out", "arctic")
VIEWS = ["stft", "cq", "causal", "wigner"]
VOWELS = {"aa", "ae", "ah", "ao", "aw", "ax", "ay", "eh", "er", "ey", "ih", "iy", "ow", "oy", "uh", "uw"}
CLASSES = {
    "stop": {"p", "b", "t", "d", "k", "g"},
    "fricative": {"f", "v", "th", "dh", "s", "z", "sh", "zh", "hh"},
    "affricate": {"ch", "jh"},
    "nasal": {"m", "n", "ng"},
    "liquid/glide": {"l", "r", "w", "y"},
}


def phone_class(name):
    if name in VOWELS:
        return "vowel"
    for label, members in CLASSES.items():
        if name in members:
            return label
    return "other"


def read_wave(path):
    handle = wave.open(path)
    rate = handle.getframerate()
    data = numpy.frombuffer(handle.readframes(handle.getnframes()), dtype="<i2").astype(float) / 32768.0
    return data.reshape(-1, handle.getnchannels()).mean(axis=1), rate


def write_wave(path, data, rate):
    handle = wave.open(path, "wb")
    handle.setnchannels(1)
    handle.setsampwidth(2)
    handle.setframerate(rate)
    handle.writeframes((numpy.clip(data, -1, 1) * 32767).astype("<i2").tobytes())


def resample(data, rate, target):
    count = int(round(len(data) * target / rate))
    spectrum = numpy.fft.rfft(data)
    keep = count // 2 + 1
    return numpy.fft.irfft(spectrum[:keep], count) * (count / len(data))


def transcripts():
    result = {}
    for line in open(os.path.join(REF, "arctic.data")):
        match = re.match(r'\( (\S+) "(.*)" \)', line.strip())
        if match:
            result[match.group(1)] = match.group(2)
    return result


def reference_phones(name):
    phones = []
    start = 0.0
    for line in open(os.path.join(REF, "lab", name + ".lab")):
        parts = line.split()
        if len(parts) != 3:
            continue
        end = float(parts[0])
        phones.append((parts[2], start, end))
        start = end
    return [p for p in phones if p[0] not in ("pau", "h#", "ssil")]


def our_phones(score):
    phones = []
    for line in open(score):
        match = re.match(r"# ([A-Z]+)(-?\d) ([\d.]+) ([\d.]+)", line)
        if match:
            start = float(match.group(3))
            phones.append((match.group(1).lower(), start, start + float(match.group(4))))
    return phones


def same_phone(a, b):
    groups = [{"ah", "ax"}, {"er", "axr"}, {"ih", "ix"}]
    return a == b or any(a in g and b in g for g in groups)


def align(ref, ours):
    """Edit-distance alignment of phone names; returns matched index pairs."""
    n, m = len(ref), len(ours)
    cost = numpy.zeros((n + 1, m + 1))
    cost[:, 0] = numpy.arange(n + 1)
    cost[0, :] = numpy.arange(m + 1)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            substitute = cost[i - 1, j - 1] + (0 if same_phone(ref[i - 1][0], ours[j - 1][0]) else 1)
            cost[i, j] = min(substitute, cost[i - 1, j] + 1, cost[i, j - 1] + 1)
    pairs = []
    i, j = n, m
    while i > 0 and j > 0:
        if cost[i, j] == cost[i - 1, j - 1] + (0 if same_phone(ref[i - 1][0], ours[j - 1][0]) else 1):
            if same_phone(ref[i - 1][0], ours[j - 1][0]):
                pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif cost[i, j] == cost[i - 1, j] + 1:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


GRID = 100.0 * 2 ** (numpy.arange(0, int(numpy.log2(7600 / 100) * 24) + 1) / 24.0)


_INTERPOLATION = {}


def interpolation_matrix(hz):
    """Linear interpolation in log frequency from a field's rows onto GRID, as
    one matrix (cached per row layout) so a whole field maps in one product."""
    key = (len(hz), float(hz[0]), float(hz[-1]))
    if key not in _INTERPOLATION:
        source = numpy.log2(hz)
        target = numpy.log2(GRID)
        matrix = numpy.zeros((len(GRID), len(hz)))
        upper = numpy.clip(numpy.searchsorted(source, target), 1, len(hz) - 1)
        lower = upper - 1
        weight = numpy.clip((target - source[lower]) / (source[upper] - source[lower]), 0.0, 1.0)
        matrix[numpy.arange(len(GRID)), lower] = 1.0 - weight
        matrix[numpy.arange(len(GRID)), upper] += weight
        _INTERPOLATION[key] = matrix
    return _INTERPOLATION[key]


def load_view(prefix, view):
    field = numpy.load("%s_%s_enhanced.npy" % (prefix, view))
    hz = numpy.load("%s_%s_enhanced_hz.npy" % (prefix, view))[0]
    seconds = numpy.load("%s_seconds.npy" % prefix)[0]
    level = 20 * numpy.log10(field + 1e-7)
    return interpolation_matrix(hz) @ level, seconds


def segment(mapped, seconds, start, end, reference_level):
    """Returns (level re the utterance's loud frames in dB, normalized shape)."""
    mask = (seconds >= start) & (seconds < end)
    if not mask.any():
        mask = numpy.abs(seconds - 0.5 * (start + end)) == numpy.abs(seconds - 0.5 * (start + end)).min()
    power = numpy.mean(10 ** (mapped[:, mask] / 10.0), axis=1)
    spectrum = 10 * numpy.log10(power + 1e-14)
    level = 10 * numpy.log10(numpy.mean(power) + 1e-14) - reference_level
    return level, spectrum - spectrum.mean()


def utterance_level(mapped):
    frame_power = numpy.mean(10 ** (mapped / 10.0), axis=0)
    return 10 * numpy.log10(numpy.percentile(frame_power, 95) + 1e-14)


def utterance_records(name, text, out, extra):
    """Render one transcript into OUT, analyse it in every view, and compare
    phone by phone with the human recording (whose analysis is cached in the
    shared ARCTIC folder)."""
    ours_wave = os.path.join(out, name + "_ours48.wav")
    score = os.path.join(out, name + ".score")
    subprocess.run([os.path.join(ROOT, "build", "say"), text, ours_wave, "--score", score] + extra,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    data, rate = read_wave(ours_wave)
    ref_wave = os.path.join(REF, "wav", name + ".wav")
    ref_data, ref_rate = read_wave(ref_wave)
    write_wave(os.path.join(out, name + "_ours.wav"), resample(data, rate, ref_rate), ref_rate)
    subprocess.run([os.path.join(ROOT, "build", "tf_analyze"), os.path.join(out, name + "_ours.wav"),
                    os.path.join(out, name + "_ours")], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    ref_prefix = os.path.join(OUT, name + "_ref")
    if not all(os.path.exists(ref_prefix + "_" + view + ".npy") for view in VIEWS):
        subprocess.run([os.path.join(ROOT, "build", "tf_analyze"), ref_wave, ref_prefix],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    ref = reference_phones(name)
    ours = our_phones(score)
    pairs = align(ref, ours)
    records = []
    low = GRID < 1000
    mid = (GRID >= 1000) & (GRID < 3000)
    high = GRID >= 3000
    for view in VIEWS:
        ref_mapped, ref_seconds = load_view(ref_prefix, view)
        our_mapped, our_seconds = load_view(os.path.join(out, name + "_ours"), view)
        ref_level = utterance_level(ref_mapped)
        our_level = utterance_level(our_mapped)
        for i, j in pairs:
            phone, rs, re_ = ref[i]
            _, os_, oe = ours[j]
            windows = (("steady", rs + 0.25 * (re_ - rs), re_ - 0.25 * (re_ - rs), os_ + 0.25 * (oe - os_), oe - 0.25 * (oe - os_)),
                       ("onset", rs, rs + 0.03, os_, os_ + 0.03),
                       ("release", re_ - 0.03, re_ + 0.01, oe - 0.03, oe + 0.01))
            for part, a0, a1, b0, b1 in windows:
                a_level, a = segment(ref_mapped, ref_seconds, a0, a1, ref_level)
                b_level, b = segment(our_mapped, our_seconds, b0, b1, our_level)
                audible = a_level > -40 and b_level > -40
                difference = b - a
                records.append(dict(utterance=name, phone=phone, cls=phone_class(phone), view=view, part=part,
                                    level=float(b_level - a_level), audible=bool(audible),
                                    rms=float(numpy.sqrt(numpy.mean(difference ** 2))) if audible else None,
                                    low=float(difference[low].mean()), mid=float(difference[mid].mean()),
                                    high=float(difference[high].mean())))
    return records


def utterance_job(arguments):
    name, text, out, extra = arguments
    return utterance_records(name, text, out, extra)


def compare(first, last, out, extra):
    """Every utterance in parallel; returns all phone records."""
    os.makedirs(out, exist_ok=True)
    texts = transcripts()
    jobs = [("arctic_a%04d" % number, texts["arctic_a%04d" % number], out, extra) for number in range(first, last + 1)]
    records = []
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as pool:
        for part in pool.map(utterance_job, jobs):
            records.extend(part)
    return records


def main():
    arguments = [a for a in sys.argv[1:] if not a.startswith("--")]
    first = int(arguments[0]) if arguments else 1
    last = int(arguments[1]) if len(arguments) > 1 else 20
    out = os.environ.get("COMPARE_OUT", OUT)
    extra = os.environ.get("SAY_ARGS", "").split()
    records = compare(first, last, out, extra)
    json.dump(records, open(os.path.join(out, "records.json"), "w"))
    if "--summary" in sys.argv:
        summary(records)
        four_view_summary(records)
    else:
        report(records)
        four_view_summary(records)


def four_view_summary(records):
    """The baseline: shape error per view (all four corners of the Gabor
    square), and their mean plus the level errors, as one total."""
    shapes = []
    line = "FOUR-VIEW shape (consonant/vowel)"
    for view in VIEWS:
        selected = [r for r in records if r["view"] == view and r["rms"] is not None]
        consonant = numpy.mean([r["rms"] for r in selected if r["cls"] in ("fricative", "stop", "affricate")])
        vowel = numpy.mean([r["rms"] for r in selected if r["cls"] in ("vowel", "nasal", "liquid/glide")])
        shapes.append(consonant + vowel)
        line += "  %s %.2f/%.2f" % (view, consonant, vowel)
    levels = [r for r in records if r["view"] == "cq"]
    fricative = numpy.mean([abs(r["level"]) for r in levels if r["cls"] == "fricative" and r["part"] == "steady"])
    stop = numpy.mean([abs(r["level"]) for r in levels if r["cls"] == "stop" and r["part"] == "release"])
    total = numpy.mean(shapes) + fricative + stop
    print(line)
    print("FOUR-VIEW total %.2f (mean shape %.2f + fricative level %.2f + stop release level %.2f)" %
          (total, numpy.mean(shapes), fricative, stop))
    return total


def summary(records):
    """One line: level and shape errors for the classes that matter most."""
    cq = [r for r in records if r["view"] == "cq"]
    def level_error(classes, parts):
        values = [abs(r["level"]) for r in cq if r["cls"] in classes and r["part"] in parts]
        return numpy.mean(values)
    def shape_error(classes):
        values = [r["rms"] for r in cq if r["cls"] in classes and r["rms"] is not None]
        return numpy.mean(values)
    fric = level_error({"fricative"}, {"steady"})
    stops = level_error({"stop"}, {"release"})
    shape = shape_error({"fricative", "stop", "affricate"})
    vowels = shape_error({"vowel", "nasal", "liquid/glide"})
    print("SUMMARY fricative|level| %.2f  stop-release|level| %.2f  consonant shape %.2f  vowel shape %.2f  total %.2f" %
          (fric, stops, shape, vowels, fric + stops + shape + vowels))


def report(records):
    print("\nours vs bdl. shape = RMS dB difference of level-normalized spectra (audible segments only);")
    print("level = ours minus ref, dB re each utterance's loud frames; tilt = ours-ref low/mid/high (cq view)")
    for part in ("steady", "onset", "release"):
        print("\n== %s ==" % part)
        print("%-14s %5s" % ("class", "n") + "".join("%8s" % v for v in VIEWS) + "   level   tilt l/m/h")
        for cls in ["vowel", "stop", "fricative", "affricate", "nasal", "liquid/glide"]:
            rows = [r for r in records if r["cls"] == cls and r["part"] == part]
            if not rows:
                continue
            cq = [r for r in rows if r["view"] == "cq"]
            line = "%-14s %5d" % (cls, len(cq))
            for view in VIEWS:
                values = [r["rms"] for r in rows if r["view"] == view and r["rms"] is not None]
                line += "%8.1f" % (numpy.mean(values) if values else float("nan"))
            line += "  %+6.1f" % numpy.median([r["level"] for r in cq])
            audible = [r for r in cq if r["audible"]]
            if audible:
                line += "   %+.0f/%+.0f/%+.0f" % (numpy.mean([r["low"] for r in audible]), numpy.mean([r["mid"] for r in audible]),
                                                numpy.mean([r["high"] for r in audible]))
            print(line)
    print("\nPer phone (cq view): shape dB steady/onset/release, level dB steady/release")
    phones = sorted({r["phone"] for r in records})
    rows = []
    for phone in phones:
        selected = [r for r in records if r["phone"] == phone and r["view"] == "cq"]
        def mean_of(part, key):
            values = [r[key] for r in selected if r["part"] == part and r[key] is not None]
            return numpy.mean(values) if values else float("nan")
        rows.append((phone, len(selected) // 3, mean_of("steady", "rms"), mean_of("onset", "rms"), mean_of("release", "rms"),
                     mean_of("steady", "level"), mean_of("release", "level")))
    rows.sort(key=lambda row: -numpy.nansum(row[2:5]))
    for row in rows:
        print("  %-4s n=%-3d shape %5.1f %5.1f %5.1f   level %+6.1f %+6.1f" % row)


if __name__ == "__main__":
    if "--report" in sys.argv:
        report(json.load(open(os.path.join(OUT, "records.json"))))
    else:
        main()
