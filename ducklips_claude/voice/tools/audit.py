"""Intent against outcome: did the body do what the plan asked?

A recording only shows what came out. Our engine also knows what was
meant: the planner says, phone by phone, what must be voiced, at what
pitch, and where the tract must close. This tool renders text with
say --audit (the body's realized state at 4 kHz: glottal area and flow,
pressures, the narrowest oral section, lips, velum, output level) and
checks every planned phone against it. Failures of the body to follow the
plan are engine defects, located in time and attributed to a cause:

  NO VOICE     a voiced phone whose folds hardly vibrate (dropout)
  BREAK        a gap in vibration inside a voiced stretch (stutter)
  ROUGH        irregular cycles in a vowel (jitter above 4%)
  PITCH        a vowel more than 2 semitones off its planned pitch
  VOICED       a voiceless phone that vibrates through its core
  NO CLOSURE   a stop or nasal whose tract never shuts
  DIP          a vowel whose level sags in its middle

usage: audit.py OUT_PREFIX "text" ["text" ...] [-- say arguments]
writes OUT_PREFIX_N.png (timeline per text) and prints the findings."""
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca

ROOT = ca.ROOT
RATE = 4000.0
VOWELS = set(ca.VOWELS)
VOICED = VOWELS | {"b", "d", "g", "v", "dh", "z", "zh", "jh", "m", "n", "ng", "l", "r", "w", "y", "dx"}
VOICELESS = {"p", "t", "k", "f", "th", "s", "sh", "ch", "hh"}
CLOSING = {"p", "t", "k", "b", "d", "g", "m", "n", "ng", "ch", "jh"}
VIBRATION = 0.02   # cm^2 peak-to-peak glottal area over 12 ms that counts as vibrating
BREAK_MS = 10.0    # shortest gap in vibration reported as a break


def render(text, prefix, extra):
    wave = prefix + ".wav"
    score = prefix + ".score"
    audit = prefix + ".audit"
    subprocess.run([os.path.join(ROOT, "build", "say"), text, wave, "--score", score, "--audit", audit] + extra,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True, cwd=ROOT)
    phones = []
    for line in open(score):
        parts = line.split()
        if line.startswith("# ") and len(parts) == 5 and parts[1] != "phones:":
            name = parts[1].rstrip("-0123456789").lower()
            start, duration = float(parts[2]), float(parts[3])
            phones.append({"name": name, "start": start, "end": start + duration, "pitch": float(parts[4])})
    words = []
    for line in open(score):
        parts = line.split()
        if len(parts) >= 5 and parts[1] == "mark" and parts[2] == "word":
            words.append((float(parts[0]), parts[4]))
    data = numpy.loadtxt(audit)
    return phones, words, data


def smooth(values, seconds):
    width = max(int(seconds * RATE), 1)
    return numpy.convolve(values, numpy.ones(width) / width, "same")


def vibration_amplitude(area):
    """Peak-to-peak oscillation of the glottal area over a sliding 12 ms
    window, around its own 8 ms moving average (so the folds moving apart
    or together for a consonant is not counted as vibration)."""
    area = area - smooth(area, 0.008)
    half = int(0.006 * RATE)
    padded = numpy.pad(area, half, mode="edge")
    windows = numpy.lib.stride_tricks.sliding_window_view(padded, 2 * half + 1)
    return windows.max(axis=1) - windows.min(axis=1)


def cycles(area, vibrating):
    """Opening instants (s) of the glottis inside vibrating stretches."""
    times = []
    for index in range(1, len(area) - 1):
        if not vibrating[index]:
            continue
        local = area[max(index - 48, 0):index + 48]
        threshold = local.min() + 0.4 * (local.max() - local.min())
        if area[index - 1] < threshold <= area[index]:
            times.append(index / RATE)
    return numpy.array(times)


def word_of(words, time):
    current = "?"
    for start, word in words:
        if start <= time + 1e-6:
            current = word
    return current


def audit_text(text, prefix, extra):
    phones, words, data = render(text, prefix, extra)
    time = data[:, 0]
    area = data[:, 1]
    level = numpy.sqrt(smooth(data[:, 9] ** 2, 0.010))  # 10 ms level
    amplitude = vibration_amplitude(area)
    vibrating = amplitude > VIBRATION
    openings = cycles(area, vibrating)
    findings = []

    def frames(start, end, core=0.8):
        margin = 0.5 * (1.0 - core) * (end - start)
        return (time >= start + margin) & (time < end - margin)

    for phone in phones:
        name, start, end = phone["name"], phone["start"], phone["end"]
        where = "%-10s %6.3f s  /%s/" % (word_of(words, start), start, name)
        core = frames(start, end)
        if not core.any():
            continue
        voiced_share = vibrating[core].mean()
        if name in VOICED and voiced_share < 0.5:
            findings.append((3.0 - voiced_share, "NO VOICE   %s  vibrating %.0f%% of its core" % (where, 100 * voiced_share)))
        if name in VOICELESS and voiced_share > 0.8 and end - start > 0.04:
            findings.append((1.0, "VOICED     %s  vibrating %.0f%% of its core" % (where, 100 * voiced_share)))
        if name in CLOSING:
            least = data[frames(start, end, 1.0), 4].min()
            if least > 0.02:
                findings.append((2.0, "NO CLOSURE %s  narrowest %.3f cm^2" % (where, least)))
        if name in VOWELS:
            inside = openings[(openings >= start + 0.25 * (end - start)) & (openings < end - 0.25 * (end - start))]
            if len(inside) >= 4:
                periods = numpy.diff(inside)
                jitter = 100.0 * numpy.mean(numpy.abs(numpy.diff(periods))) / periods.mean()
                realized = 1.0 / numpy.median(periods)
                if jitter > 4.0:
                    findings.append((1.0 + jitter / 10.0, "ROUGH      %s  jitter %.0f%%" % (where, jitter)))
                if phone["pitch"] > 0.0:
                    error = 12.0 * numpy.log2(realized / phone["pitch"])
                    if abs(error) > 2.0:
                        findings.append((abs(error) / 4.0, "PITCH      %s  %.0f Hz for %.0f planned (%+.1f st)" % (
                            where, realized, phone["pitch"], error)))
            middle = level[frames(start, end, 0.6)]
            if len(middle) > 8 and middle.max() > 0.0:
                sag = 20.0 * numpy.log10(middle.min() / middle.max() + 1e-9)
                if sag < -12.0:
                    findings.append((-sag / 10.0, "DIP        %s  level sags %.0f dB mid-vowel" % (where, sag)))
    # Breaks: gaps in vibration inside runs of voiced phones. A run ends at
    # a voiceless phone or a pause, and the last 80 ms before a pause (the
    # planned run-down of the phrase) is not searched.
    phrase_end = {}
    for index, phone in enumerate(phones):
        following = phones[index + 1]["start"] if index + 1 < len(phones) else 1e9
        phrase_end[index] = following - phone["end"] > 0.03
    runs = []
    current = []
    for index, phone in enumerate(phones):
        if phone["name"] in VOICED:
            current.append(index)
        elif current:
            runs.append(current)
            current = []
        if phrase_end[index] and current:
            runs.append(current)
            current = []
    if current:
        runs.append(current)
    for run in runs:
        first, last = phones[run[0]], phones[run[-1]]
        tail = 0.08 if phrase_end[run[-1]] else 0.02
        span = (time >= first["start"] + 0.02) & (time < last["end"] - tail)
        mask = vibrating[span]
        stamps = time[span]
        gap = 0
        for k, value in enumerate(mask):
            if not value:
                gap += 1
                continue
            if gap / RATE * 1000.0 >= BREAK_MS:
                at = stamps[k] - gap / RATE
                inside = [p["name"] for p in phones if p["start"] <= at < p["end"]]
                findings.append((2.5 + gap / RATE * 20.0, "BREAK      %-10s %6.3f s  /%s/  no vibration for %.0f ms" % (
                    word_of(words, at), at, inside[0] if inside else "?", gap / RATE * 1000.0)))
            gap = 0
    draw(prefix + ".png", text, phones, data, amplitude, vibrating, openings)
    return findings


def draw(path, text, phones, data, amplitude, vibrating, openings):
    time = data[:, 0]
    first = max(phones[0]["start"] - 0.05, 0.0)
    last = phones[-1]["end"] + 0.1
    keep = (time >= first) & (time <= last)
    figure, axes = plt.subplots(4, 1, figsize=(15, 9), sharex=True)
    axes[0].plot(time[keep], 20 * numpy.log10(data[keep, 9] + 1e-6), lw=0.6)
    axes[0].set_ylabel("output dB")
    axes[0].set_title(text, fontsize=9)
    voiced_plan = numpy.zeros(len(time))
    for phone in phones:
        if phone["name"] in VOICED:
            voiced_plan[(time >= phone["start"]) & (time < phone["end"])] = 1.0
    axes[1].fill_between(time[keep], 0, voiced_plan[keep] * amplitude[keep].max(), color="0.85", label="planned voiced")
    axes[1].plot(time[keep], amplitude[keep], lw=0.7, label="glottal swing (cm^2)")
    axes[1].axhline(VIBRATION, color="r", lw=0.5)
    axes[1].legend(fontsize=7, loc="upper right")
    if len(openings) > 2:
        mids = 0.5 * (openings[1:] + openings[:-1])
        f0 = 1.0 / numpy.diff(openings)
        ok = f0 < 400
        axes[2].plot(mids[ok], f0[ok], ".", ms=2, label="realized F0")
    for phone in phones:
        if phone["pitch"] > 0:
            axes[2].plot([phone["start"], phone["end"]], [phone["pitch"]] * 2, "r", lw=1.5)
    axes[2].set_ylabel("Hz")
    axes[2].set_ylim(60, 250)
    axes[3].semilogy(time[keep], numpy.maximum(data[keep, 4], 1e-4), lw=0.7, label="narrowest oral area")
    axes[3].set_ylabel("cm^2")
    twin = axes[3].twinx()
    twin.plot(time[keep], data[keep, 5], color="tab:orange", lw=0.5)
    twin.set_ylabel("place cm", color="tab:orange")
    for axis in axes:
        for phone in phones:
            axis.axvline(phone["start"], color="0.7", lw=0.4)
    for phone in phones:
        axes[3].text(0.5 * (phone["start"] + phone["end"]), 1.5e-4, phone["name"], ha="center", fontsize=7)
    figure.tight_layout()
    figure.savefig(path, dpi=70)
    plt.close(figure)


def main():
    arguments = sys.argv[1:]
    extra = []
    if "--" in arguments:
        extra = arguments[arguments.index("--") + 1:]
        arguments = arguments[:arguments.index("--")]
    prefix = arguments[0]
    total = 0
    for number, text in enumerate(arguments[1:], 1):
        findings = audit_text(text, "%s_%d" % (prefix, number), extra)
        total += len(findings)
        print("== %s  (%d findings)" % (text, len(findings)))
        for severity, message in sorted(findings, reverse=True):
            print("  " + message)
    print("TOTAL %d" % total)


if __name__ == "__main__":
    main()
