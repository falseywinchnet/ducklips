"""Why the pinch tube's rattle follows the pitch (owner, 2026-10-17: "no matter
what i do i cant disentangle oscillatory rattle from frequency- higher
frequency voices, no rattle").

Renders the held "ah" from a slider preset at several tension anchors (each
VARIANT is name:block.key=value,...), dumps the source's flow and opening
every step over the held vowel, and per variant prints:
  F0, jitter and shimmer of the flow pulses;
  the share of cycles whose flow peak alternates (period-2 / subharmonic);
  the flow's envelope flutter at 2-10 and 10-60 Hz (dB rms);
  contacts per cycle (walls meeting more than once a cycle = chatter);
  the wav's subharmonic level (energy at F0/2 against F0, dB).
A variant named ns... takes the default section factors.
usage: rattle_pitch.py PRESET.json VARIANT..."""
import json
import os
import subprocess
import sys

import numpy
import soundfile
from scipy.signal import butter, sosfiltfilt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import slider_server

HELD = (1.0, 2.0)  # s, inside AA1 (0.670 - 2.264)
RATE = 96000.0


def render(values, name, folder):
    base_text = open(os.path.join(folder, "base.params")).read()
    params = os.path.join(folder, name + ".params")
    open(params, "w").write(slider_server.with_values(base_text, values))
    wave = os.path.join(folder, name + ".wav")
    env = dict(os.environ, FOLDTUBE="1", PARAMS=params, WAVEDUMP="%g,%g" % HELD)
    score = os.path.join(slider_server.SCORES, "ah.score")
    done = subprocess.run([slider_server.SPEAK, score, wave, "73021"], env=env, capture_output=True, text=True)
    flow = []
    area = []
    for line in done.stdout.splitlines():
        if line.startswith("wave ") and "| cord" in line:
            cord = line.split("| cord")[1].split()
            flow.append(float(cord[0]))
            area.append(float(cord[1]))
    return numpy.array(flow), numpy.array(area), wave


# The owner's section factors as PinchSettings has them ("disregard section
# settings", 2026-10-17): the ee preset's own section factors blew the pinch
# open (flow 924 cm^3/s, no oscillation).
DEFAULT_SECTIONS = {"s1": (1.0, 2.0, 1.0), "s2": (2.0, 4.0, 0.1), "s3": (2.0, 1.0, 0.3),
                    "s4": (4.0, 0.1, 3.0), "s5": (0.1, 3.0, 0.1), "s6": (1.0, 2.0, 1.0)}


def autocorrelation_f0(flow):
    x = flow - flow.mean()
    x = x - numpy.convolve(x, numpy.ones(960) / 960.0, "same")  # off the 10 ms drift
    ac = numpy.correlate(x[:48000], x[:48000], "full")[48000 - 1:]
    low, high = int(RATE / 500), int(RATE / 50)
    lag = low + int(numpy.argmax(ac[low:high]))
    return RATE / lag, ac[lag] / ac[0]


def cycles(flow, f0):
    # Cycle starts: upward mean crossings of the flow smoothed below 1.6 F0.
    smooth = sosfiltfilt(butter(2, 1.6 * f0, "lowpass", fs=RATE, output="sos"), flow - flow.mean())
    return numpy.nonzero((smooth[:-1] < 0.0) & (smooth[1:] >= 0.0))[0]


def flutter(env_db, rate, low, high):
    band = sosfiltfilt(butter(2, [low, high], "bandpass", fs=rate, output="sos"), env_db)
    trim = int(0.05 * rate)
    return band[trim:-trim].std()


def describe(name, flow, area, wave):
    ac_f0, clarity = autocorrelation_f0(flow)
    if flow.std() < 0.1 * flow.mean():
        print("%-14s no oscillation: flow %.0f +- %.0f cm^3/s, opening %.3f cm^2" % (name, flow.mean(), flow.std(),
                                                                                    area.mean()))
        return
    starts = cycles(flow, ac_f0)
    periods = numpy.diff(starts) / RATE
    peaks = numpy.array([flow[a:b].max() for a, b in zip(starts[:-1], starts[1:])])
    f0 = 1.0 / periods.mean()
    jitter = 100.0 * numpy.mean(numpy.abs(numpy.diff(periods))) / periods.mean()
    shimmer = 100.0 * numpy.mean(numpy.abs(numpy.diff(peaks))) / peaks.mean()
    d = numpy.diff(peaks)
    alternate = 100.0 * numpy.mean(d[:-1] * d[1:] < 0) if len(d) > 2 else 0.0
    # Contacts: the opening at its floor (the narrowest section shut leaves the leak).
    floor = area.min()
    shut = area <= floor + 1e-5
    contacts = numpy.count_nonzero(shut[1:] & ~shut[:-1]) / max(len(periods), 1)
    # The flow's 1 ms peak envelope, in dB.
    hop = int(RATE * 0.001)
    frames = len(flow) // hop
    env = flow[:frames * hop].reshape(frames, hop).max(axis=1)
    env_db = 20.0 * numpy.log10(numpy.maximum(env, 1e-3))
    peak_db = 20.0 * numpy.log10(numpy.maximum(peaks, 1e-3))
    slow = peak_db.std()
    x, r = soundfile.read(wave)
    seg = x[int(HELD[0] * r):int(HELD[1] * r)]
    spectrum = numpy.abs(numpy.fft.rfft(seg * numpy.hanning(len(seg))))
    hz = numpy.fft.rfftfreq(len(seg), 1.0 / r)

    def level(at):
        band = (hz > at * 0.93) & (hz < at * 1.07)
        return 20.0 * numpy.log10(spectrum[band].max() + 1e-12)

    sub = max(level(f0 * 0.5), level(f0 * 1.5)) - level(f0)
    print("%-14s clarity %.2f F0 %6.1f  jitter %5.2f%%  shimmer %5.1f%%  alternating %4.0f%%  peak sd %4.2f dB  "
          "flutter 2-10 %4.2f 10-60 %4.2f dB  contacts/cycle %4.2f  subharmonic %6.1f dB  closed %4.0f%%" %
          (name, clarity, f0, jitter, shimmer, alternate, slow, flutter(env_db, 1000.0, 2, 10),
           flutter(env_db, 1000.0, 10, 60), contacts, sub, 100.0 * shut.mean()))


def main():
    preset = json.load(open(sys.argv[1]))
    folder = "/tmp/es/rattle_pitch"
    os.makedirs(folder, exist_ok=True)
    base = os.path.join(slider_server.PRESETS, os.path.basename(sys.argv[1]).replace(".json", ".params"))
    open(os.path.join(folder, "base.params"), "w").write(open(base).read())
    for variant in sys.argv[2:]:
        name, _, changes = variant.partition(":")
        values = dict(preset)
        if name.startswith("ns"):
            for section, factors in DEFAULT_SECTIONS.items():
                for part, factor in zip(("mass", "elasticity", "support"), factors):
                    values["pinch.%s_%s" % (section, part)] = factor
        for change in filter(None, changes.split(",")):
            key, _, value = change.partition("=")
            values[key] = float(value)
        flow, area, wave = render(values, name, folder)
        describe(name, flow, area, wave)


if __name__ == "__main__":
    main()
