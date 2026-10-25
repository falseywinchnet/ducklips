"""Equivalence gate: the in-process engine (analysis/engine.py) against the
frozen reference (the reference tree's posture_recover.render driving its
own build_dev/sustain), byte for byte, over a matrix of renders; and the
snapshot gate: render(0..N) == render(0..k) + resume(k..N), bit-identical,
also with inputs changed from step k on.

usage (on the Mini): python3 bindings/equivalence.py [--engine-only] [--cli BINARY]
  --cli BINARY: also run this tree's sustain binary through the reference
  render (its SUSTAIN replaced) and require the same bytes.

The reference module is imported read-only (no bytecode written) and its
temporary folders are redirected from /tmp/es to /tmp/esf/ref."""
import os
import sys

sys.dont_write_bytecode = True

import numpy

HERE = os.path.dirname(os.path.abspath(__file__))
FAST = os.path.dirname(HERE)
REFERENCE = os.path.expanduser("~/ducklips/tract")
SCRATCH = "/tmp/esf/equivalence"

sys.path.insert(0, os.path.join(FAST, "analysis"))
import engine  # noqa: E402

sys.path.insert(0, os.path.join(REFERENCE, "analysis"))
import posture_recover as reference  # noqa: E402
import tempfile as _tempfile  # noqa: E402


class _RedirectedTempfile:
    """The reference's tempfile, its folders in /tmp/esf/ref (never /tmp/es)."""

    @staticmethod
    def TemporaryDirectory(dir=None):  # noqa: N802 (the module's name)
        os.makedirs("/tmp/esf/ref", exist_ok=True)
        return _tempfile.TemporaryDirectory(dir="/tmp/esf/ref")


reference.tempfile = _RedirectedTempfile
assert reference.SUSTAIN == os.path.join(REFERENCE, "build_dev", "sustain"), reference.SUSTAIN
for name in ("NODE", "VISCOUS", "HOP", "JET_CLOSING", "NOSE_POLE", "NOSE_ZERO"):
    assert getattr(reference, name) == getattr(engine, name), name

STEP_RATE = 96000.0
NODE = engine.NODE


def story(name, count=46):
    """A Story (2008) vowel resampled onto `count` nodes."""
    with open(os.path.join(REFERENCE, "data", "story_2008_areas.txt")) as handle:
        for line in handle:
            parts = line.split()
            if parts and parts[0] == name:
                section = float(parts[1])
                areas = numpy.array([float(v) for v in parts[2:]])
                break
    x = (numpy.arange(count) + 0.5) * NODE
    centres = (numpy.arange(len(areas)) + 0.5) * section
    return numpy.interp(x, centres, areas)


AA = story("aa")
IY = story("iy")
UW = story("uw")
LENGTH = 46.6 * NODE


def closure(tube, node, area):
    out = tube.copy()
    out[node] = area
    out[max(node - 1, 0)] = max(area * 3.0, out[max(node - 1, 0)] * 0.5)
    return out


def word_keyframes(first, middle, last, hold=0.07, move=0.05):
    times = [0.0, hold, hold + move, 2 * hold + move, 2 * hold + 2 * move]
    frames = [first, first, middle, middle, last]
    return times, frames


def forced_waves(seconds, nan_stretches):
    steps = int(seconds * STEP_RATE) + 5
    t = numpy.arange(steps) / STEP_RATE
    wave = 2500.0 * numpy.sin(2 * numpy.pi * 118.0 * t) + 800.0 * numpy.sin(2 * numpy.pi * 410.0 * t)
    for begin, end in nan_stretches:
        wave[int(begin * STEP_RATE):int(end * STEP_RATE)] = numpy.nan
    return wave


def cases():
    """(name, arguments, files the render writes) for the matrix."""
    out = []
    out.append(("plain /aa/ 0.3 s", dict(tube_and_length=(AA, LENGTH)), ()))
    out.append(("plain /iy/ odd pitch and seconds (%g rounding)",
                dict(tube_and_length=(IY, LENGTH), seconds=0.4123456, pitch=117.123456789), ()))
    times, frames = word_keyframes(AA, closure(AA, 44, 1e-4), AA)
    out.append(("/aa b aa/ closure with seal 0.41 s",
                dict(tube_and_length=(AA, LENGTH), seconds=0.41, keyframes=(times, frames), seal=0.01), ()))
    velum = ([0.0, 0.1, 0.14, 0.2, 0.26, 0.3], [0.0, 0.0, 0.4, 0.4, 0.0, 0.0])
    times, frames = word_keyframes(AA, closure(AA, 44, 1e-4), AA)
    out.append(("nasal /aa m aa/ velum, port waves, nose dump 0.54 s",
                dict(tube_and_length=(AA, LENGTH), seconds=0.54, keyframes=(times, frames), seal=0.01, velum=velum,
                     nose_dump="NOSEDUMP", port_waves="PORTWAVES"), ("NOSEDUMP", "PORTWAVES")))
    times, frames = word_keyframes(AA, closure(AA, 40, 0.05), IY)
    gate = ([0.0, 0.1, 0.12, 0.18, 0.2], [0.0, 0.0, 1.0, 1.0, 0.0])
    out.append(("jet /aa s iy/ with gate and closing (seal) 0.41 s",
                dict(tube_and_length=(AA, LENGTH), seconds=0.41, keyframes=(times, frames), seal=0.01, jet=0.4,
                     jet_gate=gate), ()))
    times, frames = word_keyframes(AA, closure(AA, 40, 0.05), AA)
    out.append(("jet without gate, no seal 0.3 s",
                dict(tube_and_length=(AA, LENGTH), seconds=0.3, keyframes=(times, frames), jet=0.4), ()))
    folds = ([0.0, 0.09, 0.11, 0.19, 0.21],
             [(0.00495, 0.0004, 0.0002), (0.00495, 0.0004, 0.0002), (0.0501, 0.0501, 0.0002),
              (0.0501, 0.0501, 0.0002), (0.00495, 0.0004, 0.0002)])
    times, frames = word_keyframes(AA, closure(AA, 44, 1e-4), AA)
    out.append(("voiceless /aa p aa/ fold keys (FOLDKEYS) 0.41 s",
                dict(tube_and_length=(AA, LENGTH), seconds=0.41, keyframes=(times, frames), seal=0.01, folds=folds),
                ()))
    triggers = ([0.09, 0.2], [(0.0501, 0.0501, 0.0002, 0.02), (0.00495, 0.0004, 0.0002, 0.02)])
    out.append(("fold triggers (FOLDSEQ) 0.3 s",
                dict(tube_and_length=(AA, LENGTH), seconds=0.3, keyframes=(times, frames), seal=0.01,
                     folds=triggers), ()))
    out.append(("forced waves with NaN stretches 0.3 s",
                dict(tube_and_length=(UW, LENGTH), seconds=0.3,
                     forced=forced_waves(0.3, [(0.05, 0.08), (0.2, 0.31)])), ()))
    observed = numpy.full(int(0.3 * STEP_RATE) + 5, numpy.nan)
    observed[: int(0.15 * STEP_RATE)] = 150.0 + 100.0 * numpy.sin(
        2 * numpy.pi * 118.0 * numpy.arange(int(0.15 * STEP_RATE)) / STEP_RATE)
    out.append(("observer flows on forced folds 0.3 s",
                dict(tube_and_length=(AA, LENGTH), seconds=0.3, forced=forced_waves(0.3, [(0.1, 0.25)]),
                     observed=(observed, 0.5)), ()))
    out.append(("glide to /iy/ (AREAFILE2) with wave file 0.3 s",
                dict(tube_and_length=(AA, LENGTH), seconds=0.3, second=IY, glide=(0.1, 0.1), wave="WAVE"),
                ("WAVE",)))
    times, frames = word_keyframes(AA, closure(AA, 44, 1e-4), AA)
    out.append(("everything: nasal, jet, folds, seal 0.54 s",
                dict(tube_and_length=(AA, LENGTH), seconds=0.54, keyframes=(times, frames), seal=0.01, velum=velum,
                     jet=0.4, jet_gate=gate, folds=folds, port_waves="PORTWAVES"), ("PORTWAVES",)))
    return out


def same_bytes(a, b):
    return numpy.asarray(a).tobytes() == numpy.asarray(b).tobytes()


def with_paths(arguments, files, folder):
    out = dict(arguments)
    for name in files:
        key = {"NOSEDUMP": "nose_dump", "PORTWAVES": "port_waves", "WAVE": "wave"}[name]
        out[key] = os.path.join(folder, name)
    return out


def compare(name, arguments, files, cli=None):
    ok = True
    os.makedirs(SCRATCH, exist_ok=True)
    folders = {"reference": os.path.join(SCRATCH, "reference"), "engine": os.path.join(SCRATCH, "engine"),
               "cli": os.path.join(SCRATCH, "cli")}
    for folder in folders.values():
        os.makedirs(folder, exist_ok=True)
    reference.SUSTAIN = os.path.join(REFERENCE, "build_dev", "sustain")
    glottal, radiated = reference.render(**with_paths(arguments, files, folders["reference"]))
    results = {"engine": engine.render(**with_paths(arguments, files, folders["engine"]))}
    if cli is not None:
        reference.SUSTAIN = cli
        results["cli"] = reference.render(**with_paths(arguments, files, folders["cli"]))
        reference.SUSTAIN = os.path.join(REFERENCE, "build_dev", "sustain")
    line = "%-60s %6d steps" % (name, len(glottal))
    for which, (g, r) in results.items():
        good = len(g) == len(glottal) and same_bytes(g, glottal) and same_bytes(r, radiated)
        for file_name in files:
            with open(os.path.join(folders["reference"], file_name), "rb") as handle:
                expected = handle.read()
            with open(os.path.join(folders[which], file_name), "rb") as handle:
                good = good and handle.read() == expected
        if not good and len(g) == len(glottal):
            deviation = max(numpy.max(numpy.abs(g - glottal)), numpy.max(numpy.abs(r - radiated)))
            line += "  %s DIFFERS (max |d| %.3g)" % (which, deviation)
        else:
            line += "  %s %s" % (which, "identical" if good else "DIFFERS")
        ok = ok and good
    for folder in folders.values():
        for file_name in os.listdir(folder):
            os.remove(os.path.join(folder, file_name))
    print(line)
    sys.stdout.flush()
    return ok


def snapshot_gate():
    """render(0..N) == render(0..k) + resume(k..N); then with the forced
    waves and a late keyframe changed from step k on."""
    ok = True
    times, frames = word_keyframes(AA, closure(AA, 44, 1e-4), AA)
    velum = ([0.0, 0.1, 0.14, 0.2, 0.26, 0.3], [0.0, 0.0, 0.4, 0.4, 0.0, 0.0])
    gate = ([0.0, 0.1, 0.12, 0.18, 0.2], [0.0, 0.0, 1.0, 1.0, 0.0])
    base = dict(tube_and_length=(AA, LENGTH), seconds=0.41, keyframes=(times, frames), seal=0.01, velum=velum,
                jet=0.4, jet_gate=gate, forced=forced_waves(0.41, [(0.12, 0.42)]))
    settings = engine.prepare(**base)
    whole, waves, finite = engine.Run(settings).advance(port_waves=True)
    for k in (1, 12, 96, 5000, 19213, 30000):
        run = engine.Run(settings)
        head, head_waves, _ = run.advance(k, port_waves=True)
        snap = run.snapshot()
        tail, tail_waves, _ = snap.advance(port_waves=True)
        again, again_waves, _ = run.advance(port_waves=True)
        good = (same_bytes(numpy.concatenate([head, tail]), whole) and same_bytes(numpy.concatenate([head, again]), whole)
                and same_bytes(numpy.concatenate([head_waves, tail_waves]), waves))
        print("snapshot at step %6d: resume %s" % (k, "identical" if good else "DIFFERS"))
        ok = ok and good
    # Inputs changed from step k on: forced waves past k, and a keyframe
    # whose influence starts after k (the last frame, moved after k's tick).
    k = 18000
    changed = dict(base)
    changed["forced"] = base["forced"].copy()
    changed["forced"][k:] = 1500.0 * numpy.sin(2 * numpy.pi * 97.0 * numpy.arange(len(changed["forced"]) - k) / STEP_RATE)
    later_times = list(times)
    later_frames = list(frames)
    later_frames[-1] = IY
    changed["keyframes"] = (later_times, later_frames)
    # The last frame is reached from times[3] on: its change acts from the
    # tick after times[3] (step k must not precede that tick).
    assert times[3] * STEP_RATE >= k
    other = engine.prepare(**changed)
    truth = engine.Run(other).advance()[0]
    run = engine.Run(settings)
    head = run.advance(k)[0]
    tail = run.snapshot().advance(settings=other)[0]
    good = same_bytes(numpy.concatenate([head, tail]), truth)
    print("snapshot at step %6d, inputs changed from it: %s" % (k, "identical" if good else "DIFFERS"))
    return ok and good


def incremental_gate():
    """engine.Incremental over a sequence of changed renders (as a search
    makes them) against the reference, byte for byte, and how much of each
    it reused."""
    ok = True
    times, frames = word_keyframes(AA, closure(AA, 44, 1e-4), AA)
    velum = ([0.0, 0.1, 0.14, 0.2, 0.26, 0.3], [0.0, 0.0, 0.4, 0.4, 0.0, 0.0])
    gate = ([0.0, 0.1, 0.12, 0.18, 0.2], [0.0, 0.0, 1.0, 1.0, 0.0])
    folds = ([0.0, 0.09, 0.11, 0.19, 0.21],
             [(0.00495, 0.0004, 0.0002), (0.00495, 0.0004, 0.0002), (0.0501, 0.0501, 0.0002),
              (0.0501, 0.0501, 0.0002), (0.00495, 0.0004, 0.0002)])
    base = dict(tube_and_length=(AA, LENGTH), seconds=0.41, keyframes=(times, frames), seal=0.01, jet=0.4,
                jet_gate=gate, folds=folds)
    sequence = [("first render", base), ("the same again", base)]
    later = list(frames)
    later[4] = IY
    sequence.append(("last frame /iy/", dict(base, keyframes=(times, later))))
    moved = list(times)
    moved[3] = times[3] + 0.0013
    sequence.append(("fourth keyframe 1.3 ms later", dict(base, keyframes=(moved, later))))
    sequence.append(("longer (0.45 s)", dict(base, keyframes=(moved, later), seconds=0.45)))
    sequence.append(("jet gate's last key moved", dict(base, keyframes=(moved, later), seconds=0.45,
                                                      jet_gate=(gate[0][:-1] + [0.23], gate[1]))))
    forced = forced_waves(0.45, [(0.1, 0.46)])
    sequence.append(("forced to 0.1 s", dict(base, keyframes=(moved, later), seconds=0.45, forced=forced)))
    forced_more = forced.copy()
    forced_more[int(0.08 * STEP_RATE):] = numpy.nan
    sequence.append(("forced to 0.08 s", dict(base, keyframes=(moved, later), seconds=0.45, forced=forced_more)))
    sequence.append(("back to the first", base))
    nasal = dict(tube_and_length=(AA, LENGTH), seconds=0.54, keyframes=(times, frames), seal=0.01, velum=velum)
    sequence.append(("nasal word", nasal))
    sequence.append(("nasal, velum's last key later", dict(nasal, velum=(velum[0][:-1] + [0.33], velum[1]))))
    sequence.append(("nasal, initial tube changed", dict(nasal, tube_and_length=(UW, LENGTH))))
    incremental = engine.Incremental()
    for name, arguments in sequence:
        got = incremental.render(**arguments)
        expected = reference.render(**arguments)
        good = same_bytes(got[0], expected[0]) and same_bytes(got[1], expected[1])
        print("incremental %-34s reused %6d of %6d steps: %s" % (name, incremental.reused, len(expected[0]),
                                                                "identical" if good else "DIFFERS"))
        sys.stdout.flush()
        ok = ok and good
    return ok


def main():
    for name in engine.SUSTAIN_VARIABLES:
        assert name not in os.environ, "clear %s from the environment for the gate" % name
    cli = None
    if "--cli" in sys.argv:
        cli = os.path.abspath(sys.argv[sys.argv.index("--cli") + 1])
    ok = True
    if "--quick" in sys.argv:
        quick = cases()
        for name, arguments, files in (quick[0], quick[4]):
            ok = compare(name, arguments, files, cli) and ok
        print("QUICK %s" % ("PASSED" if ok else "FAILED"))
        return 0 if ok else 1
    if "--engine-only" not in sys.argv:
        for name, arguments, files in cases():
            ok = compare(name, arguments, files, cli) and ok
        # A caller's environment passes through to the tool: both must see it.
        os.environ["TREMOR"] = "0.01"
        os.environ["PITCHEND"] = "131.5"
        ok = compare("environment passed through (TREMOR, PITCHEND)", dict(tube_and_length=(AA, LENGTH)), (), cli) and ok
        del os.environ["TREMOR"]
        del os.environ["PITCHEND"]
    ok = snapshot_gate() and ok
    ok = incremental_gate() and ok
    print("EQUIVALENCE %s" % ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
