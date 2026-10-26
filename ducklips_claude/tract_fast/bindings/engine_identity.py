"""Byte identity of an engine build against another: every render of a
matrix written as raw bytes into OUT, one file a case, so two builds'
folders compare with `cmp` (an engine change is kept only where the default
renders do not move: NOTES_engine.md, NOTES_residual_speed.md).

The matrix: bindings/equivalence.py's cases (copied here) (plain, closures with seals,
nasal with port waves and nose dump, jets with gates, fold keys and
triggers, forced waves with NaN stretches, observed flows, a glide), and
the canonical equation error's final render of each word given (its
FORCEWAVE from the peel at the hinted start, JETBOUND/SEALBLEND as the
environment sets them).

usage (on the Mini; the build chosen by TRACT_ENGINE_LIB):
    TRACT_ENGINE_LIB=.../libtract_engine.dylib python3 bindings/engine_identity.py OUT [WORD ...]
then: for f in OUT_A/*; do cmp "$f" OUT_B/$(basename "$f"); done"""
import os
import sys
from typing import List

sys.dont_write_bytecode = True

import numpy

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis"))
import engine  # noqa: E402
import adjoint  # noqa: E402  (this tree's modules first: equivalence puts the reference's analysis on the path)
import puppet  # noqa: E402
import pyramid_word  # noqa: E402

STEP_RATE = 96000.0
NODE = engine.NODE


# ---- bindings/equivalence.py's matrix, copied (that module imports the
# frozen reference's analysis modules, whose names clash with this tree's) ----

def story(name, count=46):
    """A Story (2008) vowel resampled onto `count` nodes."""
    with open(os.path.join(os.path.dirname(HERE), "data", "story_2008_areas.txt")) as handle:
        for line in handle:
            parts = line.split()
            if parts and parts[0] == name:
                section = float(parts[1])
                values = []
                for value in parts[2:]:
                    values.append(float(value))
                areas = numpy.array(values)
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


def with_paths(arguments, files, folder):
    out = dict(arguments)
    for name in files:
        key = {"NOSEDUMP": "nose_dump", "PORTWAVES": "port_waves", "WAVE": "wave"}[name]
        out[key] = os.path.join(folder, name)
    return out



# ---- the identity gate ----



def file_name(index: int, name: str) -> str:
    """A case's file name: its index and its name's letters."""
    letters: List[str] = []
    for character in name:
        if character.isalnum():
            letters.append(character)
        else:
            letters.append("_")
    joined = "".join(letters)
    return "%02d_%s.bin" % (index, joined[:60])


def write_case(folder: str, index: int, name: str, arguments: dict, files: tuple) -> None:
    """One matrix case rendered; its pairs and any files it writes, bytes."""
    scratch = os.path.join(folder, "files")
    os.makedirs(scratch, exist_ok=True)
    glottal, radiated = engine.render(**with_paths(arguments, files, scratch))
    parts = [numpy.asarray(glottal, numpy.float64).tobytes(), numpy.asarray(radiated, numpy.float64).tobytes()]
    for written in files:
        path = os.path.join(scratch, written)
        with open(path, "rb") as handle:
            parts.append(handle.read())
        os.remove(path)
    with open(os.path.join(folder, file_name(index, name)), "wb") as handle:
        handle.write(b"".join(parts))


def write_word(folder: str, index: int, word: str) -> None:
    """The canonical equation error's final render of a word at its hinted
    start (adjoint.ResidualCache's own settings), its RAW pairs as bytes."""
    context = pyramid_word.Context(word.split(","), ["areas"])
    lags = puppet.friction_truth(context)
    durations = puppet.traverse_truth_durations(context)
    puppet.retime_truth_traverse(context, lags, durations)
    problem = puppet.TraverseProblem(context, lags, durations, free_lags=list(range(len(lags))))
    vector = problem.vector_of_midpoints(puppet.hinted_midpoints(context.word))
    cache = adjoint.ResidualCache(context, True)
    problem.cache = cache
    problem.residual(vector)
    settings = cache.renders.settings
    steps = settings.steps()
    raw = numpy.full((steps, 2), numpy.nan)
    run = engine.Run(settings)
    reached, finite = run.advance_into(steps, raw)
    with open(os.path.join(folder, file_name(index, word)), "wb") as handle:
        handle.write(numpy.int64(reached).tobytes())
        handle.write(raw.tobytes())


def main() -> None:
    folder = sys.argv[1]
    words = sys.argv[2:]
    os.makedirs(folder, exist_ok=True)
    index = 0
    for name, arguments, files in cases():
        write_case(folder, index, name, arguments, files)
        index += 1
    for word in words:
        write_word(folder, index, word)
        index += 1
    print("%d renders written to %s" % (index, folder))


if __name__ == "__main__":
    main()
