"""Per-render timing: the reference (posture_recover.render: the reference
tool as a subprocess, its keyframe files, its RAW file) against the
in-process engine (analysis/engine.py), for a 0.41 s and a 0.54 s word,
with and without the nose and the jet.

usage (on the Mini): python3 bindings/bench.py [REPEATS] [--library PATH]...
Reports per render the median wall time and the median CPU time (user +
sys, the subprocess's included for the reference) over REPEATS (default 9),
interleaved so the machine's load falls on both alike. Each --library adds
another build of the engine to compare (e.g. the unoptimized one)."""
import ctypes
import os
import resource
import sys
import time

sys.dont_write_bytecode = True

import numpy

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import equivalence  # noqa: E402  (the reference module, redirected; the cases' tubes)

engine = equivalence.engine
reference = equivalence.reference


def cpu_now():
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def words():
    aa = equivalence.AA
    length = equivalence.LENGTH
    times, frames = equivalence.word_keyframes(aa, equivalence.closure(aa, 44, 1e-4), aa)
    velum = ([0.0, 0.1, 0.14, 0.2, 0.26, 0.3], [0.0, 0.0, 0.4, 0.4, 0.0, 0.0])
    fric_times, fric_frames = equivalence.word_keyframes(aa, equivalence.closure(aa, 40, 0.05), equivalence.IY)
    gate = ([0.0, 0.1, 0.12, 0.18, 0.2], [0.0, 0.0, 1.0, 1.0, 0.0])
    out = []
    for seconds in (0.41, 0.54):
        out.append(("%.2f s stop word (seal)" % seconds,
                    dict(tube_and_length=(aa, length), seconds=seconds, keyframes=(times, frames), seal=0.01)))
        out.append(("%.2f s nasal word (nose)" % seconds,
                    dict(tube_and_length=(aa, length), seconds=seconds, keyframes=(times, frames), seal=0.01,
                         velum=velum)))
        out.append(("%.2f s fricative word (jet)" % seconds,
                    dict(tube_and_length=(aa, length), seconds=seconds, keyframes=(fric_times, fric_frames),
                         seal=0.01, jet=0.4, jet_gate=gate)))
    return out


def timed(function, arguments):
    cpu = cpu_now()
    wall = time.perf_counter()
    result = function(**arguments)
    return time.perf_counter() - wall, cpu_now() - cpu, result


def main():
    repeats = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 9
    libraries = [sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == "--library"]
    engines = [("engine", engine.LIBRARY)] + [("lib%d" % i, path) for i, path in enumerate(libraries)]
    print("%-28s %-10s %9s %9s %8s" % ("word", "path", "wall ms", "cpu ms", "speedup"))
    for name, arguments in words():
        samples = {"reference": []}
        for label, _ in engines:
            samples[label] = []
        expected = None
        for _ in range(repeats):
            wall, cpu, result = timed(reference.render, arguments)
            samples["reference"].append((wall, cpu))
            expected = result
            for label, path in engines:
                engine._lib = None
                engine.LIBRARY = path
                wall, cpu, result = timed(engine.render, arguments)
                samples[label].append((wall, cpu))
                assert equivalence.same_bytes(result[0], expected[0]) and equivalence.same_bytes(result[1], expected[1]), label
        base_wall = numpy.median([s[0] for s in samples["reference"]])
        for label in samples:
            wall = numpy.median([s[0] for s in samples[label]])
            cpu = numpy.median([s[1] for s in samples[label]])
            print("%-28s %-10s %9.1f %9.1f %7.1fx" % (name, label, 1e3 * wall, 1e3 * cpu, base_wall / wall))
        sys.stdout.flush()
    engine._lib = None
    engine.LIBRARY = engines[0][1]


if __name__ == "__main__":
    main()
