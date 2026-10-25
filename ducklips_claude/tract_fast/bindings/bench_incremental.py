"""A search's pattern of renders, timed three ways: the reference
(posture_recover.render, a subprocess per render), the in-process engine
(engine.render, every render whole) and engine.Incremental (each render
resumed from the last one's snapshot before its first changed input).

The search moves one phone of a 0.54 s word at a time (its tube's area at
the closure node and its keyframe time, by small steps, as the timing and
tube searches of sequence_recover do), for each of the word's later
keyframes in turn. Every result is checked byte for byte against the
in-process engine's whole render.

usage (on the Mini): python3 bindings/bench_incremental.py [STEPS_PER_PHONE] [--reference] [--phones 2,3,4]
Reports CPU and wall seconds per render (median) and the share of steps
reused."""
import os
import resource
import sys
import time

sys.dont_write_bytecode = True

import numpy

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import equivalence  # noqa: E402

engine = equivalence.engine
reference = equivalence.reference


def cpu_now():
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def search_renders(per_phone, phones):
    """The renders of a search over the later keyframes of /aa b aa/ with a
    nose (velum keyframes) and the folds' keys."""
    aa = equivalence.AA
    times, frames = equivalence.word_keyframes(aa, equivalence.closure(aa, 44, 1e-4), aa)
    times = list(times)
    velum = ([0.0, 0.1, 0.14, 0.2, 0.26, 0.3], [0.0, 0.0, 0.4, 0.4, 0.0, 0.0])
    base = dict(tube_and_length=(aa, equivalence.LENGTH), seconds=0.54, seal=0.01, velum=velum)
    renders = []
    for phone in phones:
        for trial in range(per_phone):
            moved = list(times)
            moved[phone] = times[phone] + 0.0002 * (trial + 1)
            changed = [frame.copy() for frame in frames]
            changed[phone][30] *= 1.0 + 0.01 * (trial + 1)
            renders.append(("phone %d" % phone, dict(base, keyframes=(moved, changed))))
    return renders


def main():
    per_phone = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 6
    phones = (2, 3, 4)
    if "--phones" in sys.argv:
        phones = tuple(int(k) for k in sys.argv[sys.argv.index("--phones") + 1].split(","))
    renders = search_renders(per_phone, phones)
    expected = []
    timing = {"engine": [], "incremental": [], "reference": []}
    for _, arguments in renders:
        cpu = cpu_now()
        wall = time.perf_counter()
        expected.append(engine.render(**arguments))
        timing["engine"].append((time.perf_counter() - wall, cpu_now() - cpu))
    incremental = engine.Incremental()
    reused = 0
    total = 0
    for (name, arguments), truth in zip(renders, expected):
        cpu = cpu_now()
        wall = time.perf_counter()
        got = incremental.render(**arguments)
        timing["incremental"].append((time.perf_counter() - wall, cpu_now() - cpu))
        assert equivalence.same_bytes(got[0], truth[0]) and equivalence.same_bytes(got[1], truth[1]), name
        reused += incremental.reused
        total += len(truth[0])
    if "--reference" in sys.argv:
        for (name, arguments), truth in zip(renders, expected):
            cpu = cpu_now()
            wall = time.perf_counter()
            got = reference.render(**arguments)
            timing["reference"].append((time.perf_counter() - wall, cpu_now() - cpu))
            assert equivalence.same_bytes(got[0], truth[0]) and equivalence.same_bytes(got[1], truth[1]), name
    print("%d renders of a 0.54 s nasal word, keyframes %s moved; steps reused by the incremental: %.0f%%"
          % (len(renders), ",".join(str(k) for k in phones), 100.0 * reused / total))
    for label, samples in timing.items():
        if samples:
            print("%-12s wall %7.1f ms  cpu %7.1f ms per render (median)" % (
                label, 1e3 * numpy.median([s[0] for s in samples]), 1e3 * numpy.median([s[1] for s in samples])))


if __name__ == "__main__":
    main()
