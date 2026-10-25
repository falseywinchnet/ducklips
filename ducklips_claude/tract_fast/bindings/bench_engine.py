"""In-process CPU time per render for several builds of the engine library,
interleaved (rounds of a few renders each, alternating builds) so the
machine's load falls on all of them alike.

usage (on the Mini): python3 bindings/bench_engine.py ROUNDS LIBRARY...
Prints per word and build the minimum and median CPU ms per render (process
time) and the minimum's ratio to the first build's. (The Mini has 4
performance and 6 efficiency cores under a load of 15-20: the minimum is the
render that ran on a performance core; the median mixes the two.)"""
import os
import sys
import time

sys.dont_write_bytecode = True

import numpy

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bench  # noqa: E402

engine = bench.engine


def settings_for(library, arguments):
    engine._lib = None
    engine.LIBRARY = library
    return engine.prepare(**arguments)


def cpu_per_render(settings, renders):
    steps = settings.steps()
    pairs = numpy.empty((steps, 2))
    status = engine.ctypes.c_int32(0)
    start = time.process_time()
    for _ in range(renders):
        settings.lib.tract_render(settings.handle, engine._pointer(pairs), None, None, engine.ctypes.byref(status))
    return (time.process_time() - start) / renders, pairs


def main():
    rounds = int(sys.argv[1])
    libraries = [os.path.abspath(path) for path in sys.argv[2:]]
    print("%-28s %-40s %8s %8s %7s" % ("word", "library", "min ms", "med ms", "ratio"))
    for name, arguments in bench.words():
        prepared = [settings_for(library, arguments) for library in libraries]
        times = [[] for _ in libraries]
        first = None
        for _ in range(rounds):
            for index, settings in enumerate(prepared):
                seconds, pairs = cpu_per_render(settings, 2)
                times[index].append(seconds)
                if first is None:
                    first = pairs.copy()
                assert pairs.tobytes() == first.tobytes(), "builds differ: %s" % libraries[index]
        base = numpy.min(times[0])
        for index, library in enumerate(libraries):
            least = numpy.min(times[index])
            median = numpy.median(times[index])
            print("%-28s %-40s %8.2f %8.2f %7.3f" % (name, library[-40:], 1e3 * least, 1e3 * median, least / base))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
