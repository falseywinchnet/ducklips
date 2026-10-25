"""Instructions retired and cycles per render, by /usr/bin/time -l, for
sustain binaries driven exactly as the reference posture_recover.render
drives its tool (its own code: its subprocess call is wrapped). The Mini is
shared and loaded (load 15-20 on 4 performance and 6 efficiency cores), so
wall and CPU times swing by 2x with the core a process lands on; the count
of instructions retired does not, and is the measure of work here.

usage (on the Mini): python3 bindings/bench_instructions.py REPEATS BINARY...
The first binary is the baseline for the ratios."""
import os
import subprocess as _subprocess
import sys

sys.dont_write_bytecode = True

import numpy

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bench  # noqa: E402

reference = bench.reference
MEASURED = []


class _TimedSubprocess:
    """The reference's subprocess, each run under /usr/bin/time -l."""

    @staticmethod
    def run(command, **keywords):
        result = _subprocess.run(["/usr/bin/time", "-l"] + list(command), **keywords)
        found = {}
        for line in result.stderr.decode(errors="replace").splitlines():
            parts = line.split()
            if len(parts) >= 3 and parts[1] == "instructions" and parts[2] == "retired":
                found["instructions"] = float(parts[0])
            elif len(parts) >= 3 and parts[1] == "cycles" and parts[2] == "elapsed":
                found["cycles"] = float(parts[0])
            elif len(parts) >= 2 and parts[1] == "real":
                found["real"] = float(parts[0])
                found["user"] = float(parts[2])
                found["sys"] = float(parts[4])
        MEASURED.append(found)
        return result


reference.subprocess = _TimedSubprocess


def main():
    repeats = int(sys.argv[1])
    binaries = [os.path.abspath(path) for path in sys.argv[2:]]
    print("%-28s %-44s %12s %8s %12s %8s" % ("word", "binary", "instructions", "ratio", "cycles", "ratio"))
    for name, arguments in bench.words():
        rows = []
        for binary in binaries:
            reference.SUSTAIN = binary
            del MEASURED[:]
            outputs = None
            for _ in range(repeats):
                outputs = reference.render(**arguments)
            rows.append((binary, numpy.median([m["instructions"] for m in MEASURED]),
                         numpy.median([m["cycles"] for m in MEASURED]), outputs))
        for binary, instructions, cycles, outputs in rows:
            assert bench.equivalence.same_bytes(outputs[0], rows[0][3][0]), binary
            assert bench.equivalence.same_bytes(outputs[1], rows[0][3][1]), binary
            print("%-28s %-44s %12.4g %8.3f %12.4g %8.3f" % (name, binary[-44:], instructions, instructions / rows[0][1],
                                                            cycles, cycles / rows[0][2]))
        sys.stdout.flush()
    reference.SUSTAIN = os.path.join(bench.equivalence.REFERENCE, "build_dev", "sustain")


if __name__ == "__main__":
    main()
