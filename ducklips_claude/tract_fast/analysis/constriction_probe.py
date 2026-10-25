"""The articulated geometry's narrowest area per phone: render sentences with
AREALOG and take, over the middle 40% of each phone, the smallest rest area
from node 7 up to the lip opening (the sliver past the lips, 0.02 cm^2, is
not in the log; speak's TRACE line scans two nodes past the lips and so
reported 0.020 for every fricative, 2026-10-18). Medians per phone.
usage: constriction_probe.py OUTDIR IDS PARAMS ORGANS PLAN"""
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fidelity
import render_set

ROOT = os.path.expanduser("~/ducklips")
FIRST_NODE = 7
SHARE = 0.3


def render_one(job):
    name, text, outdir, params, organs, plan = job
    score = os.path.join(outdir, name + ".score")
    areas = os.path.join(outdir, name + ".areas")
    say = os.environ.get("SAY_BIN", os.path.join(ROOT, "voice", "build_dev", "say"))
    speak = os.environ.get("SPEAK_BIN", os.path.join(ROOT, "tract", "build_dev", "speak"))
    subprocess.run([say, text, "/dev/null", "--score", score, "--plan-only", "--rate", "0.75", "--organs", organs,
                    "--plan-params", plan], cwd=os.path.join(ROOT, "voice"), capture_output=True)
    environment = dict(os.environ, PARAMS=params, KINEMATIC="1", AREALOG=areas)
    subprocess.run([speak, score, "/dev/null", render_set.SEED], env=environment, capture_output=True)
    return name


def narrowest(outdir, names):
    """{phone: [(area, node), ...]} over the named renders."""
    found = {}
    for name in names:
        rows = [list(map(float, line.split())) for line in open(os.path.join(outdir, name + ".areas"))]
        times = numpy.array([r[0] for r in rows])
        for phone, start, end in fidelity.score_segments(os.path.join(outdir, name + ".score")):
            middle_start = start + SHARE * (end - start)
            middle_end = end - SHARE * (end - start)
            for index in numpy.nonzero((times >= middle_start) & (times <= middle_end))[0]:
                nodes = rows[index][1 + FIRST_NODE:]
                if not nodes:
                    continue
                place = int(numpy.argmin(nodes))
                found.setdefault(phone.rstrip("012-").upper(), []).append((nodes[place], place + FIRST_NODE))
    return found


def render(outdir, ids, params, organs, plan):
    os.makedirs(outdir, exist_ok=True)
    table = render_set.texts()
    names = [name for name in fidelity.id_range(ids) if name in table]
    jobs = [(name, table[name], outdir, os.path.abspath(params), os.path.abspath(organs), os.path.abspath(plan))
            for name in names]
    with ThreadPoolExecutor(int(os.environ.get("WORKERS", "9"))) as pool:
        list(pool.map(render_one, jobs))
    return names


def main():
    outdir, ids, params, organs, plan = sys.argv[1:6]
    found = narrowest(outdir, render(outdir, ids, params, organs, plan))
    for phone in sorted(found):
        areas = numpy.array([a for a, _ in found[phone]])
        nodes = numpy.array([n for _, n in found[phone]])
        print("%-3s frames %5d  area median %.3f  quartiles %.3f %.3f  node %d" % (
            phone, len(areas), numpy.median(areas), *numpy.percentile(areas, [25, 75]), numpy.median(nodes)))


if __name__ == "__main__":
    main()
