"""Render ARCTIC sentences (or any id: text list) through the engine, in
parallel: say plans each (organs, plan parameters), speak renders it.
Writes OUTDIR/ID.score and OUTDIR/ID.wav.
usage: render_set.py OUTDIR IDS PARAMS ORGANS PLAN [ENV=VALUE ...]
IDS: first-last (ARCTIC ids, texts from voice/ref/arctic_bdl/arctic.data)."""
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.expanduser("~/ducklips")
SEED = "73021"


def texts():
    out = {}
    for line in open(os.path.join(ROOT, "voice", "ref", "arctic_bdl", "arctic.data")):
        line = line.strip()
        if line.startswith("( arctic_"):
            name = line.split()[1].replace("arctic_", "")
            text = line[line.index('"') + 1:line.rindex('"')]
            out[name] = text
    return out


def id_range(text):
    """"a0001-a0030", or with a step "a0001-a0593/20" (every 20th: a spread
    across the recording sessions; bdl's F0 on a0201-a0300 sits 0.7-1.6 st
    above a0101-a0120)."""
    span, _, step = text.partition("/")
    first, last = span.split("-")
    return ["%s%04d" % (first[0], n) for n in range(int(first[1:]), int(last[1:]) + 1, int(step or 1))]


def render(job):
    name, text, outdir, params, organs, plan, extra = job
    score = os.path.join(outdir, name + ".score")
    wave = os.path.join(outdir, name + ".wav")
    planned = subprocess.run([os.environ.get("SAY_BIN", os.path.join(ROOT, "voice", "build", "say")), text, "/dev/null", "--score", score,
                              "--plan-only", "--rate", "0.75", "--organs", organs, "--plan-params", plan],
                             cwd=os.path.join(ROOT, "voice"), capture_output=True, text=True)
    if planned.returncode != 0:
        return name, "plan failed"
    environment = dict(os.environ, PARAMS=params)
    environment.update(extra)
    speak = os.environ.get("SPEAK_BIN", os.path.join(ROOT, "tract", "build", "speak"))
    done = subprocess.run([speak, score, wave, SEED], env=environment,
                          capture_output=True, text=True)
    return name, "ok" if done.returncode == 0 else "speak failed"


def render_many(jobs, threads):
    """render() over the jobs, `threads` at a time (each render is its own
    process). (Called from one pool worker while the others waited, a fit's
    held-out check rendered its 40 sentences one by one: about 100 s with
    nine cores idle, several times a stage.)"""
    if threads <= 1:
        return [render(job) for job in jobs]
    with ThreadPoolExecutor(threads) as pool:
        return list(pool.map(render, jobs))


def main():
    outdir, ids, params, organs, plan = sys.argv[1:6]
    extra = dict(item.split("=", 1) for item in sys.argv[6:])
    os.makedirs(outdir, exist_ok=True)
    table = texts()
    jobs = [(name, table[name], outdir, os.path.abspath(params), os.path.abspath(organs), os.path.abspath(plan), extra)
            for name in id_range(ids) if name in table]
    with ThreadPoolExecutor(9) as pool:
        results = list(pool.map(render, jobs))
    failed = [r for r in results if r[1] != "ok"]
    print("rendered %d, failed %d %s" % (len(results) - len(failed), len(failed), failed[:5]))


if __name__ == "__main__":
    main()
