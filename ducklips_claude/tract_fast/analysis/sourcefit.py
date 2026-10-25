"""Fit the voice source to human vowels: blockwise Levenberg-Marquardt
(the owner's blockwise ALS) on fold, posture, glottal-noise and variation
settings, so the sustained /aa iy ow/ (1 s each, skeleton tongue 3) match
bdl's stressed vowels in the source measures (source_measures.py):
H1-H2, harmonics-to-noise, jitter, spectral tilt and formant valleys.
Each step is a damped, radius-capped Gauss-Newton step on one block,
kept only if the whole cost falls.
usage: sourcefit.py WORKDIR BASE.params SWEEPS"""
import json, math, os, subprocess, sys
from multiprocessing import Pool
import numpy, soundfile
sys.path.insert(0, "/tmp/es")
import source_measures

SPEAK = os.path.expanduser("~/ducklips/tract/build_x/speak")
VOWELS = ("ah", "ee", "oh")
SEED = "73021"
# bdl medians (source_measures.py arctic, 400 files): aa, iy, ow.
TARGET = {"ah": dict(h1h2=-2.69, hnr=19.13, jitter=0.42, tilt=-8.91, valley=20.83),
          "ee": dict(h1h2=-6.93, hnr=16.09, jitter=0.65, tilt=9.46, valley=21.54),
          "oh": dict(h1h2=-3.26, hnr=19.63, jitter=0.37, tilt=-8.82, valley=24.74)}
SCALE = dict(h1h2=2.0, hnr=3.0, jitter=0.15, tilt=3.0, valley=4.0)
# block.key: (kind, low, high); log kinds move multiplicatively, linear ones additively (step units).
LEAVES = {
    "posture.voicing_lower": ("linear", -0.02, 0.03, 0.004),
    "posture.voicing_upper": ("linear", -0.02, 0.03, 0.004),
    "posture.voicing_chink": ("log", 1e-5, 0.02, None),
    "fold.cover_mass": ("log", 0.02, 0.4, None),
    "fold.cover_stiffness": ("log", 5000, 200000, None),
    "fold.cover_loss": ("log", 0.01, 1.0, None),
    "fold.contact_softness": ("log", 5e-5, 0.005, None),
    "glottal_noise.aspiration": ("log", 0.5, 80, None),
    "glottal_noise.fold_turbulence": ("log", 0.05, 20, None),
    "glottal_noise.closing_turbulence": ("linear", 0.0, 20, 1.0),
    "variation.tension": ("linear", 0.0, 0.2, 0.01),
    "variation.posture": ("linear", 0.0, 0.2, 0.01),
}
BLOCKS = [["posture.voicing_lower", "posture.voicing_upper", "posture.voicing_chink"],
          ["fold.cover_mass", "fold.cover_stiffness", "fold.cover_loss", "fold.contact_softness"],
          ["glottal_noise.aspiration", "glottal_noise.fold_turbulence", "glottal_noise.closing_turbulence"],
          ["variation.tension", "variation.posture"]]

def read_tree(path):
    return open(path).read()

def value_of(text, leaf):
    block, key = leaf.split(".")
    inside = False
    for line in text.split("\n"):
        if line.startswith(block + " {"):
            inside = True
        elif inside and line.startswith("}"):
            inside = False
        elif inside and line.strip().split()[:1] == [key]:
            return float(line.split()[1])
    return None

def with_values(text, values):
    lines = text.split("\n")
    for leaf, value in values.items():
        block, key = leaf.split(".")
        inside = False
        done = False
        for i, line in enumerate(lines):
            if line.startswith(block + " {"):
                inside = True
            elif inside and line.startswith("}"):
                if not done:
                    lines.insert(i, "  %s %.10g" % (key, value))
                    done = True
                inside = False
                break
            elif inside and line.strip().split()[:1] == [key]:
                lines[i] = "  %s %.10g" % (key, value)
                done = True
        if not done:
            lines.append("%s {\n  %s %.10g\n}" % (block, key, value))
    return "\n".join(lines)

def to_unit(leaf, value):
    kind, low, high, step = LEAVES[leaf]
    return math.log(value) if kind == "log" else value / step

def from_unit(leaf, unit):
    kind, low, high, step = LEAVES[leaf]
    value = math.exp(unit) if kind == "log" else unit * step
    return min(max(value, low), high)

def render(task):
    workdir, tag, vowel, tree = task
    params = os.path.join(workdir, "p_%s_%s.params" % (tag, vowel))
    with open(params, "w") as handle:
        handle.write(tree)
    wav = os.path.join(workdir, "w_%s_%s.wav" % (tag, vowel))
    env = dict(os.environ); env["FOLDTUBE"] = "1"; env["PARAMS"] = params
    done = subprocess.run([SPEAK, os.path.join(workdir, vowel + ".score"), wav, SEED], env=env,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = None
    if done.returncode == 0:
        x, r = soundfile.read(wav)
        seg = x[int(0.85 * r):int(1.55 * r)]
        if numpy.all(numpy.isfinite(seg)) and numpy.abs(seg).max() > 1e-4:
            rows = [source_measures.measure(seg[i:i + int(0.23 * r)], r) for i in range(0, len(seg) - int(0.23 * r), int(0.23 * r))]
            rows = [row for row in rows if row]
            if rows:
                result = {k: float(numpy.median([row[k] for row in rows])) for k in rows[0]}
    for path in (params, wav):
        if os.path.exists(path):
            os.remove(path)
    return tag, vowel, result

def residuals(measures):
    out = []
    for vowel in VOWELS:
        m = measures.get(vowel)
        if m is None:
            return None
        for key, scale in SCALE.items():
            out.append((m[key] - TARGET[vowel][key]) / scale)
    return numpy.array(out)

def evaluate(pool, workdir, trees):
    tasks = [(workdir, tag, vowel, tree) for tag, tree in trees.items() for vowel in VOWELS]
    found = {tag: {} for tag in trees}
    for tag, vowel, result in pool.imap_unordered(render, tasks):
        if result is not None:
            found[tag][vowel] = result
    return {tag: (residuals(found[tag]), found[tag]) for tag in trees}

def main():
    workdir, base_path, sweeps = sys.argv[1], sys.argv[2], int(sys.argv[3])
    base = read_tree(base_path)
    units = {leaf: to_unit(leaf, value_of(base, leaf) if value_of(base, leaf) is not None else LEAVES[leaf][1] if LEAVES[leaf][0] == "log" else 0.0) for leaf in LEAVES}
    for leaf in LEAVES:
        if LEAVES[leaf][0] == "log" and value_of(base, leaf) in (None, 0.0):
            units[leaf] = math.log(max(LEAVES[leaf][1] * 2, 1e-9))
    def tree_of(u):
        return with_values(base, {leaf: from_unit(leaf, u[leaf]) for leaf in LEAVES})
    log = open(os.path.join(workdir, "sourcefit.log"), "a")
    with Pool(8) as pool:
        start = evaluate(pool, workdir, {"start": tree_of(units)})["start"]
        best_r, best_m = start
        best = float(numpy.mean(best_r ** 2))
        print("start cost %.3f %s" % (best, json.dumps({v: {k: round(x, 2) for k, x in m.items()} for v, m in best_m.items()})), flush=True)
        damping = {i: 1.0 for i in range(len(BLOCKS))}
        radius = {i: 0.5 for i in range(len(BLOCKS))}
        for sweep in range(sweeps):
            for index, block in enumerate(BLOCKS):
                delta = 0.15
                trees = {}
                for leaf in block:
                    u = dict(units); u[leaf] += delta
                    trees[leaf] = tree_of(u)
                result = evaluate(pool, workdir, trees)
                jac = []
                for leaf in block:
                    r, _ = result[leaf]
                    jac.append((r - best_r) / delta if r is not None else numpy.zeros_like(best_r))
                J = numpy.array(jac).T
                kept = False
                for attempt in range(3):
                    A = J.T @ J + damping[index] * numpy.eye(len(block))
                    step = -numpy.linalg.solve(A, J.T @ best_r)
                    step = numpy.clip(step, -radius[index], radius[index])
                    u = dict(units)
                    for leaf, s in zip(block, step):
                        u[leaf] += s
                    trial_r, trial_m = evaluate(pool, workdir, {"trial": tree_of(u)})["trial"]
                    trial = float(numpy.mean(trial_r ** 2)) if trial_r is not None else float("inf")
                    if trial < best:
                        units, best, best_r, best_m = u, trial, trial_r, trial_m
                        damping[index] = max(damping[index] / 3, 1e-3)
                        kept = True
                        break
                    damping[index] *= 4
                    radius[index] /= 2
                record = dict(sweep=sweep, block=index, kept=kept, cost=best, values={leaf: from_unit(leaf, units[leaf]) for leaf in LEAVES}, measures=best_m)
                log.write(json.dumps(record) + "\n"); log.flush()
                print("sweep %d block %d %s cost %.3f | %s" % (sweep, index, "kept" if kept else "----", best,
                      " ".join("%s:h1h2 %.1f hnr %.1f jit %.2f tilt %.1f" % (v, m["h1h2"], m["hnr"], m["jitter"], m["tilt"]) for v, m in best_m.items())), flush=True)
                with open(os.path.join(workdir, "best.params"), "w") as handle:
                    handle.write(tree_of(units))

if __name__ == "__main__":
    main()
