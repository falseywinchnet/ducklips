"""Coordinate search of vowel postures (vplace, varea, vlips) for the sheet,
all tuned vowels moved together each trial and each kept or reverted on its
own distance to bdl (vtune.py, 60 sentences). Writes the best table.
usage: vopt.py BASE_TABLE OUT_TABLE ROUNDS"""
import re, subprocess, sys

START = {  # vowel: [vplace, varea, vlips]; START env "tract" picks the segment engine's table values
    "AA": [7.0, 0.724, 0.988], "AE": [9.0, 2.1, 9.0], "AH": [8.75, 0.722, 0.599], "AO": [5.75, 0.137, 0.228],
    "OW": [5.75, 0.526, 0.356], "UW": [13.13, 3.5, 0.267], "IY": [12.25, 0.321, 4.9], "EH": [9.5, 1.4, 6.0],
    "IH": [10.75, 0.697, 4.5], "EY": [9.56, 0.325, 5.0], "UH": [8.62, 0.479, 0.769], "ER": [6.5, 1.6, 0.6]}
import os
if os.environ.get("START") == "tract":
    START = {"AA": [7.0, 0.724, 5.0], "AE": [9.0, 0.75, 9.0], "AH": [8.5, 0.722, 5.13], "AO": [5.0, 0.3, 1.0],
             "OW": [5.0, 0.526, 0.8], "UW": [13.88, 2.5, 0.6], "IY": [13.75, 0.179, 4.9], "EH": [9.5, 1.4, 6.0],
             "IH": [11.25, 0.697, 4.5], "EY": [9.56, 0.26, 5.0], "UH": [9.75, 0.5, 0.4], "ER": [6.5, 2.0, 0.6]}
STEPS = [("vplace", 0, "add", 0.6), ("varea", 1, "mul", 1.3), ("vlips", 2, "mul", 1.4)]
base_lines = [l.rstrip("\n") for l in open(sys.argv[1])]


def write(table, values):
    lines = []
    for l in base_lines:
        f = l.split()
        if f and f[0] in values:
            keep, k = [], 1
            while k + 1 < len(f):
                if f[k] not in ("vplace", "varea", "vlips"):
                    keep += [f[k], f[k + 1]]
                k += 2
            v = values[f[0]]
            l = " ".join([f[0]] + keep + ["vplace", "%.3f" % v[0], "varea", "%.3f" % v[1], "vlips", "%.3f" % v[2]])
        lines.append(l)
    names = {l.split()[0] for l in lines if l.split()}
    for n, v in values.items():
        if n not in names:
            lines.append("%s vplace %.3f varea %.3f vlips %.3f" % (n, v[0], v[1], v[2]))
    open(table, "w").write("\n".join(lines) + "\n")


def score(values, name):
    write("/tmp/es/vopt_try.txt", values)
    out = subprocess.run(["/Users/joshuahkuttenkuler/miniforge3/bin/python3", "/tmp/vtune.py", "/tmp/es/vopt_try.txt", name] +
                         [v.lower() for v in values], capture_output=True, text=True).stdout
    d = {}
    for l in out.splitlines():
        m = re.match(r"(\w+)\s+n\s+\d+.*dist ([\d.]+) st", l)
        if m:
            d[m.group(1).upper()] = float(m.group(2))
    return d


best = {k: list(v) for k, v in START.items()}
current = score(best, "opt")
print("start", {k: round(v, 1) for k, v in current.items()}, flush=True)
for r in range(int(sys.argv[3])):
    for key, i, op, step in STEPS:
        for sign in (+1, -1):
            trial = {k: list(v) for k, v in best.items()}
            for k in trial:
                trial[k][i] = trial[k][i] + sign * step if op == "add" else trial[k][i] * (step if sign > 0 else 1 / step)
            d = score(trial, "opt")
            for k in best:
                if k in d and k in current and d[k] < current[k] - 0.3:
                    best[k] = trial[k]
                    current[k] = d[k]
            print("round %d %s %+d: %s" % (r, key, sign, " ".join("%s %.1f" % (k, current[k]) for k in sorted(current))), flush=True)
    write(sys.argv[2], best)
print("final", best, "mean %.2f" % (sum(current.values()) / len(current)))
