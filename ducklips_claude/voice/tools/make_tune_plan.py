"""Builds a full tuning pass around the current phone inventory: the
built-in table in src/planner.cpp overlaid with data/phones_en.txt.
usage: make_tune_plan.py OUT.json"""
import json
import re
import sys

source = open("src/planner.cpp").read()
postures = {m.group(1): [float(m.group(i)) for i in range(2, 6)]
            for m in re.finditer(r"constexpr Posture k(\w+) = \{([\d.]+), ([\d.]+), ([\d.]+), ([\d.]+)\}", source)}
phones = {}
for m in re.finditer(r'\{"(\w+)", Kind::(\w+), (true|false), ([\d.]+), ([\d.]+), k(\w+), k(\w+), Organ::(\w+), ([\d.]+), ([\d.]+), ([\d.]+)\}', source):
    name, kind = m.group(1), m.group(2)
    vowel = postures.get(m.group(6), [0, 0, 0, 0])
    phones[name] = dict(kind=kind, organ=m.group(8), inherent=float(m.group(4)), minimum=float(m.group(5)),
                        place=float(m.group(9)), area=float(m.group(10)),
                        vplace=vowel[0], varea=vowel[1], vlips=vowel[2])
for line in open("data/phones_en.txt"):
    parts = line.split()
    if parts and not parts[0].startswith("#"):
        for i in range(1, len(parts) - 1, 2):
            phones[parts[0]][parts[i]] = float(parts[i + 1])

plan = []
for name, p in phones.items():
    if name == "DX":
        continue
    steps = [("inherent", [round(p["inherent"] * f) for f in (0.7, 0.85, 1.0, 1.2)])]
    if p["kind"] in ("vowel", "diphthong") or name in ("R", "W", "Y"):
        steps += [("varea", [round(p["varea"] * f, 3) for f in (0.75, 1.0, 1.33)]),
                  ("vplace", [round(p["vplace"] + d, 2) for d in (-0.6, 0.0, 0.6)]),
                  ("vlips", [round(p["vlips"] * f, 3) for f in (0.75, 1.0, 1.33)])]
    if p["organ"] in ("tip", "dorsum"):
        steps.append(("place", [round(p["place"] + d, 2) for d in (-0.4, 0.0, 0.4)]))
    if p["organ"] in ("tip", "dorsum", "lips") and p["area"] > 0:
        steps.append(("area", [round(p["area"] * f, 3) for f in (0.7, 1.0, 1.4)]))
    for key, values in steps:
        plan.append({"phone": name, "key": key, "values": values})
json.dump(plan, open(sys.argv[1], "w"))
print(len(plan), "steps,", sum(len(s["values"]) for s in plan), "candidates")
