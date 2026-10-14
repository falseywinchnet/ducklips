"""Merge tuned phone overrides into data/phones_en.txt, keeping its comments
and line order; keys named with --drop PHONE.KEY are left as they were.
usage: merge_overrides.py TUNED.txt [--drop PHONE.KEY ...]"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tune_phones as tp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TARGET = os.path.join(ROOT, "data", "phones_en.txt")


def main():
    tuned = tp.read_overrides(sys.argv[1])
    tuned.pop("*", None)
    drops = [sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == "--drop"]
    current = tp.read_overrides(TARGET)
    for item in drops:
        phone, key = item.split(".")
        if key in current.get(phone, {}):
            tuned[phone][key] = current[phone][key]
        else:
            tuned.get(phone, {}).pop(key, None)
    lines = open(TARGET).read().splitlines()
    seen = set()
    out = []
    for line in lines:
        parts = line.split()
        if parts and not parts[0].startswith("#") and parts[0] in tuned:
            phone = parts[0]
            seen.add(phone)
            values = tuned[phone]
            line = phone + "".join(" %s %g" % (k, values[k]) for k in sorted(values))
        out.append(line)
    for phone in sorted(set(tuned) - seen):
        if tuned[phone]:
            out.append(phone + "".join(" %s %g" % (k, tuned[phone][k]) for k in sorted(tuned[phone])))
    open(TARGET, "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
