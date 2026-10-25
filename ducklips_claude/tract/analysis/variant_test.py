"""Held-out fidelity of a few variants of one voice, with each variant's
largest per-class spectrum errors (the classes whose bands moved most).
A variant is NAME:KEY=VALUE,KEY=VALUE; a key `plan.BLOCK.X` sets plan parameter X,
`organs.PHONE.FIELD` one organ posture field, any other a tree slot.
usage: variant_test.py REFDIR IDS PARAMS ORGANS PLAN VARIANT... (the first
variant is usually `base:`)"""
import os
import pickle
import re
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fidelity
import render_set
import slider_server

WORK = "/tmp/es/variants"
# Each variant rendered with these seeds and described as one set: with one
# seed the phase group of a single voice ranged 11.3-16.0 dB (fidelity v4).
SEEDS = os.environ.get("SEEDS", "73021,11").split(",")


def organs_with(text, changes):
    lines = text.split("\n")
    for key, value in changes.items():
        phone, field = key.split(".", 1)
        for index, line in enumerate(lines):
            if line.split(" ", 1)[0] == phone:
                lines[index] = re.sub(r"\b%s=[-0-9.e]+" % re.escape(field), "%s=%s" % (field, value), line)
    return "\n".join(lines)


def build(name, spec, params, organs, plan):
    folder = os.path.join(WORK, name)
    os.makedirs(folder, exist_ok=True)
    tree, plan_changes, organ_changes = {}, {}, {}
    for item in filter(None, spec.split(",")):
        key, value = item.split("=", 1)
        if key.startswith("plan."):
            plan_changes[key[5:]] = float(value)
        elif key.startswith("organs."):
            organ_changes[key[7:]] = value
        else:
            tree[key] = float(value)
    paths = []
    for text, changes, writer, file_name in ((open(params).read(), tree, slider_server.with_values, "voice.params"),
                                             (open(plan).read(), plan_changes, slider_server.with_values, "plan.params"),
                                             (open(organs).read(), organ_changes, organs_with, "organs.txt")):
        path = os.path.join(folder, file_name)
        open(path, "w").write(writer(text, changes) if changes else text)
        paths.append(path)
    return folder, paths


def reference_of(ref_dir, ids):
    cache = os.path.join(WORK, "ref_%s_%s.pkl" % (os.path.basename(ref_dir), ids.replace("/", "s")))
    if os.path.exists(cache):
        return pickle.load(open(cache, "rb"))
    described = fidelity.describe(fidelity.arctic_set(ref_dir, fidelity.id_range(ids)))
    pickle.dump(described, open(cache, "wb"))
    return described


def class_errors(reference, ours):
    out = {}
    for key, r in reference["spectrum"].items():
        if key in ours["spectrum"]:
            c = key.rsplit("_", 1)[0]
            out.setdefault(c, []).append(ours["spectrum"][key] - r)
    return {c: numpy.array(v) for c, v in out.items()}


def main():
    ref_dir, ids, params, organs, plan = sys.argv[1:6]
    os.makedirs(WORK, exist_ok=True)
    reference = reference_of(ref_dir, ids)
    base_errors = None
    for variant in sys.argv[6:]:
        name, spec = variant.split(":", 1)
        folder, (p, plan_path, o) = build(name, spec, params, organs, plan)
        os.environ.setdefault("SPEAK_BIN", os.path.expanduser("~/ducklips/tract/build_dev/speak"))
        os.environ.setdefault("SAY_BIN", os.path.expanduser("~/ducklips/voice/build_dev/say"))
        table = render_set.texts()
        utterances = []
        for seed in SEEDS:
            render_set.SEED = seed
            seed_folder = os.path.join(folder, seed)
            os.makedirs(seed_folder, exist_ok=True)
            jobs = [(name_id, table[name_id], seed_folder, p, o, plan_path, {"KINEMATIC": "1"})
                    for name_id in fidelity.id_range(ids) if name_id in table]
            with ThreadPoolExecutor(int(os.environ.get("WORKERS", "3"))) as pool:
                list(pool.map(render_set.render, jobs))
            utterances += fidelity.ours_set(seed_folder, fidelity.id_range(ids), ref_dir)
        ours = fidelity.describe(utterances)
        errors = class_errors(reference, ours)
        line = "%-14s %s" % (name, fidelity.report(fidelity.score(reference, ours)))
        if base_errors is None:
            base_errors = errors
            worst = sorted(errors, key=lambda c: -numpy.sqrt(numpy.mean(errors[c] ** 2)))[:6]
            line += "\n    worst classes (rms dB): " + " ".join("%s %.1f" % (c, numpy.sqrt(numpy.mean(errors[c] ** 2))) for c in worst)
        else:
            moved = sorted((c for c in errors if c in base_errors),
                           key=lambda c: numpy.sqrt(numpy.mean(errors[c] ** 2)) - numpy.sqrt(numpy.mean(base_errors[c] ** 2)))
            def show(c):
                return "%s %+.1f" % (c, numpy.sqrt(numpy.mean(errors[c] ** 2)) - numpy.sqrt(numpy.mean(base_errors[c] ** 2)))
            line += "\n    better: " + " ".join(show(c) for c in moved[:5]) + " | worse: " + " ".join(show(c) for c in moved[-3:])
        print(line, flush=True)


if __name__ == "__main__":
    main()
