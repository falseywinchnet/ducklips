"""Tune the phone inventory against the four-view baseline: each candidate
override is rendered over ARCTIC a0001-a0020, analysed in all four views
(enhanced STFT, constant-Q, time-causal Gabor, Wigner-Ville) and scored by
the four-view total (mean shape error over the views plus fricative and
stop-release level errors, tools/compare_arctic.py). The best value of each
step is kept before the next step; a0021-a0040 is scored at the end as a
hold-out check.

usage: four_view_tune.py PLAN.json OVERRIDES.txt
  PLAN.json: [{"phone": "SH", "key": "area", "values": [0.1, 0.15, 0.2]}, ...]
  Global settings use the phone "*" with keys aspiration, frication, pressure.
"""
import contextlib
import io
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
import tune_phones as tp

WORK = os.path.join(ca.ROOT, "out", "four_view_tune")


def say_arguments(overrides_path, globals_):
    extra = ["--phones", overrides_path]
    for key in sorted(globals_):
        extra += ["--" + key, "%g" % globals_[key]]
    return extra


def score(settings, globals_, first, last, label):
    folder = os.path.join(WORK, label)
    os.makedirs(folder, exist_ok=True)
    overrides_path = os.path.join(folder, "phones.txt")
    tp.write_overrides(overrides_path, settings)
    records = ca.compare(first, last, folder, say_arguments(overrides_path, globals_))
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink):
        total = ca.four_view_summary(records)
    shutil.rmtree(folder, ignore_errors=True)
    return total


def main():
    plan = json.load(open(sys.argv[1]))
    output = sys.argv[2]
    settings = tp.read_overrides(output)
    globals_ = settings.pop("*", {})
    best = score(settings, globals_, 1, 20, "baseline")
    print("baseline four-view total %.3f" % best, flush=True)
    for step in plan:
        phone, key = step["phone"], step["key"]
        results = []
        for value in step["values"]:
            trial_settings = {p: dict(v) for p, v in settings.items()}
            trial_globals = dict(globals_)
            if phone == "*":
                trial_globals[key] = value
            else:
                trial_settings.setdefault(phone, {})[key] = value
            total = score(trial_settings, trial_globals, 1, 20, "%s_%s_%g" % (phone.replace("*", "all"), key, value))
            results.append((total, value))
            print("  %-3s %s=%-7g total %.3f" % (phone, key, value, total), flush=True)
        total, value = min(results)
        if total < best:
            best = total
            if phone == "*":
                globals_[key] = value
            else:
                settings.setdefault(phone, {})[key] = value
            print("%s %s -> %g (total %.3f)" % (phone, key, value, best), flush=True)
        else:
            print("%s %s unchanged (best %.3f)" % (phone, key, best), flush=True)
        saved = dict(settings)
        if globals_:
            saved["*"] = globals_
        tp.write_overrides(output, saved)
    print("hold-out a0021-a0040 four-view total %.3f" % score(settings, globals_, 21, 40, "holdout"), flush=True)


if __name__ == "__main__":
    main()
