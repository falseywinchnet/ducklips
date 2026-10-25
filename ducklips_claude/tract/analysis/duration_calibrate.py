"""Closed-form timing for a profile: plan the reference's sentences
(say --plan-only, no rendering), pair each planned phone with the speaker's
labelled phone it aligns with (the same word and place; fidelity's
aligned_names alignment), and scale the plan's duration factor of each of
our phones by the geometric mean ratio (reference / ours) over its pairs, a
few passes (the planner's other rules interact). Only phones with at least
MIN pairs. (Pooled by name instead, our weak-form AX and bdl's citation
/ae/ in "that", "and" made the two readings' /ae/ different words: held-out
/ae/ came out 0.63 of bdl's after calibration.)
usage: duration_calibrate.py REFDIR IDS ORGANS PLAN_IN PLAN_OUT [PASSES]"""
import difflib
import os
import subprocess
import sys
import tempfile

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fidelity
import render_set
import slider_server

ROOT = os.path.expanduser("~/ducklips")
MIN = 5
WEAK = "weak forms"
PHONES = ["IY", "IH", "EY", "EH", "AE", "AA", "AO", "OW", "UH", "UW", "AH", "AX", "ER", "AY", "AW", "OY", "P", "B",
          "T", "D", "DX", "K", "G", "M", "N", "NG", "F", "V", "TH", "DH", "S", "Z", "SH", "ZH", "CH", "JH", "L", "R",
          "W", "Y", "HH"]


def pairs(ids, organs, plan, folder, ref_dir):
    """{our phone: [(ours seconds, reference seconds), ...]} over aligned pairs."""
    table = render_set.texts()
    out = {}
    for name in ids:
        lab = os.path.join(ref_dir, "lab", "arctic_%s.lab" % name)
        if not os.path.exists(lab) or name not in table:
            continue
        score = os.path.join(folder, name + ".score")
        subprocess.run([os.environ.get("SAY_BIN", os.path.join(ROOT, "voice", "build", "say")), table[name], "/dev/null", "--score", score,
                        "--plan-only", "--rate", "0.75", "--organs", organs, "--plan-params", plan],
                       cwd=os.path.join(ROOT, "voice"), capture_output=True)
        ours = fidelity.score_segments(score)
        spoken = [(p, a, b) for p, a, b in fidelity.arctic_labels(lab) if p.lower() not in fidelity.SILENT_LABELS]
        matcher = difflib.SequenceMatcher(a=[p.lower() for p, _, _ in spoken], b=[fidelity.base_name(p) for p, _, _ in ours],
                                          autojunk=False)
        for tag, i1, i2, j1, j2 in matcher.get_opcodes():
            if tag == "equal" or (tag == "replace" and i2 - i1 == j2 - j1):
                for i, j in zip(range(i1, i2), range(j1, j2)):
                    phone = fidelity.base_name(ours[j][0]).upper()
                    # Our schwa against a full vowel in bdl's labels is a
                    # weak form: its own factor (timing.function_vowel_scale).
                    if phone == "AX" and spoken[i][0].lower() not in ("ax", "ah", "ih", "er"):
                        phone = WEAK
                    out.setdefault(phone, []).append((ours[j][2] - ours[j][1], spoken[i][2] - spoken[i][1]))
    return out


def main():
    ref_dir, ids_text, organs, plan_in, plan_out = sys.argv[1:6]
    passes = int(sys.argv[6]) if len(sys.argv) > 6 else 3
    ids = fidelity.id_range(ids_text)
    text = open(plan_in).read()
    folder = tempfile.mkdtemp(dir="/tmp/es")
    for step in range(passes):
        current = os.path.join(folder, "plan_%d.params" % step)
        open(current, "w").write(text)
        found = pairs(ids, os.path.abspath(organs), current, folder, ref_dir)
        changes = {}
        errors = []
        for phone in PHONES + [WEAK]:
            matched = found.get(phone, [])
            if len(matched) < MIN:
                continue
            o = numpy.array([m[0] for m in matched])
            r = numpy.array([m[1] for m in matched])
            ratio = numpy.exp(numpy.mean(numpy.log(numpy.maximum(r, 0.005))) - numpy.mean(numpy.log(numpy.maximum(o, 0.005))))
            errors.append(numpy.log(ratio) ** 2)
            key = "timing.function_vowel_scale" if phone == WEAK else "duration." + phone
            factor = slider_server.read_value(text, key) or 1.0
            changes[key] = float(numpy.clip(factor * ratio, 0.3, 3.0))
        print("pass %d: rms log ratio %.3f over %d phones" % (step, numpy.sqrt(numpy.mean(errors)), len(errors)))
        text = slider_server.with_values(text, changes)
    open(plan_out, "w").write("# duration_calibrate.py against %s %s\n" % (os.path.basename(ref_dir), ids_text) + text)
    print(" ".join("%s %.2f" % (k.split(".")[1], v) for k, v in sorted(changes.items())))


if __name__ == "__main__":
    main()
