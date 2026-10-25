"""Closed-form pitch for a profile: render the reference's sentences, take
our F0 percentiles (fidelity's: pyin frames in loud voiced speech,
semitones re 100 Hz) against the speaker's, and move the plan's
pitch.pitch_median by the median's gap, pitch.pitch_spread by the
interquartile ratio, and pitch.lowest_pitch by what is left of the 5th
percentile's gap; a few passes (renders follow the plan only roughly: the
folds' pitch, jitter and glides between vowel targets).
(bdl4 against bdl a0101-a0120: ours -1.44 / 1.66 / 3.26 / 4.26 / 5.86 st at
p5/25/50/75/95, his -0.14 / 2.26 / 3.56 / 4.86 / 6.26.)
usage: pitch_calibrate.py REFDIR IDS PARAMS ORGANS PLAN_IN PLAN_OUT [PASSES]"""
import os
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fidelity
import render_set
import slider_server


def percentiles(described):
    return numpy.array([described["pitch"]["f0_p%02d" % q] for q in (5, 25, 50, 75, 95)])


def main():
    ref_dir, ids_text, params, organs, plan_in, plan_out = sys.argv[1:7]
    passes = int(sys.argv[7]) if len(sys.argv) > 7 else 3
    ids = fidelity.id_range(ids_text)
    reference = percentiles(fidelity.describe(fidelity.arctic_set(ref_dir, ids)))
    text = open(plan_in).read()
    folder = tempfile.mkdtemp(dir="/tmp/es")
    table = render_set.texts()
    for step in range(passes + 1):
        current = os.path.join(folder, "plan_%d.params" % step)
        open(current, "w").write(text)
        jobs = [(name, table[name], folder, os.path.abspath(params), os.path.abspath(organs), current, {"KINEMATIC": "1"})
                for name in ids if name in table]
        with ThreadPoolExecutor(int(os.environ.get("WORKERS", "9"))) as pool:
            list(pool.map(render_set.render, jobs))
        ours = percentiles(fidelity.describe(fidelity.ours_set(folder, ids)))
        error = ours - reference
        print("pass %d: ours %s | gap %s st (rms %.2f)" % (step, " ".join("%.2f" % v for v in ours),
                                                           " ".join("%+.2f" % v for v in error),
                                                           numpy.sqrt(numpy.mean(error ** 2))), flush=True)
        if step == passes:
            break
        median = slider_server.read_value(text, "pitch.pitch_median") or 122.0
        spread = slider_server.read_value(text, "pitch.pitch_spread") or 1.5
        lowest = slider_server.read_value(text, "pitch.lowest_pitch") or 80.0
        ratio = (reference[3] - reference[1]) / max(ours[3] - ours[1], 0.1)
        changes = {"pitch.pitch_median": float(median * 2.0 ** (-error[2] / 12.0)),
                   "pitch.pitch_spread": float(numpy.clip(spread * ratio, 0.5, 3.0)),
                   # (The tail after the median's shift.)
                   "pitch.lowest_pitch": float(numpy.clip(lowest * 2.0 ** (-(error[0] - error[2]) / 12.0), 40.0, 200.0))}
        text = slider_server.with_values(text, changes)
        print("  " + " ".join("%s %.3f" % (k.split(".")[1], v) for k, v in sorted(changes.items())), flush=True)
    open(plan_out, "w").write("# pitch_calibrate.py against %s %s\n" % (os.path.basename(ref_dir), ids_text) + text)


if __name__ == "__main__":
    main()
