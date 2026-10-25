"""Closed-form source for a profile: render the reference's sentences and
move, per pass,
  pitch (the plan): pitch.pitch_median by the F0 median's gap,
    pitch.pitch_spread by the interquartile ratio, pitch.lowest_pitch by
    what is left of the 5th percentile's gap (as pitch_calibrate.py);
  phase (the tree): kinematic.cycle_lag (each cycle's lag drawn afresh)
    toward the speaker's phase scatter at 0-1 kHz, and
    vtl.aspiration_gain toward the scatter at 1-4 kHz.
(bdl4: scatter 0.20 / 0.35 / 0.68 rad against bdl's 0.37 / 0.49 / 0.57.
Cycle lag 0.6, shape 0.35 and aspiration gain 3 (from 11.8) brought 0.35 /
0.46 / 0.61: held-out phase 10.4 -> 24.6 dB, combined 7.5 -> 8.0, two
seeds; the low bands' scatter grows with the lag, the top band's with the
aspiration.)
usage: source_calibrate.py REFDIR IDS PARAMS_IN ORGANS PLAN_IN PARAMS_OUT PLAN_OUT [PASSES]"""
import os
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fidelity
import render_set
import slider_server

PERCENTILES = (5, 25, 50, 75, 95)
DAMPING = 0.6


def pitch_of(described):
    return numpy.array([described["pitch"]["f0_p%02d" % q] for q in PERCENTILES])


def phase_of(described):
    return numpy.array([described["phase"].get("vowels_all_%d" % k, numpy.nan) for k in range(3)])


def main():
    ref_dir, ids_text, params_in, organs, plan_in, params_out, plan_out = sys.argv[1:8]
    passes = int(sys.argv[8]) if len(sys.argv) > 8 else 3
    ids = fidelity.id_range(ids_text)
    reference = fidelity.describe(fidelity.arctic_set(ref_dir, ids))
    ref_pitch, ref_phase = pitch_of(reference), phase_of(reference)
    plan_text, params_text = open(plan_in).read(), open(params_in).read()
    folder = tempfile.mkdtemp(dir="/tmp/es")
    table = render_set.texts()
    for step in range(passes + 1):
        plan_path = os.path.join(folder, "plan_%d.params" % step)
        params_path = os.path.join(folder, "voice_%d.params" % step)
        open(plan_path, "w").write(plan_text)
        open(params_path, "w").write(params_text)
        jobs = [(name, table[name], folder, params_path, os.path.abspath(organs), plan_path, {"KINEMATIC": "1"})
                for name in ids if name in table]
        with ThreadPoolExecutor(int(os.environ.get("WORKERS", "9"))) as pool:
            list(pool.map(render_set.render, jobs))
        ours = fidelity.describe(fidelity.ours_set(folder, ids, ref_dir))
        pitch_gap = pitch_of(ours) - ref_pitch
        phase = phase_of(ours)
        print("pass %d: pitch gap %s st (rms %.2f) | phase %s against %s rad" % (
            step, " ".join("%+.2f" % v for v in pitch_gap), numpy.sqrt(numpy.mean(pitch_gap ** 2)),
            " ".join("%.3f" % v for v in phase), " ".join("%.3f" % v for v in ref_phase)), flush=True)
        if step == passes:
            break
        ours_pitch = pitch_of(ours)
        median = slider_server.read_value(plan_text, "pitch.pitch_median") or 122.0
        spread = slider_server.read_value(plan_text, "pitch.pitch_spread") or 1.5
        lowest = slider_server.read_value(plan_text, "pitch.lowest_pitch") or 80.0
        ratio = (ref_pitch[3] - ref_pitch[1]) / max(ours_pitch[3] - ours_pitch[1], 0.1)
        plan_text = slider_server.with_values(plan_text, {
            # (Damped 0.6: on 23 sentences rms's median overshot +0.5 -> -0.5 st.)
            "pitch.pitch_median": float(median * 2.0 ** (-DAMPING * pitch_gap[2] / 12.0)),
            "pitch.pitch_spread": float(numpy.clip(spread * ratio ** DAMPING, 0.5, 3.0)),
            "pitch.lowest_pitch": float(numpy.clip(lowest * 2.0 ** (-DAMPING * (pitch_gap[0] - pitch_gap[2]) / 12.0),
                                                   40.0, 200.0))})
        if numpy.all(numpy.isfinite(phase)) and numpy.all(numpy.isfinite(ref_phase)):
            lag = slider_server.read_value(params_text, "kinematic.cycle_lag") or 0.0
            shape = slider_server.read_value(params_text, "kinematic.cycle_shape") or 0.0
            gain = slider_server.read_value(params_text, "vtl.aspiration_gain") or 1.0
            # The low bands' scatter rose about 0.25 rad per unit of lag
            # (0.20 at 0, 0.30 at 0.6): a damped step toward the target.
            low_gap = ref_phase[0] - phase[0]
            # (Half that step: at 2x rms's lag went 0 -> 0.24 -> 0.08.)
            lag = float(numpy.clip(lag + 1.0 * low_gap, 0.0, 1.2))
            shape = max(shape, 0.35 if lag > 0.3 else shape)
            # The top band's excess scatter shrank with the aspiration
            # (0.76 at 11.8, 0.63 at 6, 0.61 at 3): halve it per 0.2 rad, at
            # most a factor 2 a pass (on 6 sentences bdl's own top band read
            # 0.76 against 0.57 on 20, and an undamped step went 3 -> 10.8).
            # (By the 1-4 kHz bands together: slt's start, aspiration 13.4, sat
            # at 1.38 rad in 1-2 kHz against her 0.73 while its top band
            # matched; at gain 3, 1.04 and 1.21 against 0.73 and 1.09.)
            high_gap = numpy.mean(phase[1:] - ref_phase[1:])
            step_factor = float(numpy.clip(2.0 ** (-high_gap / 0.2), 0.5, 2.0))
            gain = float(numpy.clip(gain * step_factor, 0.3, 30.0))
            params_text = slider_server.with_values(params_text, {"kinematic.cycle_lag": lag, "kinematic.cycle_shape": shape,
                                                                  "vtl.aspiration_gain": gain})
            print("  cycle_lag %.3f cycle_shape %.3f aspiration_gain %.3f" % (lag, shape, gain), flush=True)
    header = "# source_calibrate.py against %s %s\n" % (os.path.basename(ref_dir), ids_text)
    open(plan_out, "w").write(header + plan_text)
    open(params_out, "w").write(header + params_text)


if __name__ == "__main__":
    main()
