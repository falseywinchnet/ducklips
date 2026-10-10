#!/usr/bin/env python3
"""Check radiation sign/power, anti-aliasing and real native continuation."""
import copy
import json
from pathlib import Path
import tempfile

import numpy as np

from coupled_phonation import ROOT, check_join, collect, power_wave, read_trace, restart_text, sha


def rejected(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError('invalid observation or restart accepted')


def main():
    # Analytic radiation port: R=4, Q-Qrad=0.01 sin(2 pi f t).
    # The power wave must be 0.02 sin(...), including its signed phase.
    times = np.arange(10001) / 100000

    def tone(hz):
        flow = .01 * np.sin(2 * np.pi * hz * times)
        return [dict(time=float(t), radiation_resistance=4., outlet_flow=float(q),
                     radiation_flow=0.) for t, q in zip(times, flow)]

    rows = tone(1000)
    _, native, audio, rate = power_wave(rows, 48000)
    target = .02 * np.sin(2 * np.pi * 1000 * np.arange(4800) / 48000)
    assert rate == 100000 and len(audio) == 4800
    assert np.max(abs(native - .02 * np.sin(2 * np.pi * 1000 * times))) < 1e-16
    passband_error = float(np.max(abs(audio[100:-100] - target[100:-100])))
    assert passband_error < 1e-6
    energy = float(np.trapezoid(native * native, times))
    assert abs(energy - .02 ** 2 * .1 / 2) < 1e-16
    _, _, aliased, _ = power_wave(tone(35000), 48000)
    rejected_ratio = float(np.sqrt(np.mean(aliased[100:-100] ** 2)) / (.02 / np.sqrt(2)))
    assert rejected_ratio < 1e-4
    irregular = copy.deepcopy(rows)
    irregular[10]['time'] += 1e-6
    rejected(lambda: power_wave(irregular, 48000))
    invalid = copy.deepcopy(rows)
    invalid[10]['radiation_resistance'] = -1
    rejected(lambda: power_wave(invalid, 48000))

    onset_path = ROOT / 'out/joint-contact-v2/phonation-start.jsonl'
    continuation_path = ROOT / 'out/coupled-phonation-driver-control-v2/segment-0000.jsonl'
    onset, continued = read_trace(onset_path), read_trace(continuation_path)
    assert len(restart_text(onset[-1]).split()) == 415
    joined, segments, joins = collect([onset_path, continuation_path])
    assert len(joined) == len(onset) + len(continued) - 1 == 503
    assert joined[-1]['time'] == .00502
    _, _, audio, rate = power_wave(joined, 48000)
    assert rate == 100000 and len(audio) == 240
    reset = copy.deepcopy(continued[0])
    reset['restart']['thermal'][0] += 1e-15
    rejected(lambda: check_join(onset[-1], reset))
    reset = copy.deepcopy(continued[0])
    reset['restart']['source'][1] += 1
    rejected(lambda: check_join(onset[-1], reset))
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'incomplete.jsonl'
        path.write_text('\n'.join(json.dumps(r) for r in continued[:-1]))
        rejected(lambda: read_trace(path))
    run_path = continuation_path.parent / 'run.json'
    run = json.loads(run_path.read_text())
    assert run['status'] == 'complete'
    assert sha(run['binary']) == run['binary_sha256']
    assert 'PASS phonation:' in continuation_path.with_suffix('.log').read_text()
    # The 20 us control overlaps the ongoing original continuation: ensure
    # preserving state also preserves the actual coupled trajectory.
    original = [json.loads(line) for line in
                (ROOT / 'out/joint-contact-v2/phonation-continued.jsonl').read_text().splitlines()[:3]]
    original_error = {}
    for key in ['effort', 'area', 'flow', 'pressure', 'outlet_flow', 'radiation_flow']:
        error = float(np.max(abs(np.array([r[key] for r in original]) -
                                 np.array([r[key] for r in continued]))))
        original_error[key] = error
    # Segment-local arithmetic moves macro times by a few floating-point ulps.
    # Require pressure agreement five orders tighter than the 0.1 Pa scale
    # implied by this trajectory's local tolerance, rather than byte parity.
    assert original_error['pressure'] < 1e-6 and original_error['effort'] < 1e-6
    report = dict(scope='Radiation export and 20 us native continuation plumbing; not voice realism.',
                  passband_max_error=passband_error, alias_rms_ratio=rejected_ratio,
                  analytic_radiation_energy_J=energy, native_continuation_joins=joins,
                  joined_native_samples=len(joined), audio_samples=len(audio),
                  original_continuation_max_differences=original_error,
                  rejects=['irregular clock', 'negative radiation resistance', 'reset thermal memory',
                           'changed source effort', 'incomplete native segment'],
                  hashes={str(p.relative_to(ROOT)): sha(p) for p in
                          [Path(__file__), ROOT / 'tools/coupled_phonation.py', onset_path,
                           continuation_path, run_path, continuation_path.with_suffix('.input'),
                           continuation_path.with_suffix('.log')]})
    (ROOT / 'reports/coupled_phonation_controls.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
