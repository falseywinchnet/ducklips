#!/usr/bin/env python3
"""Collect completed coupled-contact controls and a bounded phonation prefix."""
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'out/joint-contact-v2'


def main():
    hashes = {}

    def keep(path):
        path = path.resolve()
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    def read(path):
        return json.loads(keep(path).read_text())

    def rows(name):
        return [json.loads(line) for line in keep(OUT / f'{name}.jsonl').read_text().splitlines()]

    build = read(OUT / 'build.json')
    assert hashlib.sha256(keep(OUT / 'probe').read_bytes()).hexdigest() == build['binary_sha256']
    for name, expected in build['source_hashes'].items():
        assert hashlib.sha256(keep(ROOT / name).read_bytes()).hexdigest() == expected, name

    controls = []
    for name in ['lower', 'upper']:
        data = rows(name)
        assert len(data) == 201 and data[-1]['time'] == .002
        assert data[-1]['closures'] == data[-1]['reopenings'] == 1
        for row in data:
            assert row['work_residual'] < 1e-12 and row['mass_residual'] < 1e-18
            assert min(row['gas_mass']) > 0
            for j in range(2):
                if row['area'][j] <= 0:
                    assert row['flow'][0 if j == 0 else 2] == 0
        for field in ['physical_loss', 'numerical_loss']:
            assert np.min(np.diff([row[field] for row in data])) >= 0
        assert f'PASS {name}:' in keep(OUT / f'{name}.log').read_text()
        controls.append(dict(name=name, **{key: data[-1][key] for key in
                        ['time', 'closures', 'reopenings', 'substeps', 'rejected',
                         'work_residual', 'mass_residual', 'physical_loss', 'numerical_loss']}))

    for path in [ROOT / '.build/contact-final-ctest.log', OUT / 'asan-ctest.log']:
        assert '100% tests passed, 0 tests failed out of 5' in keep(path).read_text()

    onset = rows('phonation-start')
    assert len(onset) == 501 and onset[-1]['time'] == .005
    assert 'PASS phonation:' in keep(OUT / 'phonation-start.log').read_text()
    restarted = read(OUT / 'phonation-restarted-initial.json')
    before, after = onset[-1]['restart'], restarted['restart']
    assert before['version'] == after['version'] == 1
    for key in ['charge', 'momentum', 'thermal', 'viscous', 'radiation']:
        assert before[key] == after[key], key
    scales = np.array([1e-6, 1000, 1000, *([1e-4] * 4), *([.03571428571428571] * 4), *([1e-4] * 3)])
    source_error = float(max(abs(np.array(before['source']) - after['source']) / scales))
    mass_error = float(max(abs(np.array(onset[-1]['gas_mass']) - restarted['gas_mass'])))
    energy_error = abs(onset[-1]['energy'] - restarted['energy'])
    assert source_error < 1e-10 and mass_error < 1e-20 and energy_error < 1e-14

    prefix = rows('phonation-observed-prefix')
    assert prefix[0] == restarted
    assert prefix[-1]['closures'] >= 1 and prefix[-1]['reopenings'] >= 1
    assert all(row['work_residual'] < 1e-12 and row['mass_residual'] < 1e-18 for row in prefix)
    assert all(min(row['gas_mass']) > 0 for row in prefix)

    actual = read(OUT / 'actual-closing-continuation.json')
    previous = read(ROOT / 'out/joint-accuracy-moving-current/tol-0.0001/closing-1-64.json')
    for row in actual['steps']:
        receipt = row['accuracy']
        assert receipt.pop('source_closures') == 0 and receipt.pop('source_reopenings') == 0
    assert actual == previous
    keep(OUT / 'actual-closing-continuation.log')
    keep(ROOT / 'out/joint-accuracy-moving-current/tol-0.0001/closing-1-64.txt')
    keep(OUT / 'phonation-continued.input')
    for path in [Path(__file__), ROOT / 'tools/planar_inertive_motion_probe.cpp',
                 ROOT / 'tools/planar_continuous_reference.py', ROOT / 'README.md',
                 ROOT / 'RESEARCH.md', ROOT / 'reports/COUPLED_SOURCE_CONTACT.md']:
        keep(path)
    report = dict(scope='Current coupled-source contact controls, zero-cap regression checks, '
                  'physical continuation and an explicitly bounded pressure-driven prefix. '
                  'The longer continued run and finer cycles are not claimed complete. No speech rendering.',
                  build=build, free_controls=controls, release_checks=5, sanitizer_checks=5,
                  continuation=dict(graph_and_loss_states_identical=True,
                                    maximum_scaled_source_difference=source_error,
                                    gas_mass_difference_m3=mass_error, energy_difference_J=energy_error),
                  pressure_driven_prefix=dict(onset_s=.005, continued_through_s=prefix[-1]['time'],
                                              absolute_end_s=.005 + prefix[-1]['time'],
                                              closures=prefix[-1]['closures'], reopenings=prefix[-1]['reopenings'],
                                              whole_requested_continuation_complete=False),
                  actual_moving_continuation_unchanged=True, hashes=hashes)
    (ROOT / 'reports/joint_contact_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS coupled contact evidence:', len(hashes), 'verified files')


if __name__ == '__main__':
    main()
