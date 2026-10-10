#!/usr/bin/env python3
"""Collect outlet correction, old regression and completed actual motion."""
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from joint_accuracy_moving_reference import ROOT, audit, energy, air, sha


def main():
    out = ROOT / 'out/lip-outlet-core-reference'
    actual = ROOT / 'out/actual-coupled-phonation-v2'
    hashes = {}

    def keep(path):
        hashes[str(path.relative_to(ROOT))] = sha(path)
        return path

    def read(path):
        return json.loads(keep(path).read_text())

    geometry = read(ROOT / 'reports/lip_outlet_core_reference.json')
    for name, expected in geometry['hashes'].items():
        assert sha(keep(ROOT / name)) == expected, name
    for name in ['fix-ctest.log', 'fix-asan-ctest.log']:
        log = keep(ROOT / 'out/actual-coupled-phonation' / name).read_text()
        assert '100% tests passed, 0 tests failed out of 5' in log
    closing = read(out / 'closing-regression.json')
    assert closing == read(ROOT / 'out/joint-contact-v2/actual-closing-continuation.json')
    keep(out / 'closing-regression.log')
    release, asan = [read(out / f'actual-open-{name}.json') for name in ['release', 'asan']]
    for name in ['release', 'asan']:
        keep(out / f'actual-open-{name}.log')
    differences = {key: float(np.max(abs(np.array(release['steps'][-1][key]) - asan['steps'][-1][key])))
                   for key in ['pressure', 'flow', 'charge', 'momentum']}
    assert differences['pressure'] < 1e-10
    short = [audit(data, 1e-4) for data in [release, asan]]
    sites = release['initial']['viscous_sites'][-9:]
    assert sites[0]['segment_inertance'] > 0
    g = release['initial']['conductance'][6]
    assert g > 0 and abs(sites[0]['segment_inertance'] * g - 1) < 1e-13
    assert all(s['segment_inertance'] == s['half_gap'] == 0 and
               s['branch_conductance'] == g and not any(s['momentum']) for s in sites[1:])

    original = ROOT / 'out/actual-coupled-phonation'
    old = read(original / 'run.json')
    assert sha(keep(original / 'probe')) == old['binary_sha256'] and old['status'] == 'failed'
    flags = ['--lip-profile', '--epilarynx', '--continuous-folds', str(ROOT / 'bodies/reference.body'),
             '--compact-geometry', '--radiation', '--thermal-slit-width', '.03',
             '--viscous-slit-width', '.03', '--lip-viscous-width', '.03', '--lip-viscous-order', '8',
             '--outlet-viscous', '--joint-relaxation', '--joint-fifth', '--joint-tolerance', '.0001',
             '--joint-max-step', '1e-5', '--source-contact-step', '1e-7']
    control_input = keep(out / 'actual-open.input').read_text()
    failed = subprocess.run([str(original / 'probe'), *flags], input=control_input,
                            text=True, capture_output=True)
    assert failed.returncode == 2 and failed.stdout == ''

    motion = [read(actual / f'{name}-100.json') for name in ['initial-hold', 'body-motion']]
    motion_checks = [audit(data, 1e-4) for data in motion]
    for key in ['volume', 'compliance', 'conductance', 'charge', 'momentum',
                'radiation_state', 'thermal_sites', 'viscous_sites']:
        assert motion[0]['steps'][-1][key] == motion[1]['initial'][key], key
    assert abs(sum(r['net_external_air_m3'] for r in motion_checks) -
               (air(motion[-1]['steps'][-1]) - air(motion[0]['initial']))) < 1e-17
    snapshot = read(actual / 'source-snapshot.json')
    for name, expected in snapshot.items():
        assert sha(keep(actual / 'source-snapshot' / name)) == expected, name
    for name in ['initial-hold', 'body-motion']:
        keep(actual / f'{name}-100.txt')
        keep(actual / f'{name}-100.log')

    continued = [json.loads(line) for line in keep(ROOT / 'out/joint-contact-v2/phonation-continued.jsonl').read_text().splitlines()]
    assert continued[-1]['time'] == .015 and 'restart' in continued[-1]
    assert continued[-1]['closures'] == continued[-1]['reopenings'] == 5
    assert all(r['work_residual'] < 1e-12 and r['mass_residual'] < 1e-18 for r in continued)
    assert 'PASS phonation:' in keep(ROOT / 'out/joint-contact-v2/phonation-continued.log').read_text()
    keep(ROOT / 'out/joint-contact-v2/completed-20ms-radiation.json')
    for path in [Path(__file__), ROOT / 'tools/actual_coupled_phonation.py',
                 ROOT / 'tools/planar_continuous_reference.py', ROOT / 'reports/OPEN_LIP_OUTLET.md',
                 out / 'planar-probe', actual / 'probe', ROOT / '.build/asan/ducklips_planar_inertive_motion_probe']:
        keep(path)
    report = dict(scope='Open-core outlet correction, completed 2 ms held/moving actual-mouth prefix, '
                        'and completed 20 ms fixed-tube reference. Longer live runs are not claimed complete.',
                  release_checks=5, sanitizer_checks=5, independent_geometry_fixtures=24,
                  old_executable_rejects_same_open_input=True, old_closing_JSON_unchanged=True,
                  current_open_build_pressure_difference_Pa=differences['pressure'],
                  current_open_short_checks=short, frozen_actual_motion_checks=motion_checks,
                  frozen_actual_source_snapshot_files=len(snapshot),
                  fixed_tube_completed_s=.020, fixed_tube_continued_closures=5,
                  fixed_tube_continued_reopenings=5, hashes=hashes)
    (ROOT / 'reports/lip_outlet_core_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS outlet correction and voiced-articulation connection:', len(hashes), 'verified files')


if __name__ == '__main__':
    main()
