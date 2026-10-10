#!/usr/bin/env python3
"""Collect completed coupled-accuracy evidence without rerunning trajectories."""
import hashlib
import json
from pathlib import Path

import numpy as np

from joint_accuracy_moving_reference import audit

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    verified = {}

    def verify(mapping):
        for name, expected in mapping.items():
            assert sha(ROOT / name) == expected, name
            verified[name] = expected

    def read(name):
        path = ROOT / name
        verified[name] = sha(path)
        return json.loads(path.read_text())

    fixed = read('reports/joint_accuracy_reference.json')
    verify(fixed['hashes'])
    assert len(fixed['runs']) == 3
    for run in fixed['runs']:
        verify(run['hashes'])
        assert run['maximum_accepted_local_difference'] <= run['tolerance']
        assert run['work_residual_J'] < 1e-13
        assert run['mass_residual_m3'] < 1e-19

    moving = read('reports/joint_accuracy_moving_current.json')
    verify(moving['hashes'])
    assert len(moving['runs']) == 4
    fresh_audits = []
    for run in moving['runs']:
        verify(run['files'])
        name = next(n for n in run['files'] if n.endswith('.json'))
        data = json.loads((ROOT / name).read_text())
        fresh_audits.append(dict(name=run['name'], **audit(data, run['tolerance'])))

    exploration = read('reports/joint_accuracy_moving_exploration.json')
    verify(exploration['hashes'])
    verify(exploration['input_poses'])
    assert len(exploration['runs']) == 12
    assert not exploration['provenance']['source_snapshot_complete']
    for run in exploration['runs']:
        verify(run['files'])
        name = next(n for n in run['files'] if n.endswith('.json'))
        audit(json.loads((ROOT / name).read_text()), run['tolerance'])
    for comparison in exploration['coarse_current_comparison']:
        verify(comparison['hashes'])
        names = list(comparison['hashes'])
        assert json.loads((ROOT / names[0]).read_text()) == json.loads((ROOT / names[1]).read_text())

    sanitizer = read('out/joint-accuracy-moving-asan.json')
    sanitizer_audit = audit(sanitizer, 1e-4)
    native = read('out/joint-accuracy-moving-current/tol-0.0001/closing-0-64.json')
    assert len(native['steps']) == len(sanitizer['steps'])
    pressure_difference = max(float(np.max(abs(np.array(a['pressure']) - b['pressure'])))
                              for a, b in zip(native['steps'], sanitizer['steps']))
    assert pressure_difference < 1e-5, pressure_difference
    sanitizer_log = ROOT / 'out/joint-accuracy-moving-asan.log'
    text = sanitizer_log.read_text()
    assert 'coupled accuracy control; no speech' in text
    assert 'ERROR:' not in text and 'runtime error:' not in text

    cold_start = read('out/joint-accuracy-first-asan-check.json')
    a = np.loadtxt(ROOT / 'out/joint-accuracy-first-asan.txt')
    b = np.loadtxt(ROOT / 'out/joint-accuracy-reference/tol-6.25e-06.txt')[0]
    indices = np.flatnonzero(a != b)
    assert a.shape == b.shape == (175,)
    assert all(abs(a[i] - b[i]) <= 2 * abs(np.spacing(b[i])) for i in indices)
    assert indices.tolist() == cold_start['different_columns']
    assert (ROOT / 'out/joint-accuracy-first-asan.log').read_text() == ''

    regression = read('out/joint-accuracy-regressions/receipt.json')
    verify(regression['hashes'])
    assert len(regression['cases']) == 31
    for case in regression['cases']:
        verify(case['hashes'])
        assert case['exact']

    paths = [Path(__file__), ROOT / 'src/continuous_glottis.cpp',
             ROOT / 'src/joint_accuracy.cpp', ROOT / 'include/ducklips/joint_accuracy.hpp',
             ROOT / 'include/ducklips/relaxation_graph.hpp',
             ROOT / 'tools/graph_continuous_glottis_probe.cpp',
             ROOT / '.build/asan/ducklips_graph_continuous_glottis_probe',
             ROOT / '.build/asan/ducklips_planar_inertive_motion_probe',
             ROOT / 'out/joint-accuracy-first.input',
             ROOT / 'out/joint-accuracy-first-asan.txt',
             ROOT / 'out/joint-accuracy-first-asan.log', sanitizer_log,
             ROOT / 'README.md', ROOT / 'RESEARCH.md',
             ROOT / 'reports/COUPLED_TIME_ACCURACY.md', ROOT / 'CMakeLists.txt',
             ROOT / 'out/nato-spatial-pressure.wav']
    for name in ['out/joint-accuracy-ctest.log', 'out/joint-accuracy-asan-ctest.log']:
        path = ROOT / name
        assert '100% tests passed, 0 tests failed out of 15' in path.read_text()
        paths.append(path)
    verified.update({str(p.relative_to(ROOT)): sha(p) for p in paths})
    assert verified['out/nato-spatial-pressure.wav'] == '8284f366f2df5f3496d3139f7cbbaf4a8947ad839737d6a7e40b78635083b307'
    report = dict(scope='Fresh fixed-source ODE comparisons, current moving controls, '
                  'Release and sanitizer checks, and historical regressions. No new speech.',
                  release_checks=15, sanitizer_checks=15, exact_historical_cases=31,
                  moving_audits=fresh_audits, moving_sanitizer_audit=sanitizer_audit,
                  moving_sanitizer_pressure_difference_Pa=pressure_difference,
                  exploratory_refinement=exploration['comparisons'],
                  exploratory_source_snapshot_complete=False,
                  cold_start_sanitizer=cold_start, hashes=verified)
    (ROOT / 'reports/joint_accuracy_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS coupled accuracy evidence:', len(verified), 'verified files')


if __name__ == '__main__':
    main()
