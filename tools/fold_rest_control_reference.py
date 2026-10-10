#!/usr/bin/env python3
"""Check physical rest-opening actuation in the current coupled mouth.

Held actual Papa vowel geometry, zero lung pressure, full passive memories.
Constant bias has an independent exact work integral F dot (x1-x0).
Ramped continuation checks state preservation and local-error refinement.
This is an actuator diagnostic, not a speech recording.
"""
import copy
import csv
import json
from pathlib import Path
import shutil
import time

import numpy as np

from actual_coupled_phonation import load_pose, save
from joint_accuracy_moving_reference import audit
from planar_continuous_reference import ROOT, run, sha, vector


def position(row):
    source = row['continuous_source']
    return np.array(source['lower_position'] + source['upper_position'])


def main():
    source = ROOT / 'out/papa-current-articulation'
    out = ROOT / 'out/fold-rest-control'
    out.mkdir(exist_ok=True)
    if (out / 'run.json').exists():
        raise ValueError('preserve the existing run; use a new output directory')
    binary = out / 'probe'
    shutil.copy2(ROOT / 'build/ducklips_planar_inertive_motion_probe', binary)
    paths = [Path(__file__), ROOT / 'tools/planar_continuous_reference.py',
             ROOT / 'tools/planar_inertive_motion_probe.cpp', ROOT / 'bodies/reference.body',
             source / 'lips.csv', source / 'airway.csv', source / 'poses/pose-0.900.json',
             *sorted((ROOT / 'src').glob('*.cpp')), *sorted((ROOT / 'include/ducklips').glob('*.hpp'))]
    report = dict(status='running', scope=__doc__, binary_sha256=sha(binary),
                  hashes={str(p.relative_to(ROOT)): sha(p) for p in paths}, runs=[])
    snapshot = out / 'source-snapshot'
    for p in paths:
        dest = snapshot / p.relative_to(ROOT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dest)
    save(out / 'run.json', report)
    pose = load_pose(source / 'poses/pose-0.900.json')
    p, q = copy.deepcopy(pose), copy.deepcopy(pose)
    p['time'], q['time'] = 0., 4e-5
    lips = list(csv.DictReader((source / 'lips.csv').open()))
    airway = list(csv.DictReader((source / 'airway.csv').open()))
    profiles = np.array([vector(lips[900])] * 2)
    epi = np.array([[[float(airway[900][f'{field}_{i}']) for i in reversed(range(8))]
                     for field in ['area', 'length']]] * 2)
    endpoints = []
    try:
        for tolerance in [1e-4, 2.5e-5, 6.25e-6]:
            previous = None
            for name, target in [('constant', (4e-6, 4e-6)), ('ramp', (4e-6, 3e-6))]:
                started = time.monotonic()
                data, receipt = run(f'{name}-{tolerance}', p, q, profiles, epi, 8, 0., previous,
                    binary, out, radiation=True, thermal_width=.03, viscous_width=.03,
                    lip_viscous_width=.03, outlet_viscous=True, joint_relaxation=True,
                    joint_fifth=True, joint_tolerance=tolerance, reservoir_pressure=0.,
                    source_contact_step=1e-7, fold_rest_target=target)
                receipt.update(audit(data, tolerance))
                receipt.update(tolerance=tolerance, wall_s=time.monotonic()-started,
                               control=data['fold_rest_control'])
                if name == 'constant':
                    force = np.array(data['fold_rest_control']['start_force_N'])
                    exact_work = float(force @ (position(data['steps'][-1])-position(data['initial'])))
                    recorded_work = sum(r['accuracy']['external_work'] for r in data['steps'])
                    assert abs(exact_work-recorded_work)<1e-15, (exact_work, recorded_work)
                    receipt.update(exact_force_work_J=exact_work,
                                   recorded_external_work_J=recorded_work,
                                   force_work_error_J=abs(exact_work-recorded_work))
                previous = data['steps'][-1]
                report['runs'].append(receipt)
                save(out / 'run.json', report)
                print(json.dumps({k:receipt[k] for k in ['name','wall_s','work_residual_J','mass_residual_m3']}), flush=True)
            endpoints.append(dict(tolerance=tolerance, position_m=position(previous).tolist(),
                                  pressure_Pa=previous['pressure']))
        report['refinement'] = endpoints
        report['status'] = 'complete'
    except Exception as error:
        report.update(status='failed', error=str(error))
        raise
    finally:
        save(out / 'run.json', report)


if __name__ == '__main__':
    main()
