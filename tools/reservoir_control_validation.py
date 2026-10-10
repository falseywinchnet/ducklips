#!/usr/bin/env python3
"""Check a streamed lung-pressure ramp against hold and split controls."""
import copy
import csv
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np

from actual_coupled_phonation import load_pose
from actual_coupled_audio import join
from joint_accuracy_moving_reference import audit
from planar_continuous_reference import ROOT, run, sha, vector


def main():
    out = ROOT/'out/reservoir-control'
    out.mkdir(exist_ok=True)
    binary = out/'probe'
    shutil.copy2(ROOT/'build/ducklips_planar_inertive_motion_probe', binary)
    seed_path = ROOT/'out/actual-coupled-phonation-50ms/run.json'
    seed = json.loads(seed_path.read_text())
    assert seed['status'] == 'complete'
    receipt = seed['segments'][-1]
    for name, expected in receipt['files'].items():assert sha(ROOT/name) == expected
    endpoint = next(ROOT/name for name in receipt['files'] if name.endswith('.json'))
    previous = json.loads(endpoint.read_text())['steps'][-1]
    pose_path = seed_path.parent/'poses/pose-1.501.json'
    pose = load_pose(pose_path)
    lips_path = ROOT/'out/lip-profile/nato.csv'
    airway_path = ROOT/'out/nato-spatial-airway.csv'
    lips = list(csv.DictReader(lips_path.open()))
    airway = list(csv.DictReader(airway_path.open()))
    profiles = np.array([vector(lips[1501])]*2)
    epi = np.array([[[float(airway[1501][f'{field}_{i}']) for i in reversed(range(8))]
                     for field in ['area', 'length']]]*2)
    receipts = []

    def case(name, duration=.0002, start=1400., end=1750., state=previous,
             subdivisions=20, maximum_step=1e-5, executable=binary):
        p, q = copy.deepcopy(pose), copy.deepcopy(pose)
        p['time'], q['time'] = 0., duration
        data, receipt = run(name, p, q, profiles, epi, subdivisions, 0., state,
                            executable, out, radiation=True, thermal_width=.03,
                            viscous_width=.03, lip_viscous_width=.03,
                            lip_viscous_order=8, outlet_viscous=True,
                            joint_relaxation=True, joint_fifth=True,
                            joint_tolerance=1e-4, joint_maximum_step=maximum_step,
                            source_contact_step=1e-7, reservoir_pressure=start,
                            reservoir_end=end)
        receipt['checks'] = audit(data, 1e-4)
        receipt['join'] = join(state, data['initial'])
        receipts.append(receipt)
        print(name, receipt['checks'], flush=True)
        return data

    old_binary = ROOT/'out/lip-outlet-core-reference/planar-probe'
    old = case('old-hold', end=None, executable=old_binary)
    hold = case('hold', end=None)
    assert hold == old, 'default native output changed'
    flat = case('flat-ramp', end=1400.)
    flat.pop('reservoir_control')
    assert flat == hold, 'constant endpoint changes the physical solve'
    ramp = case('ramp')
    first = case('split-first', duration=.0001, end=1575., subdivisions=10)
    second = case('split-second', duration=.0001, start=1575., state=first['steps'][-1], subdivisions=10)
    fine = case('fine', maximum_step=2.5e-6)
    finest = case('finest', maximum_step=6.25e-7)
    resolved = case('resolved', maximum_step=3.125e-7)

    def difference(a, b):
        return dict(pressure_Pa=float(np.max(abs(np.array(a['pressure'])-b['pressure']))),
                    source_area_m2=float(np.max(abs(np.array(a['continuous_source']['area'])-b['continuous_source']['area']))),
                    source_flow_m3_s=float(np.max(abs(np.array(a['continuous_source']['flow'])-b['continuous_source']['flow']))))

    split_error = difference(ramp['steps'][-1], second['steps'][-1])
    assert split_error['pressure_Pa'] < 1e-6 and split_error['source_area_m2'] < 1e-13
    refinement = [difference(ramp['steps'][-1], fine['steps'][-1]),
                  difference(fine['steps'][-1], finest['steps'][-1]),
                  difference(finest['steps'][-1], resolved['steps'][-1])]
    # Early adaptive patterns need not produce monotonically smaller global
    # endpoint differences. Check the last refinement and retain every result.
    print('successive endpoint differences', refinement, flush=True)
    assert refinement[2]['pressure_Pa'] < refinement[1]['pressure_Pa']
    assert refinement[2]['source_area_m2'] < refinement[1]['source_area_m2']
    effect = difference(ramp['steps'][-1], hold['steps'][-1])
    assert effect['source_flow_m3_s'] > 1e-6, 'pressure control has no physical effect'

    # At fixed geometry the only external work is reservoir effort times
    # transferred reference volume. Bound every observation's accepted work
    # using the independent gas law and its monotone pressure interval.
    def effort(pressure):
        bulk, gamma = 1.14*350**2, 1.4
        return bulk/(gamma-1)*np.expm1((gamma-1)/gamma*np.log1p(gamma*pressure/bulk))
    clock = 0.
    for row in ramp['steps']:
        a = row['accuracy'];end_time = clock+row['dt']
        assert a['reservoir_volume'] > 0
        lower = effort(1400+350*clock/.0002)*a['reservoir_volume']
        upper = effort(1400+350*end_time/.0002)*a['reservoir_volume']
        assert lower-1e-14 <= a['external_work'] <= upper+1e-14
        clock = end_time
    for flags in [['--reservoir-end', 'nan'], ['--reservoir-end', '1750'],
                  ['--reservoir-end', '1750', '--reservoir-end', '1800']]:
        result = subprocess.run([str(binary), *flags], input='', text=True, capture_output=True)
        assert result.returncode == 2
    paths = [Path(__file__), ROOT/'tools/planar_inertive_motion_probe.cpp',
             ROOT/'tools/planar_continuous_reference.py', binary, old_binary,
             seed_path, endpoint, pose_path, lips_path, airway_path]
    report = dict(scope='State-preserving lung-pressure control over a 0.2 ms held actual-mouth interval. '
                        'Accelerated validation ramp, not a physiological onset or speech recording.',
                  constant_pressure_native_output_unchanged=True, split_endpoint_error=split_error,
                  successive_endpoint_differences=refinement, pressure_control_effect=effect,
                  independent_reservoir_work_bounds=True, runs=receipts,
                  hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})
    (ROOT/'reports/reservoir_control_validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASS streamed pressure, default identity, split continuity, refinement and work bounds')


if __name__ == '__main__':main()
