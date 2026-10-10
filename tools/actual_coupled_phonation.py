#!/usr/bin/env python3
"""Pressure-driven source in a real mouth pose, then a continued body movement.

This uses the current continuous source, joint loss solve and actual contour
geometry. The two oral pressure regions remain a research discretization.
"""
import argparse
import copy
import csv
import json
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

from joint_accuracy_moving_reference import audit
from planar_continuous_reference import ROOT, run, sha, vector


def save(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def load_pose(path):
    pose = json.loads(path.read_text())
    mouth = pose['boundary_edges']['mouth_end']
    entry = pose['boundary_edges']['supraglottal_entry']
    pose.update(x_edges=[-.01, pose['domain'][mouth][0]],
                y_edges=[pose['domain'][entry][1], .09, .105])
    return pose


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'out/actual-coupled-phonation')
    parser.add_argument('--binary', type=Path, default=ROOT / 'build/ducklips_planar_inertive_motion_probe')
    parser.add_argument('--hold-segments', type=int, default=18)
    parser.add_argument('--reservoir', type=float, default=1400.)
    parser.add_argument('--continue-from', type=Path, help='completed run.json; hold its final mouth without resetting state')
    args = parser.parse_args()
    out = args.output.resolve()
    if args.hold_segments < 0 or not np.isfinite(args.reservoir):
        raise ValueError('invalid run settings')
    if (out / 'run.json').exists():
        raise ValueError('run already exists; preserve its state and provenance')
    out.mkdir(parents=True, exist_ok=True)
    binary = out / 'probe'
    seed = None
    if args.continue_from:
        seed_path = args.continue_from.resolve()
        seed = json.loads(seed_path.read_text())
        if seed['status'] != 'complete' or not seed['segments'] or seed['reservoir_Pa'] != args.reservoir:
            raise ValueError('continuation requires a completed run at the same lung pressure')
        source_binary = seed_path.parent / 'probe'
        if sha(source_binary) != seed['binary_sha256']:
            raise ValueError('seed executable changed')
        for name in ['bodies/reference.body', 'out/lip-profile/nato.csv', 'out/nato-spatial-airway.csv']:
            if sha(ROOT / name) != seed['source_hashes'][name]:
                raise ValueError(f'continuation changes body or articulation: {name}')
        (out / 'poses').mkdir(exist_ok=True)
        for t in [1.500, 1.501]:
            pose_path = seed_path.parent / f'poses/pose-{t:.3f}.json'
            if sha(pose_path) != seed['source_hashes'][str(pose_path.relative_to(ROOT))]:
                raise ValueError('seed pose changed')
            shutil.copy2(pose_path, out / 'poses' / pose_path.name)
        shutil.copy2(source_binary, binary)
    else:
        shutil.copy2(args.binary.resolve(), binary)
    pose_paths = [out / f'poses/pose-{t:.3f}.json' for t in [1.500, 1.501]]
    poses = [load_pose(p) for p in pose_paths]
    lips_path = ROOT / 'out/lip-profile/nato.csv'
    airway_path = ROOT / 'out/nato-spatial-airway.csv'
    lip_rows = list(csv.DictReader(lips_path.open()))
    airway_rows = list(csv.DictReader(airway_path.open()))
    profiles = np.array([vector(lip_rows[i]) for i in [1500, 1501]])
    epilarynx = np.array([[[float(airway_rows[t][f'{field}_{i}']) for i in reversed(range(8))]
                           for field in ['area', 'length']] for t in [1500, 1501]])
    paths = [Path(__file__), ROOT / 'tools/planar_continuous_reference.py',
             ROOT / 'tools/joint_accuracy_moving_reference.py', ROOT / 'bodies/reference.body',
             ROOT / 'tools/planar_inertive_motion_probe.cpp', ROOT / 'tools/airspace_geometry_probe.cpp',
             ROOT / 'out/nato-airspace-geometry.csv', lips_path, airway_path, *pose_paths,
             *sorted((ROOT / 'src').glob('*.cpp')), *sorted((ROOT / 'include/ducklips').glob('*.hpp'))]
    report = dict(scope='Current coupled source at declared lung pressure, held actual mouth, '
                        'one millisecond of actual articulation, then held final mouth. '
                        'No NATO word or realism claim; two oral pressure regions.',
                  binary_sha256=sha(binary), reservoir_Pa=args.reservoir,
                  base_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  source_hashes={str(p.relative_to(ROOT)): sha(p) for p in paths},
                  source_hashes_scope='Runtime files and inputs; continued executable provenance is inherited through seed_run.',
                  status='running', completed_duration_s=seed['completed_duration_s'] if seed else 0., segments=[])
    if seed:
        report['seed_run'] = dict(path=str(seed_path), sha256=sha(seed_path))
    save(out / 'run.json', report)
    previous = None
    if seed:
        receipt = seed['segments'][-1]
        for name, expected in receipt['files'].items():
            if sha(ROOT / name) != expected:raise ValueError(f'seed endpoint changed: {name}')
        endpoint = next(ROOT / name for name in receipt['files'] if name.endswith('.json'))
        previous = json.loads(endpoint.read_text())['steps'][-1]
    # Coordinates/profiles come from the existing articulatory gesture trace;
    # only the first/last holds are added. No recorded speech drives acoustics.
    intervals = [] if seed else [(0, 0, 'initial-hold'), (0, 1, 'body-motion')]
    intervals += [(1, 1, f'final-hold-{i:03d}') for i in range(args.hold_segments)]
    try:
        for a, b, name in intervals:
            p, q = copy.deepcopy(poses[a]), copy.deepcopy(poses[b])
            p['time'], q['time'] = 0., .001
            if sha(binary) != report['binary_sha256']:
                raise ValueError('native executable changed')
            report['active_segment'] = name
            save(out / 'run.json', report)
            start = time.monotonic()
            data, receipt = run(name, p, q, profiles[[a, b]], epilarynx[[a, b]],
                                100, 0., previous, binary, out,
                                radiation=True, thermal_width=.03, viscous_width=.03,
                                lip_viscous_width=.03, lip_viscous_order=8, outlet_viscous=True,
                                joint_relaxation=True, joint_fifth=True, joint_tolerance=1e-4,
                                joint_maximum_step=1e-5, source_contact_step=1e-7,
                                reservoir_pressure=args.reservoir)
            receipt.update(audit(data, 1e-4))
            receipt['wall_seconds'] = time.monotonic() - start
            receipt['source_closures'] = sum(r['accuracy']['source_closures'] for r in data['steps'])
            receipt['source_reopenings'] = sum(r['accuracy']['source_reopenings'] for r in data['steps'])
            for row in [data['initial'], *data['steps']]:
                source = row['continuous_source']
                if min(source['mass']) <= 0:
                    raise ValueError('nonpositive source gas')
                for j, port in [(0, 0), (1, 2)]:
                    if source['area'][j] <= 0 and source['flow'][port] != 0:
                        raise ValueError('flow through sealed source port')
            previous = data['steps'][-1]
            report['completed_duration_s'] += .001
            report['segments'].append(receipt)
            save(out / 'run.json', report)
            print(json.dumps(dict(name=name, completed_s=report['completed_duration_s'],
                                  wall_s=receipt['wall_seconds'],
                                  work_residual_J=receipt['work_residual_J'],
                                  mass_residual_m3=receipt['mass_residual_m3'])), flush=True)
        report['status'] = 'complete'
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = str(error)
        raise
    finally:
        save(out / 'run.json', report)


if __name__ == '__main__':
    main()
