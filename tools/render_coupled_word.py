#!/usr/bin/env python3
"""Render a prepared oral utterance through the current coupled physical source.

The first driver follows lung pressure, effort, glottal rest opening and mouth
geometry. Tension is explicitly fixed at 1. Nasal flow and output gain commands
are rejected until supported; planned pitch is recorded but not applied.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

from actual_coupled_phonation import load_pose, save
from joint_accuracy_moving_reference import audit
from planar_continuous_reference import ROOT, run, sha, vector


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path)
    parser.add_argument('--body', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--start', type=float, default=.399)
    parser.add_argument('--fixed-tension', action='store_true', required=True,
                        help='acknowledge that planned pitch is not yet applied; Q=1')
    parser.add_argument('--binary', type=Path, default=ROOT/'build/ducklips_planar_inertive_motion_probe')
    args = parser.parse_args()
    inputs, out, body = args.inputs.resolve(), args.output.resolve(), args.body.resolve()
    if not np.isfinite(args.start) or args.start<0:
        raise ValueError('finite nonnegative start time required')
    if out.exists():
        raise ValueError('use a new output directory; preserve previous renders')
    preflight = json.loads((inputs/'preflight.json').read_text())
    if preflight['incompatible_frames']:
        raise ValueError('word geometry preflight failed')
    for name, expected in preflight['hashes'].items():
        if sha(ROOT/name)!=expected:
            raise ValueError(f'preflight input changed: {name}')
    if str(body.relative_to(ROOT)) not in preflight['hashes']:
        raise ValueError('body was not included in the geometry preflight')
    controls = [{k:float(v) for k,v in row.items()} for row in csv.DictReader((inputs/'source-controls.csv').open())]
    lips = list(csv.DictReader((inputs/'lips.csv').open()))
    airway = list(csv.DictReader((inputs/'airway.csv').open()))
    if len(controls)!=len(lips) or len(lips)!=len(airway) or len(controls)<2:
        raise ValueError('inconsistent control/geometry lengths')
    for i, row in enumerate(controls):
        if not all(np.isfinite(v) for v in row.values()):raise ValueError('nonfinite control')
        if row['time']!=float(lips[i]['time']) or row['time']!=float(airway[i]['time']):raise ValueError('unsynchronized input frames')
        if row['velum_area_m2']!=0.:raise ValueError('nasal airflow is not connected in this driver')
        if row['loudness']!=1.:raise ValueError('output gain control is not applied by this driver')
        if i and abs(row['time']-controls[i-1]['time']-.001)>1e-12:raise ValueError('expected one millisecond frames')
    start = next((i for i,r in enumerate(controls) if abs(r['time']-args.start)<1e-12),None)
    if start is None or start+1>=len(controls):raise ValueError('start must select an input frame with a successor')
    if controls[start]['lung_pressure_Pa']!=0.:
        raise ValueError('fresh rest-equilibrium initialization requires zero lung pressure')
    # This is a declared prepared initial pose, not a claim that preceding
    # articulation has been acoustically simulated and discarded.
    out.mkdir(parents=True)
    frozen_inputs = out/'inputs'
    frozen_inputs.mkdir()
    for name in ['geometry.csv','lips.csv','airway.csv','source-cues.json','source-controls.csv','preflight.json']:
        shutil.copy2(inputs/name,frozen_inputs/name)
    shutil.copytree(inputs/'poses',out/'poses')
    frozen_body, binary = out/'body.body', out/'probe'
    shutil.copy2(body,frozen_body)
    shutil.copy2(args.binary.resolve(),binary)
    paths = [Path(__file__),ROOT/'tools/planar_continuous_reference.py',
             ROOT/'tools/actual_coupled_phonation.py',ROOT/'tools/joint_accuracy_moving_reference.py',
             ROOT/'tools/planar_inertive_motion_probe.cpp',ROOT/'tools/tongue_reach_audit.cpp',
             ROOT/'tools/word_geometry_preflight.py',*sorted((ROOT/'src').glob('*.cpp')),
             *sorted((ROOT/'include/ducklips').glob('*.hpp'))]
    snapshot = out/'source-snapshot'
    for p in paths:
        dest = snapshot/p.relative_to(ROOT)
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,dest)
    source_hashes = {str(p.relative_to(ROOT)):sha(p) for p in paths}
    report = dict(status='running',word=preflight['word'],
        scope='Fresh coupled oral-word attempt with physical fold actuation and mouth motion; '
              'two oral pressure regions, fixed source tension Q=1, planned pitch not applied. '
              'No intelligibility or realism claim before listening.',
        base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_hashes=source_hashes,binary_sha256=sha(binary),body_sha256=sha(frozen_body),
        input_hashes={str(p.relative_to(ROOT)):sha(p) for p in [*frozen_inputs.iterdir(),*sorted((out/'poses').glob('*.json'))]},
        input_start_s=controls[start]['time'],input_end_s=controls[-1]['time'],
        target_duration_s=controls[-1]['time']-controls[start]['time'],completed_duration_s=0.,
        applied_controls=['lung_pressure_Pa * clamp(effort,0.3,2)','glottal_rest_area_m2','oral/throat/lip geometry'],
        unapplied_controls=dict(pitch_Hz_range=[min(r['pitch_Hz'] for r in controls),max(r['pitch_Hz'] for r in controls)],source_tension=1.),
        control_interpolation='linear between synchronized 1 ms exported states; sampled at actual trial midpoints',
        initialization='ambient gas and zero flow in prepared first mouth pose; folds at commanded spring-rest equilibrium',
        segments=[])
    save(out/'run.json',report)
    previous = None
    try:
        for i in range(start,len(controls)-1):
            a,b = controls[i],controls[i+1]
            p,q = [load_pose(out/f'poses/pose-{r["time"]:.3f}.json') for r in [a,b]]
            p['time'],q['time']=0.,b['time']-a['time']
            profiles=np.array([vector(lips[j]) for j in [i,i+1]])
            epi=np.array([[[float(airway[j][f'{field}_{k}']) for k in reversed(range(8))]
                           for field in ['area','length']] for j in [i,i+1]])
            pressure=[r['lung_pressure_Pa']*np.clip(r['effort'],.3,2.) for r in [a,b]]
            name=f'word-{i:04d}'
            report.update(active_segment=name,active_input_time_s=a['time'])
            save(out/'run.json',report)

            def started(pid):
                report['active_pid']=pid
                save(out/'run.json',report)

            beginning=time.monotonic()
            data,receipt=run(name,p,q,profiles,epi,100,0.,previous,binary,out,
                radiation=True,thermal_width=.03,viscous_width=.03,lip_viscous_width=.03,
                outlet_viscous=True,joint_relaxation=True,joint_fifth=True,joint_tolerance=1e-4,
                source_contact_step=1e-7,reservoir_pressure=pressure[0],reservoir_end=pressure[1],
                fold_rest_target=(a['glottal_rest_area_m2'],b['glottal_rest_area_m2']),
                body_path=frozen_body,initial_rest_equilibrium=previous is None,
                compact_memory=True,on_start=started)
            receipt.update(audit(data,1e-4))
            receipt.update(input_start_s=a['time'],input_end_s=b['time'],
                           controls=[a,b],applied_lung_pressure_Pa=pressure,
                           wall_seconds=time.monotonic()-beginning,
                           source_closures=sum(r['accuracy']['source_closures'] for r in data['steps']),
                           source_reopenings=sum(r['accuracy']['source_reopenings'] for r in data['steps']))
            for r in [data['initial'],*data['steps']]:
                s=r['continuous_source']
                if min(s['mass'])<=0:raise ValueError('nonpositive source gas')
                for j,port in [(0,0),(1,2)]:
                    if s['area'][j]<=0 and s['flow'][port]!=0:raise ValueError('flow through closed source port')
            previous=data['steps'][-1]
            report['segments'].append(receipt)
            report['completed_duration_s']+=q['time']
            report.pop('active_pid',None)
            save(out/'run.json',report)
            print(json.dumps(dict(name=name,completed_s=report['completed_duration_s'],
                                  wall_s=receipt['wall_seconds'])),flush=True)
        report['status']='complete'
    except Exception as error:
        report.update(status='failed',error=str(error))
        raise
    finally:
        save(out/'run.json',report)


if __name__=='__main__':
    main()
