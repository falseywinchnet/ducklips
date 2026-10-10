#!/usr/bin/env python3
"""Check a fresh word's mouth frames against the current coupled adapter."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess

from actual_coupled_phonation import load_pose
from coupled_phonation import sha
from lip_profile_reference import vector
from planar_boundary_reference import check_boundaries
from planar_regions_reference import check_geometry

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--word', required=True)
    parser.add_argument('--body', type=Path, default=ROOT/'bodies/reference.body')
    args = parser.parse_args()
    out = args.directory.resolve()
    frames = list(csv.DictReader((out/'geometry.csv').open()))
    lips = list(csv.DictReader((out/'lips.csv').open()))
    airway = list(csv.DictReader((out/'airway.csv').open()))
    assert len(frames) == len(lips) == len(airway)
    binary = ROOT/'build/ducklips_planar_regions_probe'
    lip_binary = ROOT/'build/ducklips_lip_profile_probe'
    fixture = ''.join('4 '+' '.join(map(str, vector(row)))+'\n' for row in lips)
    (out/'lip-metrics.input').write_text(fixture)
    result = subprocess.run([str(lip_binary), '--metrics'], input=fixture, text=True, capture_output=True)
    (out/'lip-metrics.jsonl').write_text(result.stdout)
    (out/'lip-metrics.log').write_text(result.stderr)
    assert result.returncode == 0, result.stderr
    lip_metrics = [json.loads(s) for s in result.stdout.splitlines()]
    assert len(lip_metrics) == len(frames)
    failures, events, hashes = [], [], {}
    minimum_volume, minimum_cap = float('inf'), float('inf')
    old_closed = None
    with (out/'native-regions.jsonl').open('w') as trace:
        for frame, lip, throat, lm in zip(frames, lips, airway, lip_metrics):
            time = float(frame['time'])
            assert time == float(lip['time']) == float(throat['time'])
            path = out/f'poses/pose-{time:.3f}.json'
            p = load_pose(path)
            mouth, entry = p['boundary_edges']['mouth_end'], p['boundary_edges']['supraglottal_entry']
            lines = [f'1 2 {len(p["solids"])} {p["depth"]}',
                     ' '.join(map(str, p['x_edges'])), ' '.join(map(str, p['y_edges']))]
            for contour in [p['domain'], *p['solids']]:
                lines.append(str(len(contour)))
                lines.extend(' '.join(map(str, point)) for point in contour)
            lines.extend(['2', f'{mouth} {entry}'])
            result = subprocess.run([str(binary), '--ports'], input='\n'.join(lines)+'\n',
                                    text=True, capture_output=True)
            assert result.returncode == 0, (time, result.stderr)
            native = json.loads(result.stdout)
            check_geometry(p, native)
            check_boundaries(p, native, [mouth, entry])
            trace.write(json.dumps(dict(time=time, **native))+'\n')
            reasons = []
            if [r[0] for r in native['regions']] != [0, 1]:reasons.append('oral region ownership')
            if any(r[:4] != [0, 1, 1, 0] for r in native['links']):reasons.append('oral connection ownership')
            if {(r[0], r[1]) for r in native['boundary_ports']} != {(mouth, 0), (mouth, 1), (entry, 0)}:
                reasons.append('mouth or throat cap ownership')
            if min(lm['volume']) <= 0:reasons.append('vanishing lip storage')
            if any(float(throat[f'{field}_{i}']) <= 0 for i in range(8) for field in ['area', 'length']):
                reasons.append('vanishing throat storage')
            profile = vector(lip)
            if profile[-1] == 0 and profile[2] >= 1 and profile[3]-profile[4] > 0 and profile[3]+profile[4] < profile[0]:
                reasons.append('interior lip contact needs separate storage')
            if reasons:failures.append(dict(time=time, reasons=reasons))
            minimum_volume = min(minimum_volume, *lm['volume'], *[r[1] for r in native['regions']])
            minimum_cap = min(minimum_cap, *[r[2] for r in native['boundary_ports']])
            closed = lm['conductance'][-1] == 0
            if closed != old_closed:events.append(dict(time=time, closed=closed))
            old_closed = closed
            hashes[str(path.relative_to(ROOT))] = sha(path)
    paths = [Path(__file__), binary, lip_binary, args.body.resolve(), ROOT/'samples/nato.dict',
             *[out/name for name in ['geometry.csv', 'lips.csv', 'airway.csv', 'native-regions.jsonl',
                                     'lip-metrics.input', 'lip-metrics.jsonl', 'lip-metrics.log']],
             ROOT/'tools/tongue_reach_audit.cpp', ROOT/'tools/airspace_geometry_probe.cpp',
             ROOT/'build/ducklips_tongue_reach_audit', ROOT/'build/ducklips_airspace_geometry_probe']
    hashes.update({str(p.relative_to(ROOT)):sha(p) for p in paths})
    report = dict(word=args.word, scope='Fresh planner/body articulation; every exported 1 ms frame checked '
                  'with native geometry and independent region/cap quadrature. Endpoint preflight only; '
                  'not a continuous-time topology proof, acoustic render or intelligibility result.',
                  frames=len(frames), last_frame_s=float(frames[-1]['time']),
                  minimum_storage_m3=minimum_volume, minimum_cap_area_m2=minimum_cap,
                  lip_closure_events=events, incompatible_frames=failures, hashes=hashes)
    cue_path = out/'source-cues.json'
    if cue_path.exists():
        controls = json.loads(cue_path.read_text())
        cues = controls['cues']
        assert all(b['time_s'] >= a['time_s'] for a, b in zip(cues, cues[1:]))
        assert all(math.isfinite(r[k]) for r in cues for k in ['time_s', 'value', 'angular_rate'])
        assert all(0 <= r['time_s'] <= controls['duration_s'] for r in cues)
        report['source_cues'] = dict(count=len(cues), consumed_by_acoustic_render=False,
                                    nonzero_velum_targets=[r for r in cues if r['kind'] == 'velum_area_m2' and r['value'] != 0],
                                    initial=controls['initial'])
        report['hashes'][str(cue_path.relative_to(ROOT))] = sha(cue_path)
    (out/'preflight.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'hashes'}, indent=2))


if __name__ == '__main__':main()
