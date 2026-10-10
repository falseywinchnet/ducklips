#!/usr/bin/env python3
"""Independent checks when the outlet center lies within the flat lip core."""
import csv
import json
from pathlib import Path
import subprocess

import numpy as np

from lip_outlet_viscous_reference import outlet_geometry, outlet_segments
from lip_profile_reference import vector
from planar_inertive_motion_reference import ROOT, sha


def main():
    out = ROOT / 'out/lip-outlet-core-reference'
    out.mkdir(exist_ok=True)
    rows = list(csv.DictReader((ROOT / 'out/lip-profile/nato.csv').open()))
    fixtures, expected = [], []
    for index in [1284, 1285, 1374, 1375, 1500, 1501]:
        profile = vector(rows[index])
        _, _, last, _, _, _ = outlet_geometry(profile)
        for order in [4, 8, 16, 32]:
            fixtures.append(str(order + 1) + ' ' + ' '.join(map(str, [*profile, 0., last, .03])))
            target = outlet_segments(profile, order, .03)
            expected.append((index, order, last, profile, target))
    binary = ROOT / 'build/ducklips_lip_viscous_probe'
    input_path = out / 'fixtures.input'
    input_path.write_text('\n'.join(fixtures) + '\n')
    result = subprocess.run([str(binary), '--outlet-fixture'], input=input_path.read_text(),
                            text=True, capture_output=True, check=True)
    output_path = out / 'fixtures.txt'
    output_path.write_text(result.stdout)
    maximum = 0.
    records = []
    for line, (index, order, last, profile, target) in zip(result.stdout.splitlines(), expected, strict=True):
        fields = np.fromstring(line, sep=' ')
        actual = fields[3:].reshape(order + 1, 3)
        reference = np.array([[r[k] for k in ['position', 'segment_inertance', 'half_gap']] for r in target])
        positive = reference != 0
        assert np.all(actual[~positive] == 0)
        error = float(np.max(abs(actual[positive] / reference[positive] - 1)))
        maximum = max(maximum, error)
        assert error < 2e-9
        if target[0]['branch_conductance']:
            assert abs(fields[0] / target[0]['branch_conductance'] - 1) < 2e-9
            assert abs(sum(actual[:, 1]) * fields[0] - 1) < 1e-13
        else:
            assert fields[0] == 0
        if index >= 1500:
            assert last < profile[4] and np.count_nonzero(actual[:, 1]) == 1
            # Entire path is a constant slit: its exact distributed wall
            # impedance is one local Stokes element of length last. There
            # must be no extra shoulder inertia or profile memories.
            minimum = target[0]['half_gap'] * .06
            assert abs(actual[0, 1] / (1.14 * last / minimum) - 1) < 1e-14
        records.append(dict(pose_index=index, order=order, path_end_m=last,
                            flat_core_m=float(profile[4]), active_sites=int(np.count_nonzero(actual[:, 1])),
                            maximum_relative_geometry_error=error))
    paths = [Path(__file__), ROOT / 'tools/lip_outlet_viscous_reference.py',
             ROOT / 'src/lip_viscous.cpp', ROOT / 'include/ducklips/lip_viscous.hpp',
             ROOT / 'src/slit_viscous.cpp', ROOT / 'include/ducklips/slit_viscous.hpp',
             ROOT / 'tools/slit_viscous_probe.cpp', ROOT / 'tools/lip_viscous_probe.cpp',
             ROOT / 'out/lip-profile/nato.csv', binary, input_path, output_path]
    report = dict(scope='Independent actual lip geometry and exact constant-slit outlet; '
                        'not a moving-wall fluid reference or speech quality check.',
                  fixtures=records, maximum_relative_geometry_error=maximum,
                  hashes={str(p.relative_to(ROOT)): sha(p) for p in paths})
    (ROOT / 'reports/lip_outlet_core_reference.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS', len(records), 'outlet fixtures; maximum geometry error', maximum)


if __name__ == '__main__':
    main()
