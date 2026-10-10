#!/usr/bin/env python3
"""Independent Decimal endpoint oracle for the occupied effort-field norm."""
from decimal import Decimal, localcontext
import json
from pathlib import Path
import random
import subprocess

from planar_continuous_reference import ROOT, sha


def reference(row):
    with localcontext() as context:
        context.prec = 80
        a, b, u, v, c, d, w, z = map(Decimal, row)
        d0, d1 = u-w, v-z
        values = []
        for left, right in [(a, b), (c, d)]:
            if left > 0:
                values.append(abs(d0))
            if right > 0:
                values.append(abs(d1))
            if left <= 0 or right <= 0:
                t = -left/(right-left)
                values.append(abs((1-t)*d0+t*d1))
        return float(max(values))


def main():
    rng = random.Random(74023)
    rows = []
    for _ in range(1200):
        areas = []
        for side in range(2):
            pair = [10**rng.uniform(-12, 2), 10**rng.uniform(-12, 2)]
            kind = rng.randrange(5)
            if kind < 2:
                pair[kind] *= -1
            elif kind < 4:
                pair[kind-2] = 0.
            areas.append(pair)
        efforts = [[rng.uniform(-1e4, 1e7) for _ in range(2)] for _ in range(2)]
        row = areas[0]+efforts[0]+areas[1]+efforts[1]
        rows.extend([row, row[4:]+row[:4]])
    # Small support at either cap, coincident states, and differently occupied
    # intervals. These are independent affine fields, not solver trajectories.
    rows += [[-99.,1.,1000.,0.,-9.,1.,0.,0.],
             [1.,-99.,0.,1000.,1.,-9.,0.,0.],
             [-1.,1e-12,1e7,2.,-1.,1e-12,1e7,2.]]
    payload = ''.join(' '.join(format(v,'.17g') for v in row)+'\n' for row in rows)
    binary = ROOT/'build/ducklips_contour_pressure_probe'
    result = subprocess.run([str(binary),'--effort-difference'],input=payload,
                            text=True,capture_output=True,check=True)
    values = list(map(float,result.stdout.splitlines()))
    assert len(values)==len(rows)
    maximum_relative = 0.
    for row, actual in zip(rows,values):
        expected = reference(row)
        error = abs(actual-expected)/max(expected,1e-280)
        assert error < 5e-14,(row,actual,expected,error)
        maximum_relative = max(maximum_relative,error)
    for i in range(0,2400,2):
        assert values[i]==values[i+1],('symmetry',i)
    out = ROOT/'out/source-effort-error'
    out.mkdir(exist_ok=True)
    (out/'norm-input.txt').write_text(payload)
    (out/'norm-output.txt').write_text(result.stdout)
    report = dict(cases=len(rows),maximum_relative_error=maximum_relative,
                  swapped_state_pairs=1200,
                  method='80-digit Decimal affine endpoint supremum over both occupied supports',
                  hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),binary,
                      ROOT/'src/contour_pressure_storage.cpp',ROOT/'include/ducklips/contour_pressure_storage.hpp',
                      out/'norm-input.txt',out/'norm-output.txt']})
    (out/'norm-reference.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
