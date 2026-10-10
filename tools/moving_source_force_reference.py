#!/usr/bin/env python3
"""Independent moving-wall kinetic-force calculation in the saved gas pocket.

Integrate scalar axial kinetic energy at held physical velocity. Its negative
coordinate derivative equals the Hamiltonian derivative at held momentum.
This avoids the native inverse-matrix construction and its tissue constants.
"""
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import root

from coupled_phonation import read_trace, sha
from small_source_mass_reference import gas, flow_reference

ROOT = Path(__file__).resolve().parents[1]


def main():
    fixture = ROOT/'out/source-small-mass-4us/trace.jsonl'
    row = read_trace(fixture)[-1]
    source = row['restart']['source']
    area, effort = np.array(row['area']), np.array(row['effort'])
    # Fold length is 0.014 m in this fixed-tube fixture. Include both sides'
    # velocities; unlike the earlier small-mass fixture, walls are moving.
    velocity = np.array([*row['flow'], .014*sum(source[7:9]), .014*sum(source[9:11])])
    mass, _ = gas(area, effort, 48)
    independent = []
    for relative in [2e-5, 1e-5, 5e-6]:
        gradient = []
        for j in range(4):
            h = relative*mass[j] if j<2 else min(relative*max(1.,abs(area[j-2]/1e-6))*1e-6,.1*abs(area[j-2]))
            energy = []
            for sign in [1, -1]:
                target = mass.copy()
                trial_area = area.copy()
                if j<2:
                    target[j] += sign*h
                else:
                    trial_area[j-2] += sign*h

                def residual(u):
                    return (gas(trial_area, u, 48)[0]-target)/mass

                def jacobian(u):
                    return gas(trial_area, u, 48)[1]/mass[:, None]

                solved = root(residual, effort, jac=jacobian, tol=1e-10)
                assert max(abs(residual(solved.x))) < 1e-10
                ref = flow_reference.Reference([*trial_area, *solved.x])

                def scalar(z):
                    opening, density, _, _ = ref.field(z)
                    flux = ref.flow(z)@velocity
                    return np.array([.5*1.14*flux*flux/(density*opening)])

                energy.append(float(ref.integral(scalar, ref.lo, ref.hi, 1e-15)[0]))
            gradient.append(-(energy[0]-energy[1])/(2*h))
        independent.append(dict(relative_mass_step=relative, gradient_Pa=gradient[:2],area_gradient_Pa_m=gradient[2:]))
    assert max(abs(np.array(independent[1]['gradient_Pa'])-independent[2]['gradient_Pa'])) < 2e-6
    binary = ROOT/'build/ducklips_continuous_force_probe'
    native = []
    for mode, extra in [('inverse', []), ('total_metric', ['--total-metric-gradient']),
                        ('fluid_metric', ['--metric-gradient'])]:
        for factor in [1., .5, .25]:
            payload = ' '.join(map(str, [factor, *effort, *area, *velocity]))+'\n'
            result = subprocess.run([str(binary), '--relative-mass', *extra], input=payload,
                                    capture_output=True, text=True, check=True)
            record = dict(mode=mode, factor=factor, **json.loads(result.stdout))
            index = [1., .5, .25].index(factor)
            record['mass_gradient_error_Pa'] = (
                np.array(record['gradient'][:2])-independent[index]['gradient_Pa']).tolist()
            record['area_gradient_error_Pa_m'] = (
                np.array(record['gradient'][2:])-independent[index]['area_gradient_Pa_m']).tolist()
            native.append(record)
    fluid = [r for r in native if r['mode']=='fluid_metric']
    assert max(abs(v) for r in fluid for v in r['mass_gradient_error_Pa']) < 2e-6
    assert max(abs(v) for r in fluid for v in r['area_gradient_error_Pa_m']) < 5e-11
    assert abs(native[0]['mass_gradient_error_Pa'][0]) > .01
    paths = [Path(__file__), fixture, binary, ROOT/'tools/continuous_force_probe.cpp',
             ROOT/'tools/small_source_mass_reference.py', ROOT/'tools/contour_pressure_flow_reference.py',
             ROOT/'src/continuous_glottis.cpp']
    report = dict(scope=__doc__, area_m2=area.tolist(), effort_Pa=effort.tolist(),
                  velocity=velocity.tolist(), mass_m3=mass.tolist(), independent=independent,
                  native=native, hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})
    (ROOT/'reports/moving_source_force_reference.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASS moving-wall force; maximum fluid-metric error Pa:',
          max(abs(v) for r in fluid for v in r['mass_gradient_error_Pa']))


if __name__=='__main__':
    main()
