#!/usr/bin/env python3
"""Independent kinetic-force derivative in a small, partly closed gas pocket."""
import json
from pathlib import Path
import subprocess

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import root

import contour_pressure_flow_reference as flow_reference
from coupled_phonation import sha

ROOT = Path(__file__).resolve().parents[1]
# Select the full-span profile explicitly in the independent quadrature code.
flow_reference.C0 = 0.
flow_reference.C1 = flow_reference.H
flow_reference.ELL = flow_reference.H


def gas(area, effort, order):
    x, w = leggauss(order)
    lo = 0. if area[0] > 0 else -area[0]/(area[1]-area[0])
    hi = 1. if area[1] > 0 else area[0]/(area[0]-area[1])
    t = lo+(hi-lo)*(x+1)/2
    phi = np.array([1-t, t])
    opening = area@phi
    enthalpy = effort@phi
    density = (1+.4*enthalpy/139650.)**2.5
    weight = w*(hi-lo)*.003/2
    mass = phi@(opening*density*weight)
    compliance = (phi*(opening*density**.6*weight/139650.))@phi.T
    return mass, compliance


def main():
    out = ROOT/'out/source-closure-checkpoint'
    fixture = ROOT/'reports/coupled_phonation_long_failure.json'
    row = json.loads(fixture.read_text())['last_observation']
    area, effort = np.array(row['area']), np.array(row['effort'])
    # Isolate the gas kinetic term using the observed flows and held walls.
    # This is a force fixture, not a reconstruction of the missing tissue rates.
    velocity = np.array([*row['flow'], 0., 0.])
    mass, _ = gas(area, effort, 48)
    contact = -area[0]/(area[1]-area[0])
    occupied_effort = np.array([(1-contact)*effort[0]+contact*effort[1], effort[1]])
    occupied_pressure = 139650./1.4*np.expm1(3.5*np.log1p(.4*occupied_effort/139650.))

    def metric(target):
        def residual(u):return (gas(area, u, 48)[0]-target)/mass
        def jacobian(u):return gas(area, u, 48)[1]/mass[:, None]
        solved = root(residual, effort, jac=jacobian, tol=1e-10)
        assert np.max(abs(residual(solved.x))) < 1e-11
        ref = flow_reference.Reference([*area, *solved.x])
        return ref.matrix(1e-11)

    independent = []
    for relative in [2e-5, 1e-5]:
        derivative = []
        for j in range(2):
            h = relative*mass[j]
            plus, minus = mass.copy(), mass.copy()
            plus[j] += h;minus[j] -= h
            dg = (metric(plus)-metric(minus))/(2*h)
            derivative.append(float(-.5*velocity@dg@velocity))
        independent.append(dict(relative_mass_step=relative, gradient_Pa=derivative))
    target = np.array(independent[-1]['gradient_Pa'])
    assert np.max(abs(target-independent[0]['gradient_Pa'])) < .001
    records = []
    source = [*effort, *area, *velocity]
    binary = ROOT/'build/ducklips_continuous_force_probe'
    for mode, flags in [('floor', []), ('relative', ['--relative-mass'])]:
        for factor in [1., .25, .0625]:
            result = subprocess.run([str(binary), *flags], input=' '.join(map(str, [factor, *source]))+'\n',
                                    capture_output=True, text=True, check=True)
            measured = json.loads(result.stdout)
            records.append(dict(mode=mode, factor=factor, **measured,
                                mass_gradient_error_Pa=(np.array(measured['gradient'][:2])-target).tolist()))
    old = records[0]['mass_gradient_error_Pa'][0]
    new = records[3]['mass_gradient_error_Pa'][0]
    assert abs(old) > 1. and abs(new) < .001, (old, new)
    paths = [Path(__file__), ROOT/'tools/continuous_force_probe.cpp',
             ROOT/'tools/contour_pressure_flow_reference.py', fixture, binary]
    report = dict(scope='Independent derivative of fluid kinetic energy at fixed small gas moments; '
                        'observed gas flows, held walls, partly closed full-span contour. '
                        'This does not establish a through-contact trajectory repair.',
                  source=source, mass_m3=mass.tolist(), independent=independent, native=records,
                  occupied_length_m=float(.003*(1-contact)),
                  occupied_endpoint_effort_Pa=occupied_effort.tolist(),
                  occupied_endpoint_pressure_Pa=occupied_pressure.tolist(),
                  hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})
    (ROOT/'reports/small_source_mass_reference.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASS independent small-mass force reference', old, new)


if __name__ == '__main__':main()
