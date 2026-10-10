#!/usr/bin/env python3
"""Audit the native kinetic metric derivative at fixed mass and momentum."""
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
from continuous_phonation_reference import Reference, YSCALE

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'out/continuous-tract-contact'

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    receipt = json.loads((ROOT/'reports/continuous_tract_contact_reference.json').read_text())
    states = receipt['timescale_inputs']
    factors = [1., .25, .0625]
    fixture = ''.join(' '.join(map(str, [factor]+state['source']))+'\n'
        for state in states for factor in factors)
    executable = ROOT/'build/ducklips_continuous_force_probe'
    result = subprocess.run([str(executable)], input=fixture, text=True, capture_output=True, check=True)
    (OUT/'force.input').write_text(fixture)
    (OUT/'force.jsonl').write_text(result.stdout)
    native = list(map(json.loads, result.stdout.splitlines()))
    ref = Reference(192)
    records = []
    for i, state in enumerate(states):
        y = np.array(state['source'])
        config = ref.configuration(y[2:4], y[:2])
        velocity = y[4:]
        derivative = []
        for j in range(4):
            h = 1e-20*YSCALE[j]
            z = y[:4].astype(complex)
            z[j] += 1j*h
            dg = ref.configuration(z[2:4], z[:2], metric_only=True).imag/h
            derivative.append(-.5*velocity@dg@velocity)
        # Transform the fixed-velocity metric derivative into the kinetic
        # Hamiltonian covector at fixed gas mass and canonical momentum.
        gm = np.linalg.solve(config['C'].T, derivative[:2])
        ga = np.array(derivative[2:])-config['D'].T@gm
        target = np.r_[gm, ga]
        refinements = []
        for j, factor in enumerate(factors):
            row = native[len(factors)*i+j]
            error = np.array(row['gradient'])-target
            refinements.append(dict(step_factor=factor, **row, error=error.tolist()))
        # These bounds distinguish the measured tiny force errors from the
        # much larger trajectory discrepancy; they are not contact validity.
        assert np.max(abs(np.array(refinements[0]['error'])[:2])) < 3e-6
        assert np.max(abs(np.array(refinements[0]['error'])[2:])) < 5e-9
        records.append(dict(time=state['time'], source=state['source'],
            independent_gradient=target.tolist(), native=refinements))
    paths = [Path(__file__), ROOT/'tools/continuous_force_probe.cpp',
        ROOT/'tools/continuous_phonation_reference.py',
        ROOT/'reports/continuous_tract_contact_reference.json', executable,
        OUT/'force.input', OUT/'force.jsonl']
    report = dict(scope='Four open states in the symmetric tissue subspace; native fixed-mass inverse-metric difference versus independent complex-step chain rule. Not a closed/contact or asymmetric-state audit.',
        units=['Pa', 'Pa', 'J/m2', 'J/m2'], quadrature_order=192,
        complex_step_relative=1e-20, records=records,
        conclusion='The current difference rule has tiny errors at these states; this audit does not support changing it as a remedy for the long contact trajectory error.',
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/continuous_force_reference.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(cases=len(records), maximum_default_mass_error_Pa=max(max(abs(np.array(r['native'][0]['error'])[:2])) for r in records),
        maximum_default_area_error_J_m2=max(max(abs(np.array(r['native'][0]['error'])[2:])) for r in records)), indent=2))

if __name__ == '__main__':
    main()
