#!/usr/bin/env python3
"""Independent matrix decay and nonlinear coupled-ODE checks of fitted damping."""
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, expm
from graph_continuous_glottis_reference import ROOT, continuous, YSCALE, XSCALE, EDGES, PORT, INITIAL, XINITIAL, STAGES, metrics, NODAL, SERIES

OUT = ROOT/'out/glottal-relaxation'

def native(binary, name, text):
    (OUT/(name+'.input')).write_text(text)
    result = subprocess.run([str(ROOT/binary), '--fixture'], input=text, capture_output=True, text=True)
    (OUT/(name+'.txt')).write_text(result.stdout)
    (OUT/(name+'.stderr')).write_text(result.stderr)
    assert result.returncode == 0, result.stderr
    return np.loadtxt(OUT/(name+'.txt'))

def matrices():
    rng = np.random.default_rng(92017)
    lines, cases = [], []
    for k in range(80):
        n = 7 if k % 4 else 5
        z = rng.normal(size=(n,n))
        a = z@z.T + np.eye(n)
        l = np.linalg.cholesky(a)
        u = np.linalg.qr(rng.normal(size=(n,n)))[0]
        rates = 10**rng.uniform(-2,4,n)
        r = np.linalg.solve(l.T,(u*rates)@u.T)@np.linalg.inv(l)
        r = (r+r.T)/2
        dt = 10**rng.uniform(-7,1)
        aa, rr = np.zeros((7,7)), np.zeros((7,7))
        aa[:n,:n], rr[:n,:n] = a, r
        lines.append(' '.join(map(str,[dt,*([1]*n+[0]*(7-n)),*aa.ravel(),*rr.ravel()])))
        cases.append((a,r,dt,l,n))
    output = native('build/ducklips_glottal_relaxation_probe','bias','\n'.join(lines)+'\n').reshape(-1,7,7)
    records = []
    for bias,(a,r,dt,l,n) in zip(output,cases):
        assert np.all(bias[n:,:] == 0) and np.all(bias[:,n:] == 0)
        b = bias[:n,:n]
        rates,u = eigh(l.T@r@l)
        x = dt*rates
        weight = np.empty(n)
        small = x < .01
        weight[small] = x[small]/6-x[small]**3/360+x[small]**5/15120
        weight[~small] = 1/np.tanh(x[~small]/2)-2/x[~small]
        expected = l@(u*weight)@u.T@l.T
        error = float(np.max(abs(np.linalg.solve(l,b-expected)@np.linalg.inv(l.T))))
        transfer = np.linalg.solve(np.eye(n)+dt*r@(a+b)/2,np.eye(n)-dt*r@(a-b)/2)
        decay_error = float(np.max(abs(transfer-expm(-dt*r@a))))
        spectrum = np.linalg.eigvalsh(np.linalg.solve(l,b)@np.linalg.inv(l.T))
        assert error < 1e-8 and decay_error < 1e-8
        assert spectrum.min() > -1e-12 and spectrum.max() < 1+1e-12
        records.append(dict(active=n, step=dt, scaled_bias_error=error, exponential_error=decay_error,
            smallest_bias_eigenvalue=float(spectrum.min()), largest_bias_eigenvalue=float(spectrum.max())))
    return records

def coupled():
    previous = json.loads((ROOT/'reports/graph_continuous_glottis_reference.json').read_text())
    cases = []
    for moving in [False, True]:
        path = ROOT/f'out/graph-continuous-glottis/continuous-{moving}.npz'
        if path.exists():
            assert hashlib.sha256(path.read_bytes()).hexdigest() == previous['hashes'][str(path.relative_to(ROOT))]
            for name in ['tools/continuous_phonation_reference.py','tools/graph_continuous_glottis_reference.py']:
                assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == previous['hashes'][name]
            reference = np.load(path)['fine']
            provenance = 'hash-verified independent ODE artifact from previous coupling audit'
        else:
            solution = continuous(moving,48,2e-10)
            reference = solution.sol(np.linspace(0.,.002,12801))
            provenance = 'fresh independent ODE integration'
        runs = []
        for steps in [200,400,800]:
            dt = .002/steps
            lines = [f'3 4 {steps} {PORT}',*[f'{i} {j}' for i,j in EDGES],
                ' '.join(map(str,INITIAL)),'1',' '.join(map(str,XINITIAL))]
            for k in range(steps):
                c,g = zip(*(metrics((k+s)*dt,moving) for s in STAGES))
                lines.append(' '.join(map(str,np.r_[dt,1400.,2.,c[0],c[1],g[0],g[1],NODAL,SERIES])))
            data = native('build/ducklips_graph_continuous_glottis_probe',f'{moving}-{steps}','\n'.join(lines)+'\n')
            assert data.shape == (steps,42)
            error = abs(data[:,:16]/np.r_[YSCALE,XSCALE]-reference[:,np.arange(1,steps+1)*(12800//steps)].T)
            assert abs(data[:,32]).max() < 1e-13 and abs(data[:,35]).max() < 1e-20
            assert data[:,31].min() >= 0
            runs.append(dict(steps=steps, dt=dt, maximum_trajectory_scaled_error=float(error.max()),
                endpoint_scaled_error=float(error[-1].max()), maximum_source_step_work_residual_J=float(abs(data[:,32]).max()),
                maximum_port_residual_J=float(abs(data[:,35]).max()), work_defect_J=float(data[:,31].sum())))
        ratios = [a['endpoint_scaled_error']/b['endpoint_scaled_error'] for a,b in zip(runs,runs[1:])]
        assert 3.7 < ratios[-1] < 4.3, ratios
        cases.append(dict(moving=moving, reference_provenance=provenance, native=runs, endpoint_refinement_ratios=ratios))
    return cases

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    matrix_cases = matrices()
    trajectory_cases = coupled()
    paths = [Path(__file__), ROOT/'build/ducklips_glottal_relaxation_probe', ROOT/'build/ducklips_graph_continuous_glottis_probe',
        ROOT/'src/glottal_relaxation.cpp',ROOT/'include/ducklips/glottal_relaxation.hpp',ROOT/'src/continuous_glottis.cpp',
        ROOT/'include/ducklips/continuous_glottis.hpp',ROOT/'tools/glottal_relaxation_probe.cpp',
        ROOT/'tools/graph_continuous_glottis_probe.cpp',ROOT/'tools/continuous_phonation_reference.py',
        ROOT/'tools/graph_continuous_glottis_reference.py',ROOT/'reports/graph_continuous_glottis_reference.json',
        *OUT.glob('*.input'),*OUT.glob('*.txt'),*OUT.glob('*.stderr')]
    report = dict(scope='Exact linear viscous relaxation in coupled momentum coordinates; nonlinear open-domain consistency, not a contact or speech validation',
        matrix_cases=matrix_cases, coupled_cases=trajectory_cases,
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/glottal_relaxation_reference.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS 80 independent matrix-decay cases and fixed/moving coupled endpoint refinement',flush=True)

if __name__ == '__main__':
    main()
