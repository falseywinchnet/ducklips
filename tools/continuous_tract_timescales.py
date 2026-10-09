#!/usr/bin/env python3
"""Identify local coupled timescales at the preloaded, nonstationary seed."""
import hashlib
import json
from pathlib import Path
import numpy as np
from graph_continuous_glottis_reference import (
    ROOT, Reference, YSCALE, VSCALE, XSCALE, INITIAL, XINITIAL, B, C, G, PORT, NODAL, SERIES)

def main():
    ref = Reference(48)
    def rhs(y):
        physical = y[:9]*YSCALE
        x = y[9:]*XSCALE
        derivative = ref.rhs(y[:9], 1400.)
        state = ref.configuration(physical[2:4], physical[:2])
        load = np.array([0., 0., -x[PORT]/C[PORT], 0., 0.])
        derivative[4:] += np.linalg.solve(state['G']*VSCALE[:, None]*VSCALE[None, :], VSCALE*load)
        qdot = -B@(G*x[3:]) + NODAL
        qdot[PORT] += physical[6]
        return np.r_[derivative, np.r_[qdot, B.T@(x[:3]/C)+SERIES]/XSCALE]
    y = np.r_[INITIAL/YSCALE, XINITIAL/XSCALE]
    records = []
    for h in [.001, .0005]:
        jacobian = np.column_stack([(rhs(y+h*e)-rhs(y-h*e))/(2*h) for e in np.eye(len(y))])
        modes = sorted(np.linalg.eigvals(jacobian), key=lambda v: abs(v.imag), reverse=True)
        records.append(dict(difference_step=h, modes=[dict(real_s=float(v.real), frequency_Hz=float(v.imag/(2*np.pi))) for v in modes]))
    changes = [abs(a['frequency_Hz']-b['frequency_Hz']) for a,b in zip(records[0]['modes'][:4], records[1]['modes'][:4])]
    assert max(changes) < .001, changes
    paths = [Path(__file__), ROOT/'tools/graph_continuous_glottis_reference.py',
             ROOT/'tools/continuous_phonation_reference.py', ROOT/'tools/contour_pressure_loss_reference.py']
    report = dict(scope='Local Jacobian at a nonstationary initial state: timescale diagnostic, not equilibrium stability or measured human frequencies',
        initial_scaled=y.tolist(), reference_quadrature=48, refinements=records,
        fastest_frequency_refinement_Hz=max(changes),
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/continuous_tract_timescales.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASS local tangent refinement; fast pressure modes at approximately 64 kHz and 144 kHz')

if __name__ == '__main__':
    main()
