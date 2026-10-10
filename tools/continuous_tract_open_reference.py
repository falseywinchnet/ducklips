#!/usr/bin/env python3
"""Independent stiff ODE of the source and eight-cell tube before first contact."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from continuous_phonation_reference import Reference, YSCALE, VSCALE

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--order',type=int,default=96)
    parser.add_argument('--tolerance',type=float,default=1e-8)
    parser.add_argument('--duration',type=float,default=.000425)
    parser.add_argument('--name',default='open-96')
    args = parser.parse_args()
    out = ROOT/'out/continuous-tract-contact'
    out.mkdir(parents=True,exist_ok=True)
    receipt = json.loads((ROOT/'reports/continuous_tract_contact_reference.json').read_text())
    seed = receipt['seed']
    c = np.full(8,2e-4*.175/8/139650.)
    g = np.full(8,2e-4/(1.14*.175/8))
    g[-1] *= 2
    incidence = np.eye(8)-np.eye(8,k=-1)
    xscale = np.r_[np.full(8,1e-8),np.full(8,.1)]
    scale = np.r_[YSCALE,xscale]
    initial = np.r_[seed['source'],seed['charge'],seed['momentum']]/scale
    ref = Reference(args.order)
    def rhs(t,y):
        physical = y[:9]*YSCALE
        x = y[9:]*xscale
        derivative = ref.rhs(y[:9],1400.)
        state = ref.configuration(physical[2:4],physical[:2])
        load = np.array([0.,0.,-x[0]/c[0],0.,0.])
        derivative[4:] += np.linalg.solve(state['G']*VSCALE[:,None]*VSCALE[None,:],VSCALE*load)
        qdot = -incidence@(g*x[8:])
        qdot[0] += physical[6]
        pdot = incidence.T@(x[:8]/c)
        return np.r_[derivative,np.r_[qdot,pdot]/xscale]
    start = time.monotonic()
    solution = solve_ivp(rhs,(0.,args.duration),initial,method='Radau',rtol=args.tolerance,
        atol=args.tolerance*.01,max_step=2e-6,dense_output=True)
    assert solution.success,solution.message
    times = np.linspace(0.,args.duration,1361)
    states = solution.sol(times).T
    assert np.min(states[:,2:4]) > 0
    target = out/(args.name+'.npz')
    np.savez_compressed(target,times=times,scaled_states=states,scale=scale)
    paths = [Path(__file__),ROOT/'tools/continuous_phonation_reference.py',
        ROOT/'tools/contour_pressure_loss_reference.py',ROOT/'reports/continuous_tract_contact_reference.json',target]
    report = dict(scope='Independent open continuous-field source plus lossless eight-cell tube; ends before contact without flooring area',
        method='Radau',quadrature_order=args.order,rtol=args.tolerance,atol=args.tolerance*.01,
        maximum_step_s=2e-6,duration_s=args.duration,original_start_time_s=seed['time'],
        evaluations=solution.nfev,jacobians=solution.njev,accepted_steps=len(solution.t)-1,
        elapsed_s=time.monotonic()-start,minimum_area_m2=float(np.min(states[:,2:4])*1e-6),
        initial_scaled=initial.tolist(),final_scaled=states[-1].tolist(),
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (out/(args.name+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)

if __name__ == '__main__':
    main()
