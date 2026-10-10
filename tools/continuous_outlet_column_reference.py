#!/usr/bin/env python3
"""Independent continuous motion for a shared outlet-column kinetic metric."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from continuous_phonation_reference import Reference, YSCALE, VSCALE, STIFF, REST
from graph_continuous_glottis_reference import EDGES, B, C, G, NODAL, SERIES, PORT, XSCALE, INITIAL, XINITIAL, collocation

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'out/continuous-outlet-column'
DURATION = .002

def column(t, moving):
    if not moving:
        return 80., 0.
    phase = 2*np.pi*400*t
    return 80.*(1+.2*np.sin(phase)), 16.*2*np.pi*400*np.cos(phase)

class ColumnReference(Reference):
    def __init__(self, order):
        super().__init__(order)
        self.inertance = 80.
    def configuration(self, area, effort, metric_only=False):
        result = super().configuration(area, effort, metric_only)
        if metric_only:
            result[2,2] += self.inertance
        else:
            result['G'][2,2] += self.inertance
        return result

def continuous(moving, order, tolerance):
    ref = ColumnReference(order)
    def rhs(t,y):
        ref.inertance, rate = column(t,moving)
        physical = y[:9]*YSCALE
        x = y[9:16]*XSCALE
        derivative = ref.rhs(y[:9],1400.)
        state = ref.configuration(physical[2:4],physical[:2])
        # Gdot includes the prescribed column change as well as the four
        # configuration derivatives already calculated in Reference.rhs.
        force = np.array([0.,0.,-x[PORT]/C[PORT]-rate*physical[6],0.,0.])
        derivative[4:] += np.linalg.solve(state['G']*VSCALE[:,None]*VSCALE[None,:],VSCALE*force)
        qdot = -B@(G*x[3:])+NODAL
        qdot[PORT] += physical[6]
        pdot = B.T@(x[:3]/C)+SERIES
        geometry_power = -.5*physical[6]**2*rate
        return np.r_[derivative,np.r_[qdot,pdot]/XSCALE,geometry_power/1e-6]
    sol = solve_ivp(rhs,(0.,DURATION),np.r_[INITIAL/YSCALE,XINITIAL/XSCALE,0.],
        method='Radau',rtol=tolerance,atol=tolerance*.01,max_step=1e-5,dense_output=True)
    assert sol.success,sol.message
    return sol

def run(moving,n,oracle,binary):
    dt = DURATION/n
    lines = [f'3 4 {n} {PORT}',*[f'{i} {j}' for i,j in EDGES],
        ' '.join(map(str,INITIAL)),'1',str(column(0,moving)[0]),' '.join(map(str,XINITIAL))]
    for k in range(n):
        lines.append(' '.join(map(str,np.r_[dt,1400.,2.,column((k+1)*dt,moving)[0],C,C,G,G,NODAL,SERIES])))
    name = f'{"moving" if moving else "fixed"}-{n}'
    (OUT/(name+'.input')).write_text('\n'.join(lines)+'\n')
    result = subprocess.run([str(binary),'--column-fixture'],input='\n'.join(lines)+'\n',text=True,capture_output=True)
    (OUT/(name+'.txt')).write_text(result.stdout)
    (OUT/(name+'.stderr')).write_text(result.stderr)
    assert result.returncode==0,result.stderr
    data = np.loadtxt(OUT/(name+'.txt'))
    assert data.shape==(n,44),data.shape
    ledger = data[:,23:42]
    geom = data[:,42]
    assert np.max(abs(ledger[:,9])) < 1e-13
    assert np.max(abs(ledger[:,12])) < 1e-20
    assert np.array_equal(data[:,43],np.array([column((k+1)*dt,moving)[0] for k in range(n)]))
    native = np.c_[data[:,:9]/YSCALE,data[:,9:16]/XSCALE]
    times = np.arange(1,n+1)*dt
    error = native-oracle.sol(times).T[:,:16]
    ref = ColumnReference(48)
    state = ref.configuration(INITIAL[2:4],INITIAL[:2])
    p = ref.fields(ref.t,INITIAL[2:4],INITIAL[:2])[3]
    energy = INITIAL[:2]@state['mass']-np.sum(state['w']*state['a']*p)
    energy += .5*INITIAL[4:]@state['G']@INITIAL[4:]+.5*(INITIAL[2:4]-REST)@STIFF@(INITIAL[2:4]-REST)
    energy += .5*np.sum(XINITIAL[:3]**2/C)+.5*np.sum(G*XINITIAL[3:]**2)
    graph = data[:,9:16]
    energies = ledger[:,0]+.5*np.sum(graph[:,:3]**2/C,axis=1)+.5*np.sum(G*graph[:,3:]**2,axis=1)
    work = np.cumsum(ledger[:,3:5].sum(axis=1)+ledger[:,11]+geom)
    loss = np.cumsum(ledger[:,5:9].sum(axis=1))
    residual = float(np.max(abs(energies-energy+loss-work)))
    assert residual<1e-13,residual
    initial_mass = state['mass'].sum()+XINITIAL[:3].sum()
    mass = ledger[:,1:3].sum(axis=1)+graph[:,:3].sum(axis=1)
    transfer = np.cumsum(dt*(ledger[:,17]+NODAL.sum()-data[:,22]))
    mass_error = float(np.max(abs(mass-initial_mass-transfer)))
    assert mass_error<1e-19,mass_error
    # Check the tract independently with its held source flux. Column work
    # belongs to the source, not a duplicate series-pressure graph drive.
    x=XINITIAL.copy();graph_error=0.
    for k in range(n):
        x,pressure,flow=collocation(x,np.array([C,C]),np.array([G,G]),dt,ledger[k,13])
        graph_error=max(graph_error,float(np.max(abs(x[:3]-graph[k,:3])/C)))
    assert graph_error<1e-8,graph_error
    return dict(steps=n,dt_s=dt,endpoint_scaled_error=float(np.max(abs(error[-1]))),
        maximum_sampled_scaled_error=float(np.max(abs(error))),joint_work_residual_J=residual,
        mass_residual_m3=mass_error,graph_collocation_error_Pa=graph_error,
        geometry_work_J=float(geom.sum()),independent_geometry_work_J=float(oracle.y[-1,-1]*1e-6),
        geometry_work_error_J=float(abs(geom.sum()-oracle.y[-1,-1]*1e-6)),
        port_residual_J=float(np.max(abs(ledger[:,12]))))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--steps',nargs='+',type=int,default=[200,400,800,1600])
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    binary=ROOT/'build/ducklips_graph_continuous_glottis_probe'
    cases=[]
    for moving in [False,True]:
        fine=continuous(moving,48,2e-10);coarse=continuous(moving,32,2e-9)
        times=np.linspace(0,DURATION,6401)
        discrepancy=float(np.max(abs(fine.sol(times)-coarse.sol(times))))
        assert discrepancy<2e-7,discrepancy
        np.savez_compressed(OUT/f'continuous-{moving}.npz',times=times,fine=fine.sol(times),coarse=coarse.sol(times))
        runs=[]
        for n in args.steps:
            r=run(moving,n,fine,binary);runs.append(r);print(json.dumps(dict(moving=moving,**r)),flush=True)
        ratios=[a['endpoint_scaled_error']/b['endpoint_scaled_error'] for a,b in zip(runs,runs[1:])]
        assert 3.5<ratios[-1]<4.5,ratios
        if moving:
            work_ratios=[a['geometry_work_error_J']/b['geometry_work_error_J'] for a,b in zip(runs,runs[1:])]
            assert 3.5<work_ratios[-1]<4.5,work_ratios
        else:
            work_ratios=[]
            assert all(r['geometry_work_J']==0 for r in runs)
        cases.append(dict(moving_column=moving,oracle_refinement_scaled_error=discrepancy,runs=runs,endpoint_refinement_ratios=ratios,geometry_work_refinement_ratios=work_ratios))
    paths=[Path(__file__),ROOT/'tools/continuous_phonation_reference.py',ROOT/'tools/graph_continuous_glottis_reference.py',
        ROOT/'tools/graph_continuous_glottis_probe.cpp',ROOT/'src/continuous_glottis.cpp',ROOT/'include/ducklips/continuous_glottis.hpp',
        ROOT/'src/graph_continuous_glottis.cpp',ROOT/'include/ducklips/graph_continuous_glottis.hpp',binary,*sorted(OUT.glob('*.input')),*sorted(OUT.glob('*.txt')),*sorted(OUT.glob('*.stderr')),*sorted(OUT.glob('*.npz'))]
    report=dict(scope='Open symmetric continuous source and cyclic fixed tract with constant or prescribed changing outlet-column inertia. Independent time-continuous trajectory, geometric work, mass and graph coupling; not a speech or through-contact validation.',cases=cases,
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/continuous_outlet_column_reference.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
