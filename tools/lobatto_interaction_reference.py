#!/usr/bin/env python3
"""Four-stage RK4 with Lobatto nodes for oscillatory geometry interaction.

Same number of matrix products as classical RK4; six-order quadrature of the
first-order-in-coupling term. The nonlinear/state-dependent method remains
fourth order. Uses exact diagonal acoustic phase from a local frequency fit.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from nonuniform_metric_reference import MovingNetwork,integrate
from moving_modal_reference import operator
from interaction_modal_reference import rotate

U=(5-np.sqrt(5))/10;V=1-U
NODES=np.array([0,U,V,1])
INVERSE=np.linalg.inv(np.array([[c**k for k in range(4)] for c in NODES]))
PHASE=np.array([[c**(k+1)/(k+1) for k in range(4)] for c in NODES])@INVERSE
A32=1/(10*U*U);A31=V-A32;A43=5*U
A42=(1/6-(5/12)*A32*U)/((1/12)*U)-A43*V/U;A41=1-A42-A43


def step(x,samples,dt):
    matrices=[s[0] for s in samples];frequency=np.array([s[1] for s in samples])
    phase=dt*(PHASE@frequency)
    def rhs(index,state):return rotate(matrices[index]@rotate(state,phase[index]),-phase[index])
    k1=rhs(0,x);k2=rhs(1,x+dt*U*k1)
    k3=rhs(2,x+dt*(A31*k1+A32*k2))
    k4=rhs(3,x+dt*(A41*k1+A42*k2+A43*k3))
    return rotate(x+dt*(k1+5*k2+5*k3+k4)/12,phase[3])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--frequency-scale',type=float,default=3);args=p.parse_args()
    network=MovingNetwork();scale=args.frequency_scale
    network.c0/=scale;network.g0*=scale;network.k0*=scale;network.omega*=scale
    n=network.n;dt=1/48000;frames=round(network.duration/dt)
    # Keep physical inputs separate from the numerical time integrator.
    samples=[];observations=[operator(network,0)[1]]
    for i in range(frames):
        row=[]
        for c in NODES:
            unused,v,m,w=operator(network,(i+c)*dt);row.append((m,w))
        samples.append(row);observations.append(v)
    cases=[]
    for label in ['high','mixed']:
        initial=np.zeros(2*n);initial[n-1]=.3;initial[-1]=.7
        if label=='mixed':initial=np.random.default_rng(7109).normal(0,.2,2*n)
        network.initial=observations[0]@initial
        times,reference,unused=integrate(network);x=initial.copy();trajectory=[network.initial]
        for i in range(frames):
            x=step(x,samples[i],dt);trajectory.append(observations[i+1]@x)
        error=float(np.max(np.linalg.norm(np.array(trajectory)-reference,axis=1))/np.linalg.norm(network.initial))
        cases.append({'initial_condition':label,'max_relative_state_error':error})
    report={'schema':1,'scope':__doc__,'frequency_scale':scale,'coefficient_supply':'Exact physical samples at all four RK nodes; no control-rate interpolation',
            'nodes':NODES.tolist(),'stage_rows':[[U],[A31,A32],[A41,A42,A43]],'phase_weights':PHASE.tolist(),'cases':cases}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if any(c['max_relative_state_error']>1e-6 for c in cases):raise SystemExit('Lobatto interaction gate failed')

if __name__=='__main__':main()
