#!/usr/bin/env python3
"""Interaction-coordinate RK4: exact diagonal acoustic phase, moving coupling.

A numerical experiment for the existing nonuniform canonical-network fixture.
No speech, damping model, eigenvalue crossing, or complete engine speed claim.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from moving_modal_reference import operator
from nonuniform_metric_reference import MovingNetwork,integrate


def rotate(x,phase):
    n=len(phase);c=np.cos(phase);s=np.sin(phase)
    return np.r_[c*x[:n]-s*x[n:],s*x[:n]+c*x[n:]]


def step(x,first,middle,last,dt):
    m0,w0=first;m1,w1=middle;m2,w2=last
    half=dt*(5*w0+8*w1-w2)/24
    whole=dt*(w0+4*w1+w2)/6
    k1=m0@x
    k2=rotate(m1@rotate(x+dt*k1/2,half),-half)
    k3=rotate(m1@rotate(x+dt*k2/2,half),-half)
    k4=rotate(m2@rotate(x+dt*k3,whole),-whole)
    return rotate(x+dt*(k1+2*k2+2*k3+k4)/6,whole)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();network=MovingNetwork();results=[]
    for highest in [False,True]:
        mode=network.n-1 if highest else 0
        network.initial=np.r_[.3*network.pressure_modes[:,mode],.7*network.flow_modes[:,mode]]
        times,reference,unused=integrate(network)
        cases=[]
        for subdivisions in [1,2,4]:
            dt=1/(48000*subdivisions);count=round(network.duration/dt)
            x=np.zeros(network.n*2);x[mode]=.3;x[network.n+mode]=.7
            trajectory=[network.initial];unused,v,m,w=operator(network,0);first=(m,w)
            for i in range(count):
                unused,unused2,m,w=operator(network,(i+.5)*dt);middle=(m,w)
                unused,v,m,w=operator(network,(i+1)*dt);last=(m,w)
                x=step(x,first,middle,last,dt);first=last
                if (i+1)%subdivisions==0:trajectory.append(v@x)
            error=float(np.max(np.linalg.norm(np.array(trajectory)-reference,axis=1))/np.linalg.norm(network.initial))
            cases.append({'substeps':subdivisions,'max_relative_state_error':error})
        results.append({'initial_mode':mode,'cases':cases})
    report={'schema':1,'scope':__doc__,'comparisons':results}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if any(r['cases'][0]['max_relative_state_error']>1e-6 for r in results):
        raise SystemExit('One-step interaction accuracy gate failed')

if __name__=='__main__':main()
