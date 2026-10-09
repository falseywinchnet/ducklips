#!/usr/bin/env python3
"""Measure motion-coefficient sampling error separately from acoustic stepping.

Local polynomial interpolation of continuous, isolated-mode coefficients.
Requires future control knots (bounded lookahead), not an online predictor.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from moving_modal_reference import operator
from nonuniform_metric_reference import MovingNetwork,integrate
from interaction_modal_reference import step


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    network=MovingNetwork();n=network.n
    network.initial=np.r_[.3*network.pressure_modes[:,-1],.7*network.flow_modes[:,-1]]
    times,reference,unused=integrate(network);dt=1/48000;frames=len(times)-1
    observations=[operator(network,t)[1] for t in times]
    cases=[]
    for order,offsets in [(3,np.arange(-1,3)),(5,np.arange(-2,4))]:
        for rate in [500,1000,2000]:
            count=round(network.duration*rate)
            knots={i:operator(network,i/rate) for i in range(int(offsets[0]),count+int(offsets[-1])+1)}
            coefficients=[]
            for index in range(2*frames+1):
                position=index*dt*.5*rate;cell=int(np.floor(position));fraction=position-cell
                weights=np.array([np.prod([(fraction-b)/(a-b) for b in offsets if a!=b]) for a in offsets])
                moving=sum(weight*knots[cell+int(offset)][2] for weight,offset in zip(weights,offsets))
                omega=sum(weight*knots[cell+int(offset)][3] for weight,offset in zip(weights,offsets))
                coefficients.append((moving,omega))
            x=np.zeros(2*n);x[n-1]=.3;x[-1]=.7;trajectory=[network.initial]
            for i in range(frames):
                x=step(x,*coefficients[2*i:2*i+3],dt);trajectory.append(observations[i+1]@x)
            error=float(np.max(np.linalg.norm(np.array(trajectory)-reference,axis=1))/np.linalg.norm(network.initial))
            cases.append({'polynomial_degree':order,'control_rate_hz':rate,
                          'required_future_knots':int(offsets[-1]),'lookahead_seconds':float(offsets[-1]/rate),
                          'max_relative_state_error':error})
    report={'schema':1,'scope':__doc__,'initial_condition':'highest nominal mode','cases':cases,
            'native_interpolator':'Not implemented or timed; does not establish complete speech throughput'}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
