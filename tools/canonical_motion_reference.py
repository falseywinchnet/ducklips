#!/usr/bin/env python3
"""Use canonical dual transport and reuse physical frequency derivatives.

One dense pressure matrix determines the flow matrix. Frequencies and their
rates supply a three-knot quintic Hermite interpolant; its own derivative
provides the diagonal flow connection. This is a numerical scheduling oracle,
not yet an online/native interpolator or whole speech benchmark.
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
    constraints=[]
    for position in [-1,0,1]:
        constraints.extend([[position**k for k in range(6)],[0]+[k*position**(k-1) for k in range(1,6)]])
    hermite_inverse=np.linalg.inv(np.array(constraints,dtype=float))
    cases=[]
    for order,offsets in [(3,np.arange(-1,3)),(5,np.arange(-2,4))]:
        for rate in [500,1000,2000]:
            count=round(network.duration*rate)
            knots={i:operator(network,i/rate) for i in range(int(offsets[0]),count+int(offsets[-1])+1)}
            polynomials={}
            for cell in range(count+1):
                data=[]
                for i in [cell-1,cell,cell+1]:
                    unused,unused2,m,w=knots[i]
                    derivative=w*(np.diag(m[:n,:n])+np.diag(m[n:,n:]))
                    data.extend([w,derivative/rate])
                polynomials[cell]=hermite_inverse@np.array(data)
            coefficients=[]
            for index in range(2*frames+1):
                position=index*dt*.5*rate;cell=int(np.floor(position));fraction=position-cell
                weights=np.array([np.prod([(fraction-b)/(a-b) for b in offsets if a!=b]) for a in offsets])
                pressure=sum(weight*knots[cell+int(offset)][2][:n,:n] for weight,offset in zip(weights,offsets))
                polynomial=polynomials[cell]
                omega=np.polynomial.polynomial.polyval(fraction,polynomial)
                derivative=rate*np.polynomial.polynomial.polyval(fraction,np.arange(1,6)[:,None]*polynomial[1:])
                if np.min(omega)<=0:raise RuntimeError('Nonpositive interpolated frequency')
                flow=-omega[:,None]*pressure.T/omega[None,:]+np.diag(derivative/omega)
                moving=np.block([[pressure,np.zeros((n,n))],[np.zeros((n,n)),flow]])
                coefficients.append((moving,omega))
            x=np.zeros(2*n);x[n-1]=.3;x[-1]=.7;trajectory=[network.initial]
            for i in range(frames):
                x=step(x,*coefficients[2*i:2*i+3],dt);trajectory.append(observations[i+1]@x)
            error=float(np.max(np.linalg.norm(np.array(trajectory)-reference,axis=1))/np.linalg.norm(network.initial))
            cases.append({'pressure_polynomial_degree':order,'frequency_interpolation':'quintic Hermite, three knots with exact rates',
                          'control_rate_hz':rate,'required_future_knots':int(offsets[-1]),
                          'lookahead_seconds':float(offsets[-1]/rate),'max_relative_state_error':error})
    report={'schema':1,'scope':__doc__,'initial_condition':'highest nominal mode','cases':cases}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
