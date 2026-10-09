#!/usr/bin/env python3
"""Integrate each moving modal 2x2 block before inter-mode coupling.

Fourth-order two-point Magnus block reference, with quadratic input samples.
No native speed or complete speech claim. Compares against canonical physics.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from nonuniform_metric_reference import MovingNetwork,integrate
from moving_modal_reference import operator


def block_samples(first,middle,last):
    samples=[]
    for m,w in [first,middle,last]:
        n=len(w);a=np.zeros((n,2,2));a[:,0,0]=np.diag(m[:n,:n]);a[:,1,1]=np.diag(m[n:,n:])
        a[:,0,1]=-w;a[:,1,0]=w;samples.append(a)
    return samples


def propagator(samples,dt,fraction,start=0):
    first,middle,last=samples
    stages=[]
    for q in [.5-np.sqrt(3)/6,.5+np.sqrt(3)/6]:
        x=start+fraction*q
        stages.append((2*x*x-3*x+1)*first+(4*x-4*x*x)*middle+(2*x*x-x)*last)
    a,b=stages;h=dt*fraction
    log=.5*h*(a+b)-np.sqrt(3)*h*h/12*(a@b-b@a)
    trace=(log[:,0,0]+log[:,1,1])/2
    residual=log.copy();residual[:,0,0]-=trace;residual[:,1,1]-=trace
    square=residual[:,0,0]**2+residual[:,0,1]*residual[:,1,0]
    if np.max(square)>0:raise RuntimeError('Fixture left oscillatory block regime')
    angle=np.sqrt(-square);gain=np.sinc(angle/np.pi)
    result=residual*gain[:,None,None]
    result[:,0,0]+=np.cos(angle);result[:,1,1]+=np.cos(angle)
    return result*np.exp(trace)[:,None,None]


def apply(block,x):
    n=len(block)
    return np.r_[block[:,0,0]*x[:n]+block[:,0,1]*x[n:],block[:,1,0]*x[:n]+block[:,1,1]*x[n:]]


def inverse(block):
    determinant=block[:,0,0]*block[:,1,1]-block[:,0,1]*block[:,1,0]
    result=block.copy();result[:,0,0]=block[:,1,1];result[:,1,1]=block[:,0,0]
    result[:,0,1]*=-1;result[:,1,0]*=-1
    return result/determinant[:,None,None]


def step(x,first,middle,last,dt,nodes=2,corrections=1,phase_only=False):
    diagonal=block_samples(first,middle,last)
    gauss,weights=np.polynomial.legendre.leggauss(nodes)
    gauss=(gauss+1)/2;weights/=2
    # Collocation integral coefficients, computed here only for this oracle.
    tableau=np.zeros((nodes,nodes))
    for j in range(nodes):
        polynomial=np.polynomial.Polynomial.fromroots(np.delete(gauss,j))
        polynomial/=polynomial(gauss[j]);primitive=polynomial.integ()
        tableau[:,j]=primitive(gauss)-primitive(0)
    transforms=[];position=0.;accumulated=np.tile(np.eye(2),(len(first[1]),1,1))
    for c in gauss:
        accumulated=propagator(diagonal,dt,c-position,start=position)@accumulated
        transforms.append(accumulated);position=c
    whole=propagator(diagonal,dt,1-position,start=position)@accumulated
    if phase_only:
        def rotation(c):
            w0=first[1];wm=middle[1];w1=last[1]
            phase=dt*((2*c**3/3-1.5*c*c+c)*w0+(2*c*c-4*c**3/3)*wm+(2*c**3/3-.5*c*c)*w1)
            result=np.zeros((len(w0),2,2));co=np.cos(phase);si=np.sin(phase)
            result[:,0,0]=co;result[:,1,1]=co;result[:,0,1]=-si;result[:,1,0]=si
            return result
        transforms=[rotation(c) for c in gauss];whole=rotation(1.)
    m0=first[0];mm=middle[0];m2=last[0]
    off=[]
    for c in gauss:
        m=(2*c*c-3*c+1)*m0+(4*c-4*c*c)*mm+(2*c*c-c)*m2
        off.append(m if phase_only else m-np.diag(np.diag(m)))
    inverses=[inverse(u) for u in transforms]
    def action(i,y):return apply(inverses[i],off[i]@apply(transforms[i],y))
    stages=np.array([action(i,x) for i in range(nodes)])
    # Each correction adds an order in the small coupling. It does not
    # eliminate the separate oscillatory quadrature error.
    for _ in range(corrections):
        increments=tableau@stages
        stages=np.array([action(i,x+dt*increments[i]) for i in range(nodes)])
    return apply(whole,x+dt*(weights@stages))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--frequency-scale',type=float,default=3)
    p.add_argument('--phase-only',action='store_true')
    p.add_argument('--nodes',type=int,choices=[2,3,4],default=2)
    p.add_argument('--corrections',type=int,choices=[1,2,3,4],default=1);args=p.parse_args()
    network=MovingNetwork();scale=args.frequency_scale
    network.c0/=scale;network.g0*=scale;network.k0*=scale;network.omega*=scale
    n=network.n;dt=1/48000;frames=round(network.duration/dt)
    operators=[operator(network,i*dt/2) for i in range(2*frames+1)]
    reports=[]
    for label in ['high','mixed']:
        initial=np.zeros(2*n);initial[n-1]=.3;initial[-1]=.7
        if label=='mixed':initial=np.random.default_rng(7109).normal(0,.2,2*n)
        network.initial=operators[0][1]@initial
        times,reference,unused=integrate(network);x=initial.copy();trajectory=[network.initial]
        for i in range(frames):
            samples=[(entry[2],entry[3]) for entry in operators[2*i:2*i+3]]
            x=step(x,*samples,dt,nodes=args.nodes,corrections=args.corrections,phase_only=args.phase_only);trajectory.append(operators[2*i+2][1]@x)
        error=float(np.max(np.linalg.norm(np.array(trajectory)-reference,axis=1))/np.linalg.norm(network.initial))
        reports.append({'initial_condition':label,'max_relative_state_error':error})
    report={'schema':1,'scope':__doc__,'frequency_scale':scale,'phase_only':args.phase_only,'maximum_frequency_hz':float(max(np.max(o[3]) for o in operators)/(2*np.pi)),
            'coefficient_supply':'Every half audio sample; no control interpolation in this experiment','off_diagonal_integrator':{'gauss_nodes':args.nodes,'picard_corrections':args.corrections,'matrix_products':args.nodes*(args.corrections+1)},'cases':reports}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if any(r['max_relative_state_error']>1e-6 for r in reports):raise SystemExit('Diagonal transport accuracy gate failed')

if __name__=='__main__':main()
