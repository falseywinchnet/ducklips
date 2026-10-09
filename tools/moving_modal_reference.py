#!/usr/bin/env python3
"""Derive continuous modal transport from fixed canonical nodal potentials.

Restricted to gradient branch momentum in the nonuniform lossless fixture.
The eigenvalue gaps stay open here; this is not a degenerate-mode algorithm.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
from nonuniform_metric_reference import MovingNetwork,integrate


def operator(network,t):
    c,g,rc,rg=network.parameters(t)
    k=network.k0*np.sqrt(g)[None,:]/np.sqrt(c)[:,None]
    dk=.5*(-rc[:,None]+rg[None,:])*k
    values,u=eigh(k@k.T)
    # Choose signs relative to the nominal eigenvectors. Verified nonzero
    # overlaps below; no eigenvalue crossing occurs on this fixture path.
    overlap=np.sum(u*network.pressure_modes,axis=0)
    if np.min(abs(overlap))<.1:raise RuntimeError('Fixture mode orientation is ambiguous')
    u*=np.sign(overlap)[None,:]
    omega=np.sqrt(values)
    derivative=u.T@(dk@k.T+k@dk.T)@u
    gaps=values[None,:]-values[:,None]
    connection=np.zeros_like(gaps)
    off_diagonal=~np.eye(network.n,dtype=bool)
    if np.min(abs(gaps[off_diagonal]))<1e-6*np.max(values):raise RuntimeError('Near-degenerate fixture modes')
    connection[off_diagonal]=derivative[off_diagonal]/gaps[off_diagonal]
    f=.5*(u.T@(rc[:,None]*u))+connection
    # R = C^1/2 U, q=R a, eta=R^-T h; F=R^-1 R_dot.
    # a_dot=-Omega^2 h-F a; h_dot=a+F^T h; b=Omega h.
    flow=np.diag(np.diag(derivative)/(2*values))+omega[:,None]*f.T/omega[None,:]
    a=np.block([[-f,-np.diag(omega)],[np.diag(omega),flow]])
    reconstruction=np.zeros((network.n+network.m,2*network.n))
    reconstruction[:network.n,:network.n]=u
    reconstruction[network.n:,network.n:]=k.T@u/omega[None,:]
    transport=np.block([[-f,np.zeros_like(f)],[np.zeros_like(f),flow]])
    return a,reconstruction,transport,omega


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();network=MovingNetwork();times,full,reference=integrate(network)
    initial=np.zeros(2*network.n);initial[0]=.3;initial[network.n]=.7
    def rhs(t,x):return operator(network,t)[0]@x
    result=solve_ivp(rhs,[0,network.duration],initial,method='DOP853',t_eval=times,rtol=2e-11,atol=2e-13)
    if not result.success:raise RuntimeError(result.message)
    reconstructed=[];work_error=[];transport_norm=[];frequency=[]
    for t,x in zip(times,result.y.T):
        a,v,moving,omega=operator(network,t);physical=v@x;reconstructed.append(physical)
        work_error.append(abs(x@(a@x)-network.work_rate(t,physical)))
        transport_norm.append(np.linalg.norm(moving,2));frequency.append(omega)
    reconstructed=np.array(reconstructed)
    error=float(np.max(np.linalg.norm(reconstructed-full,axis=1))/np.linalg.norm(network.initial))
    network.initial=np.r_[.3*network.pressure_modes[:,-1],.7*network.flow_modes[:,-1]]
    unused,high_reference,high_report=integrate(network)
    high_initial=np.zeros(2*network.n);high_initial[network.n-1]=.3;high_initial[-1]=.7
    high=solve_ivp(rhs,[0,network.duration],high_initial,method='DOP853',t_eval=times,rtol=2e-11,atol=2e-13)
    if not high.success:raise RuntimeError(high.message)
    high_reconstructed=np.array([operator(network,t)[1]@x for t,x in zip(times,high.y.T)])
    high_error=float(np.max(np.linalg.norm(high_reconstructed-high_reference,axis=1))/np.linalg.norm(network.initial))
    report={'schema':1,'scope':'Continuous changing-mode coordinates; isolated eigenvalues; gradient momentum; lossless fixture only',
            'states':2*network.n,'max_relative_state_error':error,'highest_mode_max_relative_state_error':high_error,
            'max_instantaneous_wall_power_residual':float(max(work_error)),
            'max_geometry_transport_operator_norm_per_second':float(max(transport_norm)),
            'min_frozen_angular_frequency':float(np.min(frequency)),
            'max_frozen_angular_frequency':float(np.max(frequency)),
            'native_speed':'Not implemented or measured; this numerical oracle is not a speed claim',
            'direct_reference':reference}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if error>1e-8 or high_error>1e-8 or max(work_error)>1e-9:raise SystemExit('Continuous modal transport gate failed')

if __name__=='__main__':main()
