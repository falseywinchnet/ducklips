#!/usr/bin/env python3
"""Continuous nonuniform acoustic metric motion in canonical and energy states.

This small network is an analytic/numerical fixture, not speech or an anatomic
body fit. No audio is produced. A fixed orthonormal energy subspace preserves
the projected wall-work identity without any per-block state reprojection.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh,qr,block_diag,null_space
from moving_metric_reference import propagator


class MovingNetwork:
    def __init__(self, uniform=False):
        self.n=5
        edges=[(0,1),(1,2),(2,3),(2,4),(3,-1)] + [(i,-1) for i in range(self.n)]
        self.m=len(edges)
        incidence=np.zeros((self.n,self.m))
        for j,(a,b) in enumerate(edges):
            incidence[a,j]=1
            if b>=0:incidence[b,j]=-1
        rho=1.2;speed=343.0;length=.035;area=2e-4
        self.c0=np.array([1.,1.2,.8,1.1,.6])*area*length/(rho*speed*speed)
        self.g0=np.array([1.,.9,1.1,.4,2.,.04,.045,.035,.05,.03])*area/(rho*length)
        self.k0=incidence*np.sqrt(self.g0)[None,:]/np.sqrt(self.c0)[:,None]
        self.c_shape=np.ones(self.n) if uniform else np.array([.15,-.4,-1.,-.2,.6])
        self.g_shape=np.ones(self.m) if uniform else np.array([-.1,-.7,-.5,.3,.1,.08,-.2,-.5,-.1,.3])
        self.uniform=uniform
        self.duration=.04
        values,self.pressure_modes=eigh(self.k0@self.k0.T)
        self.omega=np.sqrt(values)
        self.flow_modes=self.k0.T@self.pressure_modes/self.omega[None,:]
        self.initial=np.r_[.3*self.pressure_modes[:,0],.7*self.flow_modes[:,0]]

    def parameters(self,t):
        if self.uniform:
            position=np.log(4)*t/self.duration;velocity=np.log(4)/self.duration
        else:
            phase=2*np.pi*t/self.duration
            position=.5*(1-np.cos(phase));velocity=np.pi/self.duration*np.sin(phase)
        c=np.exp(self.c_shape*position);g=np.exp(self.g_shape*position)
        return c,g,self.c_shape*velocity,self.g_shape*velocity

    def energy_generator(self,t):
        c,g,rc,rg=self.parameters(t)
        k=self.k0*np.sqrt(g)[None,:]/np.sqrt(c)[:,None]
        return np.block([[-np.diag(rc)/2,-k],[k.T,np.diag(rg)/2]])

    def work_rate(self,t,x):
        unused,unused2,rc,rg=self.parameters(t)
        return (-rc@(x[:self.n]**2)+rg@(x[self.n:]**2))/2

    def canonical_rhs(self,t,x):
        c,g,unused,unused2=self.parameters(t)
        return np.r_[-self.k0@(g*x[self.n:]),self.k0.T@(x[:self.n]/c)]

    def energy_from_canonical(self,t,x):
        c,g,unused,unused2=self.parameters(t)
        return np.r_[x[:self.n]/np.sqrt(c),x[self.n:]*np.sqrt(g)]

    def potential_rhs(self,t,x):
        # lambda = B^T eta remains exact when branch momentum is initially a
        # gradient and has no branch-specific damping or series forcing.
        # h = sqrt(C0) eta, so sqrt(G0) lambda = K0^T h.
        c,g,unused,unused2=self.parameters(t)
        return np.r_[-self.k0@(g*(self.k0.T@x[self.n:])),x[:self.n]/c]

    def energy_from_potential(self,t,x):
        return self.energy_from_canonical(t,np.r_[x[:self.n],self.k0.T@x[self.n:]])

    def fixed_basis(self,positions):
        # All pressure modes are retained in this fixture. Flow modes span
        # only the gradients at a frozen pose; circulation lies outside them.
        snapshots=[]
        for position in positions:
            c=np.exp(self.c_shape*position);g=np.exp(self.g_shape*position)
            k=self.k0*np.sqrt(g)[None,:]/np.sqrt(c)[:,None]
            values,u=eigh(k@k.T)
            snapshots.append(k.T@u/np.sqrt(values)[None,:])
        flow,r,unused=qr(np.hstack(snapshots),mode='economic',pivoting=True)
        rank=int(np.sum(abs(np.diag(r))>abs(r[0,0])*1e-11))
        return block_diag(self.pressure_modes,flow[:,:rank])


def integrate(network,basis=None):
    size=network.n+network.m
    v=np.eye(size) if basis is None else basis
    x0=v.T@network.initial
    def rhs(t,y):
        x=v@y[:-1]
        return np.r_[v.T@(network.energy_generator(t)@x),network.work_rate(t,x)]
    times=np.arange(round(network.duration*48000)+1)/48000
    solution=solve_ivp(rhs,[0,network.duration],np.r_[x0,0.],method='DOP853',t_eval=times,rtol=2e-11,atol=2e-13)
    if not solution.success:raise RuntimeError(solution.message)
    state=solution.y[:-1].T@v.T
    energy=np.sum(state*state,axis=1)/2
    residual=energy-energy[0]-solution.y[-1]
    return times,state,{'states':v.shape[1],'rhs_evaluations':solution.nfev,
                       'max_work_balance_residual':float(abs(residual).max()),
                       'final_wall_work':float(solution.y[-1,-1]),'initial_energy':float(energy[0]),
                       'final_energy':float(energy[-1])}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();network=MovingNetwork();times,full,full_report=integrate(network)
    canonical=solve_ivp(network.canonical_rhs,[0,network.duration],network.initial,method='DOP853',
                        t_eval=times,rtol=2e-11,atol=2e-13)
    if not canonical.success:raise RuntimeError(canonical.message)
    converted=np.array([network.energy_from_canonical(t,x) for t,x in zip(times,canonical.y.T)])
    canonical_error=float(np.max(np.linalg.norm(converted-full,axis=1))/np.linalg.norm(network.initial))
    potential_initial=np.r_[network.initial[:network.n],.7*network.pressure_modes[:,0]/network.omega[0]]
    potential=solve_ivp(network.potential_rhs,[0,network.duration],potential_initial,method='DOP853',
                        t_eval=times,rtol=2e-11,atol=2e-15)
    if not potential.success:raise RuntimeError(potential.message)
    potential_energy=np.array([network.energy_from_potential(t,x) for t,x in zip(times,potential.y.T)])
    potential_error=float(np.max(np.linalg.norm(potential_energy-full,axis=1))/np.linalg.norm(network.initial))
    # A nominal energy-null component need not be canonical circulation.
    # In this fixture all branch momenta still lie in range(K0^T).
    canonical_null=null_space(network.k0)
    canonical_circulation=float(np.max(np.linalg.norm(canonical.y[network.n:].T@canonical_null,axis=1)))
    reductions=[];arrays={'seconds':times,'full_energy_state':full,'canonical_energy_state':converted,'potential_energy_state':potential_energy}
    for label,poses in [('single_pose',[0.]),('two_poses',[0.,1.]),('three_poses',[0.,.5,1.])]:
        basis=network.fixed_basis(poses);unused,reduced,r=integrate(network,basis)
        r.update({'basis':label,'pose_parameters':poses,'relative_state_rms_error':float(np.linalg.norm(reduced-full)/np.linalg.norm(full)),
                  'max_relative_state_error':float(np.max(np.linalg.norm(reduced-full,axis=1))/np.linalg.norm(network.initial))})
        reductions.append(r);arrays[label]=reduced
    null=null_space(network.k0)
    circulation=full[:,network.n:]@null
    uniform=MovingNetwork(uniform=True);ut,ux,ur=integrate(uniform)
    exact=np.array([np.r_[uniform.pressure_modes[:,0]*(propagator(uniform.omega[0],np.log(4)/uniform.duration,t)@np.array([.3,.7]))[0],
                         uniform.flow_modes[:,0]*(propagator(uniform.omega[0],np.log(4)/uniform.duration,t)@np.array([.3,.7]))[1]] for t in ut])
    uniform_error=float(np.max(np.linalg.norm(ux-exact,axis=1))/np.linalg.norm(uniform.initial))
    report={'schema':1,'scope':'Five-node, ten-branch lossless network; smooth nonuniform positive metrics; no speech or closure claim',
            'frozen_frequencies_hz':(network.omega/(2*np.pi)).tolist(),'full':full_report,
            'canonical_vs_energy_max_relative_state_error':canonical_error,
            'maximum_energy_in_nominal_frozen_circulation_subspace':float(np.max(np.sum(circulation*circulation,axis=1)/2)),
            'uniform_exact_max_relative_state_error':uniform_error,
            'canonical_potential':{'states':2*network.n,'max_relative_state_error':potential_error,
                                   'canonical_circulation_norm':canonical_circulation,
                                   'scope':'Exact gradient momentum invariant only; no arbitrary branch source or branch-specific damping'},
            'reductions':reductions}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(args.out.with_suffix('.npz'),**arrays)
    print(json.dumps(report,indent=2))
    if canonical_error>1e-8 or uniform_error>1e-8 or potential_error>1e-8 or full_report['max_work_balance_residual']>1e-9:
        raise SystemExit('Moving-network reference gate failed')

if __name__=='__main__':main()
