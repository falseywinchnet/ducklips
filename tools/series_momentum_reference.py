#!/usr/bin/env python3
"""Exact source-history augmentation for a moving lossless acoustic network.

Reference only: fixed series-source locations, no branch-specific losses,
closure or speech quality claim. Generates numerical receipts, not audio.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve
from nonuniform_metric_reference import MovingNetwork


class Sources:
    def __init__(self,network):
        self.network=network;self.locations=[0,3]
        self.e=np.eye(network.m)[:,self.locations]
        self.d=solve(network.k0@network.k0.T,network.k0@self.e,assume_a='pos')
        self.residual=self.e-network.k0.T@self.d
        self.nodal=np.eye(network.n)[:,0]

    def values(self,t):
        # Bounded pressure/flow probes, not human reference recordings.
        envelope=np.sin(np.pi*t/self.network.duration)**2
        series=envelope*np.array([250*np.sin(2*np.pi*1234*t),180*np.sin(2*np.pi*2173*t+.3)])
        nodal=50*envelope*np.sin(2*np.pi*140*t)
        return nodal,series


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    reports=[]
    for moving in [False,True]:
        network=MovingNetwork()
        if not moving:network.c_shape*=0;network.g_shape*=0
        source=Sources(network);n=network.n;m=network.m;k=len(source.locations)
        times=np.arange(round(network.duration*48000)+1)/48000
        def full_rhs(t,y):
            q=y[:n];momentum=y[n:n+m]
            c,g,rc,rg=network.parameters(t);flow,series=source.values(t)
            derivative=np.r_[-network.k0@(g*momentum)+source.nodal*flow,network.k0.T@(q/c)+source.e@series]
            work=(-np.sum(rc*q*q/c)+np.sum(rg*g*momentum*momentum))/2
            power=(q/c)@(source.nodal*flow)+(g*momentum)@(source.e@series)
            return np.r_[derivative,work,power]
        full=solve_ivp(full_rhs,[0,network.duration],np.zeros(n+m+2),t_eval=times,method='DOP853',rtol=2e-11,atol=2e-14)
        if not full.success:raise RuntimeError(full.message)
        def augmented_rhs(t,y):
            q=y[:n];h=y[n:2*n];impulse=y[2*n:]
            c,g,_,_=network.parameters(t);flow,series=source.values(t)
            momentum=network.k0.T@h+source.residual@impulse
            return np.r_[-network.k0@(g*momentum)+source.nodal*flow,q/c+source.d@series,series]
        augmented=solve_ivp(augmented_rhs,[0,network.duration],np.zeros(2*n+k),t_eval=times,method='DOP853',rtol=2e-11,atol=2e-16)
        if not augmented.success:raise RuntimeError(augmented.message)
        def dropped_rhs(t,y):
            q=y[:n];h=y[n:]
            c,g,_,_=network.parameters(t);flow,series=source.values(t)
            return np.r_[-network.k0@(g*(network.k0.T@h))+source.nodal*flow,q/c+source.d@series]
        dropped=solve_ivp(dropped_rhs,[0,network.duration],np.zeros(2*n),t_eval=times,method='DOP853',rtol=2e-11,atol=2e-16)
        if not dropped.success:raise RuntimeError(dropped.message)
        exact_energy=[];augmented_energy=[];dropped_energy=[]
        for i,t in enumerate(times):
            c,g,_,_=network.parameters(t)
            q=augmented.y[:n,i];h=augmented.y[n:2*n,i];impulse=augmented.y[2*n:,i]
            canonical=np.r_[q,network.k0.T@h+source.residual@impulse]
            bad=np.r_[dropped.y[:n,i],network.k0.T@dropped.y[n:,i]]
            exact_energy.append(network.energy_from_canonical(t,full.y[:n+m,i]))
            augmented_energy.append(network.energy_from_canonical(t,canonical))
            dropped_energy.append(network.energy_from_canonical(t,bad))
        exact=np.array(exact_energy);aug=np.array(augmented_energy);bad=np.array(dropped_energy)
        peak=float(np.max(np.linalg.norm(exact,axis=1)))
        error=float(np.max(np.linalg.norm(aug-exact,axis=1))/peak)
        pressure_error=float(np.linalg.norm(bad[:,:n]-exact[:,:n])/np.linalg.norm(exact[:,:n]))
        energy=np.sum(exact*exact,axis=1)/2
        balance=float(np.max(abs(energy-full.y[-2]-full.y[-1])))
        reports.append({'moving':moving,'full_states':n+m,'augmented_states':2*n+k,'series_source_count':k,
                        'residual_rank':int(np.linalg.matrix_rank(source.residual)),
                        'nominal_divergence_of_source_residual':float(np.linalg.norm(network.k0@source.residual)),
                        'augmented_max_relative_state_error':error,'gradient_only_pressure_rms_relative_error':pressure_error,
                        'max_energy_balance_residual':balance,'maximum_energy':float(np.max(energy)),
                        'final_wall_work':float(full.y[-2,-1]),'final_source_work':float(full.y[-1,-1])})
    report={'schema':1,'scope':__doc__,'cases':reports}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if any(r['augmented_max_relative_state_error']>1e-8 or r['max_energy_balance_residual']>1e-10 for r in reports):
        raise SystemExit('Source-history reference gate failed')

if __name__=='__main__':main()
