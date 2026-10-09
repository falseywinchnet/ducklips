#!/usr/bin/env python3
"""Exact uniform-area motion reference for acoustic energy/state transport.

A lossless tube has C(t)=s(t) C0 and G(t)=s(t) G0, with G=1/M.
Its frozen modal frequencies stay fixed. In normalized energy coordinates
(a,b), charge/momentum conservation gives

    a' = -r*a/2 - omega*b
    b' =  omega*a + r*b/2,       r = d(log s)/dt
    E' = r*(b*b-a*a)/2,         E=(a*a+b*b)/2.

For exponential area motion r is constant, so A*A=-nu*nu I and its
propagator is available in closed form. This is an analytic physics fixture,
not a speech waveform or a replacement tract model. No audio is emitted.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm


def propagator(omega, rate, dt):
    a=np.array([[-rate/2,-omega],[omega,rate/2]])
    square=omega*omega-rate*rate/4
    if square>0:
        nu=np.sqrt(square)
        return np.eye(2)*np.cos(nu*dt)+np.sin(nu*dt)/nu*a
    if square<0:
        nu=np.sqrt(-square)
        return np.eye(2)*np.cosh(nu*dt)+np.sinh(nu*dt)/nu*a
    return np.eye(2)+dt*a


def run_case(ratio, frequency=1234.0, duration=.1, fs=48000):
    omega=2*np.pi*frequency;rate=np.log(ratio)/duration
    initial=np.array([.3,.7]);energy0=initial@initial/2
    count=round(duration*fs);times=np.arange(count+1)/fs
    exact=np.array([propagator(omega,rate,t)@initial for t in times])
    energies=np.sum(exact*exact,axis=1)/2
    # An independent numerical matrix exponential checks the closed form.
    a=np.array([[-rate/2,-omega],[omega,rate/2]])
    matrix_error=float(np.max(abs(expm(a*duration)-propagator(omega,rate,duration))))
    # Integrate the stated work rate with Gauss-Legendre quadrature in each
    # audio interval, independently from endpoint energy subtraction.
    nodes,weights=np.polynomial.legendre.leggauss(4)
    work=0.0
    for node,weight in zip(nodes,weights):
        values=exact[:-1]@propagator(omega,rate,(node+1)/(2*fs)).T
        power=rate*(values[:,1]**2-values[:,0]**2)/2
        work+=float(np.sum(power))*weight/(2*fs)
    rows=[]
    for block in [480,96,48,24,8,1]:
        for rule in ['charge_momentum','passive','pressure_flow']:
            x=initial.copy();trajectory=[x.copy()]
            step=propagator(omega,0,1/fs)
            # The geometry change is lumped at the end of each block. This
            # keeps total geometry motion and sample time identical.
            local_rate=rate*block/fs
            if rule=='charge_momentum':mapping=np.diag(np.exp([-local_rate/2,local_rate/2]))
            elif rule=='pressure_flow':mapping=np.diag(np.exp([local_rate/2,-local_rate/2]))
            else:mapping=np.eye(2)*np.exp(-abs(local_rate)/2)
            for sample in range(count):
                x=step@x
                if (sample+1)%block==0:x=mapping@x
                trajectory.append(x.copy())
            trajectory=np.array(trajectory)
            rows.append({'block_samples':block,'control_hz':fs/block,'rule':rule,
                         'relative_state_rms_error':float(np.linalg.norm(trajectory-exact)/np.linalg.norm(exact)),
                         'final_energy':float(x@x/2),
                         'final_energy_ratio':float((x@x/2)/energy0)})
    return {'area_ratio':ratio,'frequency_hz':frequency,'duration_seconds':duration,'sample_rate':fs,
            'log_area_rate_per_second':rate,'initial_energy':float(energy0),
            'exact_final_energy':float(energies[-1]),'integrated_wall_work':work,
            'energy_minus_initial_minus_work':float(energies[-1]-energy0-work),
            'matrix_exponential_absolute_error':matrix_error,'stepped_transports':rows}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();cases=[run_case(ratio) for ratio in [1.,4.,.25]]
    report={'schema':1,'model':'Lossless uniform-area tube mode; exponential area motion; no source',
            'scope':'Analytic transport reference only; excludes nonuniform motion, losses, closures and speech',
            'cases':cases}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(report,indent=2)+'\n')
    for case in cases:
        print('Area ratio',case['area_ratio'],'work residual',case['energy_minus_initial_minus_work'])
        for row in case['stepped_transports']:
            if row['block_samples'] in [48,1]:print(row)
        if abs(case['energy_minus_initial_minus_work'])>1e-10 or case['matrix_exponential_absolute_error']>1e-10:
            raise SystemExit('Analytic moving-energy reference failed')

if __name__=='__main__':main()
