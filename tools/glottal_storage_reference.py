#!/usr/bin/env python3
"""Independent gas-energy and moving-tissue checks for two glottal chambers."""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
BULK=1.14*350**2
GAMMA=1.4


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def gas_audit(binary,out):
    rng=np.random.default_rng(137);fixtures=[]
    for i in range(180):
        v0=10**rng.uniform(-15,-5);v1=v0*10**rng.uniform(-1,1)
        r0,r1=np.exp(rng.uniform(np.log(.05),np.log(20),2)) if i<80 else 1+rng.uniform(-.01,.01,2)
        fixtures.append([v0*r0,v0,v1*r1,v1])
    for offset in [0.,1e-13,1e-9,1e-5,-1e-13,-1e-9,-1e-5]:
        for increment in [0.,1e-13,-1e-13,1e-5]:
            fixtures.append([1e-8*(1+offset),1e-8,3e-9*(1+offset+increment),3e-9])
    fixtures.extend([[1e-109,1e-9,1e-107,1e-8],[1e-8,1e-8,1e-12,1e-12]])
    invalid=[[0.,1e-8,1e-8,1e-8],[1e-8,0.,1e-8,1e-8],[-1e-8,1e-8,1e-8,1e-8]]
    path=out/'gas-input.txt';path.write_text(''.join(' '.join(map(str,v))+'\n' for v in fixtures+invalid))
    r=subprocess.run([str(binary),'--gas'],input=path.read_text(),capture_output=True,text=True);assert r.returncode==0,r.stderr
    output=out/'gas-output.jsonl';output.write_text(r.stdout);native=[json.loads(line) for line in r.stdout.splitlines()]
    assert len(native)==len(fixtures)+len(invalid) and all(not v['accepted'] for v in native[len(fixtures):])
    mp.mp.dps=120;k=mp.mpf(BULK);gamma=mp.mpf(GAMMA)
    def f(r):return k/(gamma*(gamma-1))*(mp.power(r,gamma)-gamma*r+gamma-1)
    def derivative(r):return k/(gamma-1)*(mp.power(r,gamma-1)-1)
    maxima=dict(energy_scaled=0.,pressure_scaled=0.,mass_effort_scaled=0.,pressure_slope_scaled=0.,chain_residual_scaled=0.)
    for values,actual in zip(fixtures,native):
        assert actual['accepted'],values
        m0,v0,m1,v1=map(mp.mpf,values);r0=m0/v0;r1=m1/v1
        mu=(f(r1)-f(r0))/(r1-r0) if r1!=r0 else derivative(r0)
        pressure=(r0+r1)*mu/2-(f(r0)+f(r1))/2
        slope_mu=(derivative(r1)-mu)/(r1-r0) if r1!=r0 else k*mp.power(r0,gamma-2)/2
        slope=-r0*r1/v1*slope_mu
        energy0=v0*f(r0);energy1=v1*f(r1)
        density_scale=max(mp.mpf(1),r0,r1)**gamma
        maxima['energy_scaled']=max(maxima['energy_scaled'],float(max(abs(mp.mpf(actual['energy0'])-energy0)/(k*v0*density_scale),abs(mp.mpf(actual['energy1'])-energy1)/(k*v1*density_scale))))
        maxima['pressure_scaled']=max(maxima['pressure_scaled'],float(abs(mp.mpf(actual['pressure'])-pressure)/(k*density_scale)))
        maxima['mass_effort_scaled']=max(maxima['mass_effort_scaled'],float(abs(mp.mpf(actual['mass_effort'])-mu)/(k*density_scale)))
        maxima['pressure_slope_scaled']=max(maxima['pressure_slope_scaled'],float(abs(mp.mpf(actual['pressure_volume_slope'])-slope)*v1/(k*density_scale)))
        residual=mp.mpf(actual['energy1'])-actual['energy0']-mp.mpf(actual['mass_effort'])*(m1-m0)+mp.mpf(actual['pressure'])*(v1-v0)
        maxima['chain_residual_scaled']=max(maxima['chain_residual_scaled'],float(abs(residual)/(k*max(v0,v1)*density_scale)))
    assert max(maxima.values())<1e-11,maxima
    return dict(valid_cases=len(fixtures),invalid_cases=len(invalid),mpmath_digits=mp.mp.dps,maxima=maxima,
                input=str(path.relative_to(ROOT)),input_sha256=sha(path),output=str(output.relative_to(ROOT)),output_sha256=sha(output))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_glottal_storage_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/glottal-storage')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/glottal_storage_reference.json')
    args=parser.parse_args();args.binary=args.binary.resolve();args.output=args.output.resolve();args.report=args.report.resolve();args.output.mkdir(parents=True,exist_ok=True)
    gas=gas_audit(args.binary,args.output);print(json.dumps(gas),flush=True)
    runs=[];oracles=[];length=.014;thickness=np.array([.0025,.0005]);wall=length*thickness;duration=.005
    for mode in ['sealed','filling','rest_motion','asymmetric','throughflow']:
        q=np.array([1.,1.15 if mode=='asymmetric' else 1.]);mass=np.r_[1.25e-4/q,2.5e-5/q];stiff=np.r_[80*q,8*q];coupling=25*q
        damping=2*np.r_[[.1]*2,[.6]*2]*np.sqrt(mass*stiff)
        pressure=np.array([700.,0.]) if mode=='sealed' else np.zeros(2)
        air_mass=3e-6*thickness*np.power(1+GAMMA*pressure/BULK,1/GAMMA)
        rates=np.array([.4e-6,.4e-6]) if mode=='throughflow' else np.array([1e-6,0.]) if mode in ['filling','asymmetric'] else np.zeros(2)
        y0=np.r_[air_mass,np.zeros(8)]
        def derivative(t,y):
            stored=y[:2];x=y[2:6];v=y[6:];rest=3e-6-1e-6*t/duration if mode=='rest_motion' else 3e-6
            volume=thickness*(rest+length*np.array([x[0]+x[1],x[2]+x[3]]))
            p=BULK/GAMMA*np.expm1(GAMMA*np.log(stored/volume))
            force=np.repeat(p*wall,2)-stiff*x-damping*v
            if mode=='throughflow':force+=np.array([1e-3,-2e-4,3e-4,-1e-4])
            difference=x[:2]-x[2:];force[:2]-=coupling*difference;force[2:]+=coupling*difference
            return np.r_[rates,v,force/mass]
        options=dict(method='DOP853',rtol=2e-12,atol=np.r_[[1e-23]*2,[1e-16]*4,[1e-11]*4])
        ode=solve_ivp(derivative,(0,duration),y0,max_step=1e-6,**options);assert ode.success
        refined=solve_ivp(derivative,(0,duration),y0,max_step=5e-7,**options);assert refined.success
        difference=abs(ode.y[:,-1]-refined.y[:,-1]);assert max(difference[2:6])<1e-11 and max(difference[6:])<1e-7
        oracles.append(dict(mode=mode,final_state=refined.y[:,-1].tolist(),refinement_difference=difference.tolist()))
        last=None
        for steps in [120,240,480,960,1920,3840,7680]:
            command=[str(args.binary),'--trajectory',str(steps),str(duration),mode]
            r=subprocess.run(command,capture_output=True,text=True);assert r.returncode==0,r.stderr
            path=args.output/f'{mode}-{steps}.json';path.write_text(r.stdout);native=json.loads(r.stdout);state=native['final']
            y=np.r_[state['mass'],state['position'],state['velocity']];error=abs(y-refined.y[:,-1])
            norm=float(max(error/np.r_[[1e-8]*2,[1e-4]*4,[1.]*4]))
            receipt=dict(mode=mode,steps=steps,command=command,maximum_mass_error_m3=float(max(error[:2])),maximum_position_error_m=float(max(error[2:6])),maximum_velocity_error_m_s=float(max(error[6:])),
                         normalized_error=norm,refinement_ratio=None if last is None else last/norm,
                         energy_balance_J=native['balance'],maximum_step_balance_J=native['maximum_step_balance'],minimum_volume_m3=native['minimum_volume'],
                         output=str(path.relative_to(ROOT)),output_sha256=sha(path))
            runs.append(receipt);last=norm
        print(json.dumps(runs[-7:]),flush=True)
        assert all(3.5<r['refinement_ratio']<4.5 for r in runs[-2:]),runs[-5:]
    checks=subprocess.run([str(args.binary)],capture_output=True,text=True);assert checks.returncode==0,checks.stderr
    log=args.output/'native-checks.log';log.write_text(checks.stdout+checks.stderr)
    sources=['src/gas_storage.cpp','include/ducklips/gas_storage.hpp','src/glottal_storage.cpp','include/ducklips/glottal_storage.hpp','include/ducklips/sources.hpp','include/ducklips/physics.hpp','tools/glottal_storage_probe.cpp','tools/glottal_storage_reference.py']
    args.report.write_text(json.dumps(dict(use='two positive-volume gas chambers coupled to tissue; supplied flow ports, no jet or speech path',gas=gas,runs=runs,oracles=oracles,
        native_checks=str(log.relative_to(ROOT)),native_checks_sha256=sha(log),hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),args.binary]},
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__)),indent=2)+'\n')


if __name__=='__main__':main()
