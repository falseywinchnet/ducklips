#!/usr/bin/env python3
"""Continuous equations for pressure-driven chambers, tissue, passages and graph."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_graph_glottal_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/graph-glottal')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/graph_glottal_reference.json')
    args=parser.parse_args();args.binary=args.binary.resolve();args.output=args.output.resolve();args.report=args.report.resolve();args.output.mkdir(parents=True,exist_ok=True)
    bulk=1.14*350**2;gamma=1.4;duration=.005;length=.014;thickness=np.array([.0025,.0005]);wall=length*thickness
    tissue_mass=np.array([1.25e-4]*2+[2.5e-5]*2);stiff=np.array([80.]*2+[8.]*2);damping=2*np.array([.1]*2+[.6]*2)*np.sqrt(tissue_mass*stiff)
    inertance=np.array([150.,200.,100.]);resistance=np.array([1e6,2e6,1e6]);quadratic=np.array([6e10,4e10,6e10])
    c=np.array([2e-10,3e-10,1e-10]);g=np.array([.002,.003,.004,.002]);b=np.array([[1,0,1,1],[-1,1,0,0],[0,-1,-1,0]])
    runs=[];oracles=[]
    # Coordinates: gas mass, tissue position/velocity, passage flow,
    # graph pressure/flow. No native finite-change efforts or implicit solver.
    for mode in ['forward','reverse','relaxation','rest_motion']:
        reservoir=0. if mode=='relaxation' else -200. if mode=='reverse' else 700.
        reservoir_mu=bulk/(gamma-1)*np.expm1((gamma-1)/gamma*np.log1p(gamma*reservoir/bulk))
        pressure=np.array([700.,300.]) if mode=='relaxation' else np.zeros(2)
        m0=thickness*3e-6*np.power(1+gamma*pressure/bulk,1/gamma)
        y0=np.r_[m0,np.zeros(18)]
        if mode=='relaxation':y0[10:13]=np.array([.001,-.0005,.0004])/inertance;y0[13:16]=np.array([2e-8,-1e-8,3e-8])/c
        def derivative(t,y):
            m=y[:2];x=y[2:6];v=y[6:10];q=y[10:13];p=y[13:16];flow=y[16:]
            rest=3e-6-5e-7*t/duration if mode=='rest_motion' else 3e-6
            volume=thickness*(rest+length*np.array([x[0]+x[1],x[2]+x[3]]));log_r=np.log(m/volume)
            pressure=bulk/gamma*np.expm1(gamma*log_r);mu=bulk/(gamma-1)*np.expm1((gamma-1)*log_r)
            force=np.repeat(pressure*wall,2)-stiff*x-damping*v;coupled=25*(x[:2]-x[2:]);force[:2]-=coupled;force[2:]+=coupled
            qdot=(np.array([reservoir_mu-mu[0],mu[0]-mu[1],mu[1]-p[2]])-resistance*q-quadratic*q*abs(q))/inertance
            drive=-b@flow;drive[2]+=q[2]
            return np.r_[q[0]-q[1],q[1]-q[2],v,force/tissue_mass,qdot,drive/c,g*(b.T@p)]
        options=dict(method='DOP853',rtol=2e-12,atol=np.r_[[1e-23]*2,[1e-16]*4,[1e-10]*4,[1e-15]*3,[1e-8]*3,[1e-15]*4])
        ode=solve_ivp(derivative,(0,duration),y0,max_step=1e-6,**options);assert ode.success
        refined=solve_ivp(derivative,(0,duration),y0,max_step=5e-7,**options);assert refined.success
        difference=abs(ode.y[:,-1]-refined.y[:,-1]);scale=np.r_[[1e-8]*2,[1e-4]*4,[1.]*4,[1e-4]*3,[700.]*3,[1e-4]*4]
        assert max(difference/scale)<1e-8,difference
        oracles.append(dict(mode=mode,final_state=refined.y[:,-1].tolist(),normalized_refinement_difference=float(max(difference/scale))))
        previous=None
        for steps in [120,240,480,960,1920,3840,7680,15360]:
            command=[str(args.binary),'--trajectory',str(steps),str(duration),mode]
            r=subprocess.run(command,capture_output=True,text=True);assert r.returncode==0,r.stderr
            path=args.output/f'{mode}-{steps}.json';path.write_text(r.stdout);native=json.loads(r.stdout)
            y=np.r_[native['mass'],native['position'],native['velocity'],np.array(native['passage_momentum'])/inertance,np.array(native['charge'])/c,g*np.array(native['network_momentum'])]
            error=abs(y-refined.y[:,-1]);norm=float(max(error/scale))
            receipt=dict(mode=mode,steps=steps,normalized_error=norm,refinement_ratio=None if previous is None else previous/norm,
                maximum_mass_error_m3=float(max(error[:2])),maximum_position_error_m=float(max(error[2:6])),maximum_velocity_error_m_s=float(max(error[6:10])),
                maximum_passage_flow_error_m3_s=float(max(error[10:13])),maximum_network_pressure_error_Pa=float(max(error[13:16])),maximum_network_flow_error_m3_s=float(max(error[16:])),
                maximum_mass_balance_m3=native['maximum_mass_balance'],maximum_work_balance_J=native['maximum_work_balance'],maximum_law_residual_Pa=native['maximum_law_residual'],
                command=command,output=str(path.relative_to(ROOT)),output_sha256=sha(path))
            runs.append(receipt);previous=norm
        print(json.dumps(runs[-1]),flush=True)
        assert all(3.5<r['refinement_ratio']<4.5 for r in runs[-2:]),runs[-4:]
        assert norm<1e-4,runs[-1]
    checks=subprocess.run([str(args.binary)],capture_output=True,text=True);assert checks.returncode==0,checks.stderr
    log=args.output/'native-checks.log';log.write_text(checks.stdout+checks.stderr)
    sources=['include/ducklips/graph_glottal_storage.hpp','src/graph_glottal_storage.cpp','include/ducklips/glottal_storage.hpp','src/glottal_storage.cpp','include/ducklips/gas_storage.hpp','src/gas_storage.cpp','include/ducklips/moving_graph.hpp','src/moving_graph.cpp','include/ducklips/sources.hpp','include/ducklips/physics.hpp','tools/graph_glottal_probe.cpp','tools/graph_glottal_reference.py']
    args.report.write_text(json.dumps(dict(scope='pressure-driven chamber/tissue/graph coupling with fixed passage coefficients; no moving jet geometry or speech',runs=runs,oracles=oracles,
        native_checks=str(log.relative_to(ROOT)),native_checks_sha256=sha(log),hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),args.binary]},versions=dict(numpy=np.__version__,scipy=scipy.__version__)),indent=2)+'\n')

if __name__=='__main__':main()
