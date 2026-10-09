#!/usr/bin/env python3
"""Independent continuous and discrete checks of coupled one-strip reference folds."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from planar_inertive_motion_reference import ROOT,sha,dense_step

RHO=1.14
REST=3e-6
LENGTH=.014
MASS=np.array([1.25e-4,2.5e-5])
STIFF=np.array([80.,8.])
COUPLING=25.
DAMP=np.array([.1,.6])
SOFT=3e-8
ROUND=2.5e-7


def aperture(fold,rest=REST):
    gaps=rest+2*LENGTH*np.array([fold[0],fold[2]])
    narrow=min(gaps)-ROUND*np.logaddexp(0.,-abs(gaps[0]-gaps[1])/ROUND)
    return float(SOFT*np.logaddexp(0.,narrow/SOFT))


def force(fold,drop,rest=REST,tension=1.):
    x=fold[[0,2]];v=fold[[1,3]];gaps=rest+2*LENGTH*x
    lower=0.
    if min(gaps)>0:lower=drop*(1-(min(gaps)/gaps[0])**2)
    elif gaps[0]>0:lower=drop
    k=STIFF*tension;m=MASS/tension
    forces=np.array([LENGTH*.0025*lower,0.])-k*x-COUPLING*tension*np.array([x[0]-x[1],x[1]-x[0]])
    forces-=3*k*np.minimum(gaps,0)/(2*LENGTH)
    forces-=(DAMP+(gaps<0))*2*np.sqrt(m*k)*v
    return forces/m


def half_step(fold,drop,dt):
    result=fold.copy();result[[1,3]]+=dt*force(fold,drop);result[[0,2]]+=dt*result[[1,3]]
    return result


def coupled_step(fold,free_drop,network_impedance,old_momentum,inertance,dt):
    effective=free_drop+2*old_momentum/dt
    impedance=network_impedance+2*inertance/dt
    def response(drop):
        mid=half_step(fold,drop,.5*dt);area=aperture(mid)
        if area==0:return drop-free_drop,mid,area,0.
        # Solve in velocity coordinates so tiny aperture does not square
        # into numerical zero. This is independent of the native U root.
        linear=impedance*area
        velocity=np.sign(effective)*2*abs(effective)/(linear+np.hypot(linear,np.sqrt(2*RHO*abs(effective)))) if effective else 0.
        actual=.5*RHO*velocity*abs(velocity)
        return drop-actual,mid,area,area*velocity
    lo=min(0.,free_drop,effective);hi=max(0.,free_drop,effective)
    flo=response(lo)[0];fhi=response(hi)[0];tolerance=1e-11*(1+max(abs(lo),abs(hi)))
    if abs(flo)<=tolerance:drop=lo
    elif abs(fhi)<=tolerance:drop=hi
    else:drop=brentq(lambda p:response(p)[0],lo,hi,xtol=1e-10,rtol=1e-14)
    residual,mid,area,u=response(drop);assert abs(residual)<1e-6
    return half_step(mid,drop,.5*dt),area,u,2*inertance*u-old_momentum if area>0 else 0.,drop


def native_fold(record):
    f=np.array(record['fold_state']);assert np.array_equal(f[[0,2,4,6]],f[[1,3,5,7]])
    assert f[8]==REST
    return f[[0,2,4,6]]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_graph_phonation_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/graph-phonation')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/graph_phonation_reference.json')
    args=parser.parse_args();args.binary=args.binary.resolve();args.output=args.output.resolve();args.report=args.report.resolve();args.output.mkdir(parents=True,exist_ok=True)
    c=np.array([2e-10,3e-10,1e-10]);g=np.array([.002,.003,.004,.002]);b=np.array([[1,0,1,1],[-1,1,0,0],[0,-1,-1,0]])
    initial_fold=np.array([1e-6,0.,0.,0.]);inertance=25.;reservoir=700.
    def derivative(t,y):
        p=y[:3];flows=y[3:7];u=y[7];fold=y[8:];area=aperture(fold)
        drop=.5*RHO*(u/area)*abs(u/area)
        charge_rate=-b@flows;charge_rate[-1]+=u
        acceleration=force(fold,drop)
        return np.r_[charge_rate/c,g*(b.T@p),(reservoir-p[-1]-drop)/inertance,
                     fold[1],acceleration[0],fold[3],acceleration[1]]
    y0=np.r_[np.zeros(8),initial_fold];duration=.003
    # This smooth opening prefix avoids claiming a high-order explicit ODE
    # integration of the singular near-contact fluid relaxation later on.
    options=dict(method='DOP853',rtol=2e-12,atol=np.r_[[1e-10]*3,[1e-16]*5,[1e-15]*4])
    oracle=solve_ivp(derivative,(0,duration),y0,max_step=1e-6,**options);assert oracle.success
    refined=solve_ivp(derivative,(0,duration),y0,max_step=5e-7,**options);assert refined.success
    discrepancy=abs(oracle.y[:,-1]-refined.y[:,-1]);assert max(discrepancy[:3])<1e-6 and max(discrepancy[8:])<1e-10
    runs=[];last=None
    def run(name,steps,duration):
        command=[str(args.binary),'--trajectory',str(steps),str(duration)]
        r=subprocess.run(command,capture_output=True,text=True);assert r.returncode==0,r.stderr+r.stdout[-1500:]
        path=args.output/f'{name}-{steps}.json';path.write_text(r.stdout)
        return json.loads(r.stdout),dict(steps=steps,duration_s=duration,command=command,output=str(path.relative_to(ROOT)),output_sha256=sha(path))
    for steps in [144,288,576,1152,2304,4608]:
        native,receipt=run('continuous',steps,duration);v=native['final']
        y=np.r_[np.array(v['charge'])/c,g*np.array(v['momentum']),v['source_momentum']/inertance,native_fold(v)]
        error=abs(y-refined.y[:,-1]);norm=max(error/np.r_[[700.]*3,[1e-4]*5,[1e-4,1.,1e-4,1.]])
        receipt.update(maximum_pressure_error_Pa=float(max(error[:3])),maximum_flow_error_m3_s=float(max(error[3:8])),
                       maximum_fold_position_error_m=float(max(error[[8,10]])),maximum_fold_velocity_error_m_s=float(max(error[[9,11]])),
                       normalized_error=float(norm),refinement_ratio=None if last is None else float(last/norm))
        runs.append(receipt);last=norm
    print(json.dumps(runs),flush=True)
    assert all(1.8<r['refinement_ratio']<2.2 for r in runs[-3:])
    # Full 25 ms, including very small contact apertures: independent stage
    # solve and scalar root, carrying its own fold and air states throughout.
    steps=2400;duration=.025;dt=duration/steps
    z=np.zeros(7);lm=0.;fold=initial_fold.copy();unit=np.array([0.,0.,1.]);zero=np.zeros(3)
    response,response_flow,response_p=dense_step(np.zeros(7),c,c,g,g,unit,dt,b)
    for step in range(steps):
        free,_,free_p=dense_step(z,c,c,g,g,zero,dt,b)
        fold,area,u,lm,drop=coupled_step(fold,reservoir-free_p[-1],response_p[-1],lm,inertance,dt)
        z=free+u*response
    native,discrete=run('discrete',steps,duration);v=native['final']
    discrete.update(maximum_pressure_error_Pa=float(max(abs(np.array(v['charge'])-z[:3])/c)),
                    maximum_momentum_error_Pa_s=float(max(abs(np.array(v['momentum'])-z[3:]))),
                    source_momentum_error_Pa_s=abs(v['source_momentum']-lm),
                    maximum_fold_state_error=float(max(abs(native_fold(v)-fold))),
                    transglottal_error_Pa=abs(v['transglottal_pressure']-drop),aperture_range_m2=native['aperture_range'])
    print(json.dumps(discrete),flush=True)
    assert discrete['maximum_pressure_error_Pa']<1e-3 and discrete['maximum_fold_state_error']<1e-7
    checks=subprocess.run([str(args.binary)],capture_output=True,text=True);assert checks.returncode==0,checks.stderr
    log=args.output/'native-checks.log';log.write_text(checks.stdout+checks.stderr)
    sources=['src/graph_phonation.cpp','include/ducklips/graph_phonation.hpp','src/graph_orifice.cpp','include/ducklips/graph_orifice.hpp',
             'src/sources.cpp','include/ducklips/sources.hpp','src/moving_graph.cpp','tools/graph_phonation_probe.cpp','tools/graph_phonation_reference.py','tools/planar_inertive_motion_reference.py']
    args.report.write_text(json.dumps(dict(use='dynamic fold and total-flow pressure coupling; no speech-quality claim',runs=runs,discrete=discrete,
        continuous_prefix_s=.003,oracle_final=refined.y[:,-1].tolist(),oracle_refinement=discrepancy.tolist(),
        native_checks=str(log.relative_to(ROOT)),native_checks_sha256=sha(log),
        hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),args.binary]},versions=dict(numpy=np.__version__,scipy=scipy.__version__)),indent=2)+'\n')


if __name__=='__main__':main()
