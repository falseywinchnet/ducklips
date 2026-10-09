#!/usr/bin/env python3
"""Continuous ODE reference for the implicit graph source, including inlet inertia."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from planar_inertive_motion_reference import ROOT,sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_graph_orifice_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/graph-orifice')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/graph_orifice_reference.json')
    args=parser.parse_args();args.binary=args.binary.resolve();args.output=args.output.resolve();args.report=args.report.resolve();args.output.mkdir(parents=True,exist_ok=True)
    c=np.array([2e-10,3e-10,1e-10]);g=np.array([.002,.003,.004]);b=np.array([[1,0,1],[-1,1,0],[0,-1,-1]])
    area=3e-6;rho=1.14;inertance=25.;reservoir=700.
    # Integrate in pressure and flow coordinates, independently of the
    # native charge/momentum stage elimination and Bernoulli root.
    def derivative(t,y):
        p=y[:3];flow=y[3:6];u=y[-1]
        drive=-b@flow;drive[-1]+=u
        return np.r_[drive/c,g*(b.T@p),(reservoir-p[-1]-rho*u*abs(u)/(2*area*area))/inertance]
    reference=solve_ivp(derivative,(0,.004),np.zeros(7),method='DOP853',rtol=2.5e-13,atol=np.r_[[1e-11]*3,[1e-17]*4],max_step=1e-6)
    assert reference.success
    tight=solve_ivp(derivative,(0,.004),np.zeros(7),method='DOP853',rtol=5e-14,atol=np.r_[[2e-12]*3,[2e-18]*4],max_step=5e-7)
    assert tight.success
    discrepancy=abs(reference.y[:,-1]-tight.y[:,-1]);assert max(discrepancy[:3])<1e-7 and max(discrepancy[3:])<1e-12
    runs=[];previous=None
    for count in [32,64,128,256,512,1024,2048,4096,8192,16384]:
        command=[str(args.binary),'--trajectory',str(count)];r=subprocess.run(command,capture_output=True,text=True)
        assert r.returncode==0,r.stderr
        path=args.output/f'trajectory-{count}.json';path.write_text(r.stdout);native=json.loads(r.stdout)
        y=np.r_[np.array(native['charge'])/c,g*np.array(native['momentum']),native['source_momentum']/inertance]
        error=abs(y-tight.y[:,-1]);normalized=float(max(error/np.r_[[700.]*3,[1e-4]*4]))
        ratio=None if previous is None else previous/normalized
        runs.append(dict(steps=count,command=command,maximum_pressure_error_Pa=float(max(error[:3])),maximum_flow_error_m3_s=float(max(error[3:])),
                         normalized_maximum_error=normalized,refinement_ratio=ratio,output=str(path.relative_to(ROOT)),output_sha256=sha(path)))
        previous=normalized
    print(json.dumps(runs),flush=True)
    assert runs[-1]['maximum_pressure_error_Pa']<.001 and runs[-1]['maximum_flow_error_m3_s']<1e-9
    assert all(3.5<r['refinement_ratio']<4.5 for r in runs[-3:]),runs
    checks=subprocess.run([str(args.binary)],capture_output=True,text=True);assert checks.returncode==0,checks.stderr
    log=args.output/'native-checks.log';log.write_text(checks.stdout+checks.stderr)
    sources=['include/ducklips/graph_orifice.hpp','src/graph_orifice.cpp','include/ducklips/moving_graph.hpp','src/moving_graph.cpp',
             'src/sources.cpp','tools/graph_orifice_probe.cpp','tools/graph_orifice_reference.py']
    report=dict(use='continuous nonlinear source/inlet ODE and reservoir energy balance; no speech',command=[sys.executable,*sys.argv],runs=runs,
                oracle_final_state=tight.y[:,-1].tolist(),oracle_pressure_refinement_Pa=float(max(discrepancy[:3])),oracle_flow_refinement_m3_s=float(max(discrepancy[3:])),
                native_checks=str(log.relative_to(ROOT)),native_checks_sha256=sha(log),
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),args.binary]},versions=dict(numpy=np.__version__,scipy=scipy.__version__))
    args.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
