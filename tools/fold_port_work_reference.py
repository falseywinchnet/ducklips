#!/usr/bin/env python3
"""Audit available glottal port work against the inherited aerodynamic fold force.

Reports an unresolved physical work-accounting gap; it is not a speech test.
"""
import csv
import json
import subprocess
from pathlib import Path
import numpy as np
from graph_phonation_reference import coupled_step,half_step,native_fold,REST,LENGTH
from planar_inertive_motion_reference import ROOT,sha,dense_step


def main():
    out=ROOT/'out/fold-port-work';out.mkdir(exist_ok=True)
    c=np.array([2e-10,3e-10,1e-10]);g=np.array([.002,.003,.004,.002]);b=np.array([[1,0,1,1],[-1,1,0,0],[0,-1,-1,0]])
    binary=ROOT/'build/ducklips_graph_phonation_probe';cases=[]
    for count in [2400,4800,9600]:
        dt=.025/count;z=np.zeros(7);lm=0.;fold=np.array([1e-6,0.,0.,0.]);unit=np.array([0.,0.,1.]);zero=np.zeros(3)
        response,_,rp=dense_step(np.zeros(7),c,c,g,g,unit,dt,b)
        aero=0.;available=0.;deficit=0.;peak=0.;negative_steps=0;rows=[]
        for i in range(count):
            free,_,fp=dense_step(z,c,c,g,g,zero,dt,b);old=fold.copy()
            fold,area,u,lm,drop=coupled_step(fold,700-fp[-1],rp[-1],lm,25,dt)
            mid=half_step(old,drop,.5*dt);work=0.
            for start,end in [(old,mid),(mid,fold)]:
                gap=REST+2*LENGTH*start[[0,2]]
                drive=drop*(1-(min(gap)/gap[0])**2) if min(gap)>0 else drop if gap[0]>0 else 0.
                # Work of the aerodynamic force over the displacement that
                # its symplectic half-step accepts, summed over both folds.
                work+=2*LENGTH*.0025*drive*(end[0]-start[0])
            port=dt*drop*u;shortfall=work-port;aero+=work;available+=port
            if shortfall>0:negative_steps+=1;deficit+=shortfall;peak=max(peak,shortfall/dt)
            rows.append([(i+1)*dt,area,u,drop,work/dt,port/dt,shortfall/dt])
            z=free+u*response
        path=out/f'work-{count}.csv'
        with path.open('w') as f:
            writer=csv.writer(f);writer.writerow(['time_s','midpoint_aperture_m2','source_flow_m3_s','glottal_drop_Pa','fold_aerodynamic_power_W','orifice_port_power_W','shortfall_W']);writer.writerows(rows)
        command=[str(binary),'--trajectory',str(count),'.025'];r=subprocess.run(command,capture_output=True,text=True);assert r.returncode==0,r.stderr
        output=out/f'native-{count}.json';output.write_text(r.stdout);native=json.loads(r.stdout)['final']
        pressure_error=float(max(abs(np.array(native['charge'])-z[:3])/c));fold_error=float(max(abs(native_fold(native)-fold)))
        assert pressure_error<1e-3 and fold_error<1e-7
        case=dict(steps=count,dt_s=dt,duration_s=.025,fold_aerodynamic_work_J=aero,orifice_port_work_J=available,
                  integrated_positive_shortfall_J=deficit,maximum_shortfall_W=peak,shortfall_steps=negative_steps,
                  native_pressure_difference_Pa=pressure_error,native_fold_state_difference=fold_error,
                  command=command,trace=str(path.relative_to(ROOT)),trace_sha256=sha(path),native=str(output.relative_to(ROOT)),native_sha256=sha(output))
        cases.append(case);print(json.dumps(case),flush=True)
    sources=['tools/fold_port_work_reference.py','tools/graph_phonation_reference.py','tools/planar_inertive_motion_reference.py',
             'src/graph_phonation.cpp','src/graph_orifice.cpp','src/sources.cpp','include/ducklips/sources.hpp','tools/graph_phonation_probe.cpp']
    report=dict(finding='unresolved: inherited fold aerodynamic work exceeds instantaneous orifice port work during parts of closure',
                interpretation='The current algebraic glottis has no internal air-storage/displacement state to account for this interval work difference. Total net orifice work remains larger over the whole control; this is not a claim of net energy creation over 25 ms.',
                work_accounting_complete=False,cases=cases,hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),binary]})
    (ROOT/'reports/fold_port_work_reference.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
