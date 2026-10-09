#!/usr/bin/env python3
"""Independent stiff integration while the full-span passage remains open."""
import json,sys,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from continuous_phonation_reference import Reference,YSCALE
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-phonation'
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 duration=float(sys.argv[1]) if len(sys.argv)>1 else .06;tol=float(sys.argv[2]) if len(sys.argv)>2 else 2e-7;name=sys.argv[3] if len(sys.argv)>3 else 'open-ode';order=int(sys.argv[4]) if len(sys.argv)>4 else 32;ref=Reference(order);initial=np.array(list(map(float,(OUT/'mature.input').read_text().split()))) / YSCALE if (OUT/'mature.input').exists() else np.array(json.loads((ROOT/'reports/continuous_phonation_motion.json').read_text())['seeds']['mature_state'])/YSCALE
 def near_contact(t,y):return min(y[2:4])-.01
 near_contact.terminal=True;near_contact.direction=-1;start=time.monotonic();sol=solve_ivp(lambda t,y:ref.rhs(y,1400),(0.,duration),initial,method='Radau',rtol=tol,atol=tol*.01,max_step=.0001,dense_output=True,events=near_contact);elapsed=time.monotonic()-start
 times=np.r_[np.arange(0.,sol.t[-1],.00005),sol.t[-1]];states=(sol.sol(times).T*YSCALE).tolist();report=dict(scope='independent open-domain ODE; stops at 1e-8 m2 area before contact, without flooring area or changing equations',method='Radau',quadrature_order=order,rtol=tol,atol=tol*.01,maximum_step_s=.0001,pressure_Pa=1400,initial_scaled=initial.tolist(),success=sol.success,message=sol.message,final_time=float(sol.t[-1]),near_contact_event=sol.t_events[0].tolist(),seconds=elapsed,nfev=sol.nfev,njev=sol.njev,times=times.tolist(),states=states);(OUT/(name+'.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['times','states','initial_scaled']}),flush=True);assert sol.success
if __name__=='__main__':main()
