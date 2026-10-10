#!/usr/bin/env python3
"""Independent blocked-profile diffusion and absent-state fixtures."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from slit_viscous_reference import spectrum,ROOT

def main():
    w,p=spectrum();rng=np.random.default_rng(901017);out=ROOT/'out/blocked-viscous-reference';out.mkdir(parents=True,exist_ok=True)
    fixtures=[];references=[]
    # Independent stiff ODE: zero mean flow, remaining profile energy decays.
    # Integrate normalized time so gap/inertia scales cannot weaken tolerance.
    for i in range(24):
        l=10**rng.uniform(0,5);gap=10**rng.uniform(-5,-2);time=10**rng.uniform(-8,-2)
        dt=time*gap*gap/(1.86e-5/1.14);primary=float(rng.normal());metric=np.sqrt(p/(l*w))
        initial=rng.normal(size=24);energy=.5*float(initial@initial)
        sol=solve_ivp(lambda t,y:-p*y,(0,time),initial,method='Radau',jac=-np.diag(p),rtol=2e-10,atol=2e-12)
        assert sol.success;final=sol.y[:,-1]
        fixtures.append(' '.join(map(str,np.r_[l,0.,gap,dt,primary,initial/metric])))
        references.append((primary,final/metric,energy-.5*float(final@final),metric,energy))
    fixtures.append(' '.join(map(str,[0,0,0,1,.17,*([0]*24)])))
    inp=out/'states.input';inp.write_text('\n'.join(fixtures)+'\n')
    result=subprocess.run([str(ROOT/'build/ducklips_slit_viscous_probe'),'--fixture'],input=inp.read_text(),text=True,capture_output=True)
    assert result.returncode==0,result.stderr;dest=out/'states.txt';dest.write_text(result.stdout);data=np.loadtxt(dest)
    state_error=heat_error=0.
    for row,(primary,expected,heat,metric,energy) in zip(data[:-1],references,strict=True):
        assert row[0]==primary
        state_error=max(state_error,float(max(abs(row[1:25]-expected)*metric)))
        heat_error=max(heat_error,abs(row[25]-heat)/energy)
    assert state_error<1e-9 and heat_error<1e-9,(state_error,heat_error)
    assert data[-1,0]==.17 and np.all(data[-1,1:]==0)
    paths=[Path(__file__),ROOT/'tools/slit_viscous_reference.py',ROOT/'src/slit_viscous.cpp',ROOT/'include/ducklips/slit_viscous.hpp',ROOT/'include/ducklips/viscous_spectrum.hpp',ROOT/'tools/slit_viscous_probe.cpp',ROOT/'build/ducklips_slit_viscous_probe',inp,dest]
    report=dict(scope='24 independent stiff ODE profiles at zero through-flow and one exactly absent state. Native transition controls separately test retirement and reopening. No moving-wall PDE or speech.',ode_cases=24,absent_cases=1,maximum_scaled_state_error=state_error,maximum_relative_dissipation_error=heat_error,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/blocked_viscous_reference.json').write_text(json.dumps(report,indent=2)+'\n');print(state_error,heat_error)
if __name__=='__main__':main()
