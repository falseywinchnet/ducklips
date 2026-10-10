#!/usr/bin/env python3
"""Expose steady-throughput bias in separated pressure/viscous updates."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from planar_continuous_reference import decay_exponential
from slit_viscous_reference import spectrum,ROOT

def main():
    out=ROOT/'out/viscous-forced-flow';out.mkdir(parents=True,exist_ok=True);w,poles=spectrum()
    inertia=100.;gap=1e-5;pressure=700.;nu=1.86e-5/1.14;g=1/inertia
    matrix=np.zeros((25,25));matrix[0,0]=-(3+sum(w));matrix[0,1:]=matrix[1:,0]=np.sqrt(w*poles);matrix[1:,1:]=-np.diag(poles)
    parameters=[];expected=[]
    for stiffness in [.01,.1,1.,10.,100.,10000.]:
        dt=stiffness*gap**2/(3*nu);half=decay_exponential(matrix*(dt/2*nu/gap**2))
        force=np.zeros(25);force[0]=np.sqrt(g)*pressure*dt
        # Fixed point of E (E y + dt f), independently of time marching.
        state=np.linalg.solve(np.eye(25)-half@half,half@force)
        flow=np.sqrt(g)*((half@state)[0]+.5*force[0])
        expected.append(flow);parameters.append((stiffness,dt))
    inp=out/'states.input';inp.write_text('\n'.join(f'{inertia} {gap} {pressure} {dt} 8000' for _,dt in parameters)+'\n')
    result=subprocess.run([str(ROOT/'build/ducklips_viscous_forced_flow_probe')],input=inp.read_text(),text=True,capture_output=True)
    assert result.returncode==0,result.stderr;dest=out/'states.txt';dest.write_text(result.stdout);data=np.loadtxt(dest);cases=[]
    for (stiffness,dt),row,predicted in zip(parameters,data,expected,strict=True):
        assert abs(row[0]/predicted-1)<1e-10
        assert abs(row[3])<1e-12*max(abs(dt*pressure*row[0]),1e-30)
        cases.append(dict(dc_decay_times_step=stiffness,step_seconds=dt,split_average_flow_m3_s=row[0],split_endpoint_flow_m3_s=row[1],continuous_steady_flow_m3_s=row[2],average_flow_ratio=row[0]/row[2],last_step_work_residual_J=row[3]))
    paths=[Path(__file__),ROOT/'tools/viscous_forced_flow_probe.cpp',ROOT/'tools/planar_continuous_reference.py',ROOT/'tools/slit_viscous_reference.py',ROOT/'src/slit_viscous.cpp',ROOT/'include/ducklips/slit_viscous.hpp',ROOT/'build/ducklips_viscous_forced_flow_probe',inp,dest]
    report=dict(scope='Uniform fixed slit under constant pressure: same separated drag/pressure/drag composition used by the graph adapter. A diagnostic of forced-flow consistency, not a full airway error estimate.',finding='Exact unforced decay and a closed work ledger do not make the separated scheme asymptotically correct for stiff driven flow. Pressure and drag need a common update.',cases=cases,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/viscous_forced_flow_reference.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(cases,indent=2))
if __name__=='__main__':main()
