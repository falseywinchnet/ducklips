#!/usr/bin/env python3
"""Measure the constitutive closing-cap limit, not a dynamic event solution."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from contour_sealed_reference import metrics
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=ROOT/'out/contour-sealed';out.mkdir(parents=True,exist_ok=True);binary=ROOT/'build/ducklips_contour_sealed_probe';fixtures=[];labels=[]
    for cell in range(2):
        closed=np.array([3e-6,3e-6]);closed[cell]=0;mass=np.array(metrics(np.r_[1e-9,1e-9,closed])[0])
        for area in np.geomspace(1e-8,1e-13,11):
            shape=closed.copy();shape[cell]=area;fixtures.append([*mass,*shape,0,0,0]);labels.append((cell,float(area)))
        fixtures.append([*mass,*closed,0,0,0]);labels.append((cell,0.))
    inp=out/'sealing-limit-input.txt';inp.write_text(''.join(' '.join(map(str,r))+'\n' for r in fixtures));run=subprocess.run([str(binary),'--kernel'],input=inp.read_text(),capture_output=True,text=True);assert run.returncode==0,run.stderr
    output=out/'sealing-limit-output.jsonl';output.write_text(run.stdout);native=[json.loads(s) for s in run.stdout.splitlines()];assert len(native)==len(fixtures);rows=[];checks=[]
    for (cell,area),state in zip(labels,native):
        assert state['accepted'];j=0 if cell==0 else 2;g=state['inverse'][j];r=state['resistance'][j]
        if area==0:
            assert g==r==0 and state['volume'][cell]>0;rows.append(dict(cell=cell,area_m2=0.,sealed=True,retained_volume_m3=state['volume'][cell]));continue
        tau=1/(g*r);pressure_effort=700.;momentum=tau*pressure_effort;flow=pressure_effort/r;energy=.5*g*momentum**2
        rows.append(dict(cell=cell,area_m2=area,inverse_inertia=g,resistance=r,viscous_relaxation_s=tau,quasisteady_flow_m3_s=flow,quasisteady_momentum_Pa_s=momentum,quasisteady_energy_J=energy))
    for cell in range(2):
        data=[r for r in rows if r['cell']==cell and r['area_m2']>0][-4:];x=np.log([r['area_m2'] for r in data]);exponents={key:float(np.polyfit(x,np.log([r[key] for r in data]),1)[0]) for key in ['inverse_inertia','resistance','viscous_relaxation_s','quasisteady_flow_m3_s','quasisteady_momentum_Pa_s','quasisteady_energy_J']};expected=[1,-3,2,3,2,5]
        checks.append(dict(cell=cell,exponents=exponents,passed=all(abs(v-e)<1e-4 for v,e in zip(exponents.values(),expected))))
    report=dict(scope='constitutive and quasisteady limiting powers at fixed positive chamber masses and 700 Pa mass effort drop; not a simulated contact event',conditions=['positive viscosity','positive retained chamber volume','finite nonzero chamber density','finite passage driving effort'],rows=rows,checks=checks,hashes={str(p.relative_to(ROOT)):sha(p) for p in [binary,Path(__file__),ROOT/'tools/contour_sealed_reference.py',inp,output]})
    (ROOT/'reports/contour_sealing_limit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(checks));assert all(c['passed'] for c in checks)
if __name__=='__main__':main()
