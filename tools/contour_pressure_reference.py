#!/usr/bin/env python3
"""Independent thermodynamic audit of continuous contour pressure storage."""
import hashlib,json,subprocess
from pathlib import Path
import mpmath as mp
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/contour-pressure'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def oracle(row):
    a,b,u0,u1=map(mp.mpf,row);c0=mp.mpf('.00125');c1=mp.mpf('.00275');end=mp.mpf('.003');bulk=mp.mpf('1.14')*350**2;gamma=mp.mpf('1.4');edges=[0,c0,c1,end];total=[mp.mpf(0)]*13
    def theta(z):return max(mp.mpf(0),min(mp.mpf(1),(z-c0)/(c1-c0)))
    def opening(z):t=theta(z);return (1-t)*a+t*b
    def f(z,k):
        t=theta(z);u=1-t;area=opening(z);mu=u*u0+t*u1;r=(1+(gamma-1)*mu/bulk)**(1/(gamma-1));p=bulk/gamma*(r**gamma-1);energy=bulk/(gamma*(gamma-1))*(r**gamma-gamma*r+gamma-1);slope=r**(2-gamma)/bulk
        return [area,area*r*u,area*r*t,area*energy,area*p,-p*u,-p*t,area*slope*u*u,area*slope*u*t,area*slope*t*t,r*u*u,r*u*t,r*t*t][k]
    for lo,hi in zip(edges,edges[1:]):
        x,y=opening(lo),opening(hi)
        if x<=0 and y<=0:continue
        if x<=0:lo+=(hi-lo)*(-x)/(y-x)
        elif y<=0:hi=lo+(hi-lo)*x/(x-y)
        for k in range(13):total[k]+=mp.quad(lambda z:f(z,k),[lo,hi])
    return total

def unpack(s):return [s['volume'],*s['mass'],s['energy'],s['coenergy'],*s['area_effort'],s['compliance'][0],s['compliance'][1],s['compliance'][3],s['mass_area'][0],s['mass_area'][1],s['mass_area'][3]]
def main():
    OUT.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(92117);rows=[]
    for n in range(72):
        a,b=10**rng.uniform(-7,-5,2)
        if 24<=n<48:a=-b*rng.uniform(.1,12)
        if n>=48:b=-a*rng.uniform(.1,3)
        if n in [24,25,26]:a=-b*[4.999999,5.,5.000001][n-24]
        effort=rng.uniform(-1200,1800,2) if n<60 else rng.uniform(-300000,200000,2)
        if n<3:effort=np.array([0.,0.]) if n==0 else np.array([700.,700.]) if n==1 else np.array([-100.,-100.])
        rows.append([a,b,*effort])
    inp=OUT/'input.txt';inp.write_text(''.join(' '.join(map(str,row))+'\n' for row in rows));binary=ROOT/'build/ducklips_contour_pressure_probe';run=subprocess.run([str(binary),'--batch'],input=inp.read_text(),capture_output=True,text=True);assert run.returncode==0,run.stderr
    output=OUT/'output.jsonl';output.write_text(run.stdout);native=[json.loads(line) for line in run.stdout.splitlines()];assert len(native)==len(rows)
    maxima=np.zeros(13);recover_rows=[];mp.mp.dps=80;records=[]
    for index,(row,s) in enumerate(zip(rows,native)):
        expected=oracle(row);volume=float(expected[0]);scale=np.array([volume,volume,volume,*([139650*volume]*2),*([139650*.003]*2),*([volume/139650]*3),*([.003]*3)])
        errors=np.array([float(abs(mp.mpf(v)-e)) for v,e in zip(unpack(s),expected)])/scale;maxima=np.maximum(maxima,errors)
        recover_rows.append([row[0],row[1],float(expected[1]),float(expected[2])]);records.append(dict(input=row,maximum_scaled_error=float(max(errors))))
        assert s['compliance'][0]*s['compliance'][3]>s['compliance'][1]**2
        if index%12==11:print('checked',index+1,'maximum',float(max(maxima)),flush=True)
    assert max(maxima)<1e-10,maxima
    recovery_input=OUT/'recover-input.txt';recovery_input.write_text(''.join(' '.join(map(str,row))+'\n' for row in recover_rows));run=subprocess.run([str(binary),'--recover'],input=recovery_input.read_text(),capture_output=True,text=True);assert run.returncode==0,run.stderr
    recovery_output=OUT/'recover-output.jsonl';recovery_output.write_text(run.stdout);recovered=[json.loads(line) for line in run.stdout.splitlines()];assert len(recovered)==len(rows)
    max_effort=max(max(abs(np.array(s['effort'])-row[2:])) for s,row in zip(recovered,rows));assert max_effort<1e-4,max_effort
    sweeps=[];previous=None
    for steps in [100,200,400,800]:
        command=[str(binary),'--sweep',str(steps)];run=subprocess.run(command,capture_output=True,text=True);assert run.returncode==0,run.stderr
        path=OUT/f'sweep-{steps}.json';path.write_text(run.stdout);s=json.loads(run.stdout);error=abs(s['work_residual_J']);ratio=None if previous is None else previous/error
        assert s['maximum_mass_error_m3']<1e-20 and s['maximum_relative_density']<1.32
        sweeps.append(dict(steps=steps,work_residual_J=s['work_residual_J'],work_refinement_ratio=ratio,maximum_mass_error_m3=s['maximum_mass_error_m3'],maximum_relative_density=s['maximum_relative_density'],output=str(path.relative_to(ROOT)),output_sha256=sha(path)));previous=error
    assert all(3.5<s['work_refinement_ratio']<4.5 for s in sweeps[-2:]) and previous<1e-14,sweeps
    report=dict(scope='continuous enthalpy field storage and inverse mass moments; no momentum or speech integration',seam_sweeps=sweeps,cases=len(rows),digits=80,maximum_scaled_error=float(max(maxima)),scaled_errors=dict(zip(['volume','mass0','mass1','energy','coenergy','force0','force1','C00','C01','C11','D00','D01','D11'],maxima)),maximum_recovered_effort_error_Pa=float(max_effort),cases_detail=records,hashes={str(p.relative_to(ROOT)):sha(p) for p in [inp,output,recovery_input,recovery_output,binary,ROOT/'include/ducklips/contour_pressure_storage.hpp',ROOT/'src/contour_pressure_storage.cpp',ROOT/'tools/contour_pressure_probe.cpp',Path(__file__)]})
    (ROOT/'reports/contour_pressure_reference.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['hashes','cases_detail']},indent=2))
if __name__=='__main__':main()
