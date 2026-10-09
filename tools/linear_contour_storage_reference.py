#!/usr/bin/env python3
"""Independent 80-digit storage check for the explicit linear-span profile."""
import json,subprocess,hashlib
from pathlib import Path
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-contact'
def oracle(row):
 a,b,u0,u1=map(mp.mpf,row);H=mp.mpf('.003');B=mp.mpf('139650');g=mp.mpf('1.4');lo=max(mp.mpf(0),-a/(b-a)) if a<=0 else mp.mpf(0);hi=min(mp.mpf(1),a/(a-b)) if b<=0 else mp.mpf(1)
 def f(t,k):
  s=1-t;A=a*s+b*t;u=u0*s+u1*t;lr=mp.log1p((g-1)*u/B)/(g-1);r=mp.exp(lr);p=B/g*mp.expm1(g*lr);e=B/(g*(g-1))*(mp.expm1(g*lr)-g*mp.expm1(lr));C=A*r**(2-g)/B
  return H*[A,A*r*s,A*r*t,A*e,A*p,-p*s,-p*t,C*s*s,C*s*t,C*t*t,r*s*s,r*s*t,r*t*t][k]
 return [mp.quad(lambda t:f(t,k),[lo,hi]) for k in range(13)]
def main():
 rows=[[3e-6,3e-6,0,0],[1e-6,4e-6,700,-200],[5e-6,1e-6,-100,900],[0,3e-6,700,-100],[3e-6,0,-100,700],[-1e-6,3e-6,700,-200],[3e-6,-1e-6,-200,700],[-5e-6,1e-6,-100,900],[1e-6,-3e-6,900,-100],[2e-6,4e-6,-20000,10000],[3e-6,3e-6,1e-9,-1e-9]]
 inp=''.join(' '.join(map(str,r))+'\n' for r in rows);(OUT/'linear-storage-input.txt').write_text(inp);binary=ROOT/'build/ducklips_contour_profile_limit_probe';p=subprocess.run([str(binary),'linear'],input=inp,capture_output=True,text=True);assert p.returncode==0,p.stderr;(OUT/'linear-storage-output.jsonl').write_text(p.stdout);mp.mp.dps=80;records=[]
 for row,line in zip(rows,p.stdout.splitlines()):
  s=json.loads(line);assert s['accepted'];native=[s['volume'],*s['mass'],s['energy'],s['coenergy'],*s['area_effort'],s['compliance'][0],s['compliance'][1],s['compliance'][3],s['mass_area'][0],s['mass_area'][1],s['mass_area'][3]];truth=oracle(row);v=truth[0];scales=[v,v,v,139650*v,139650*v,139650*.003,139650*.003,v/139650,v/139650,v/139650,.003,.003,.003];error=float(max(abs(mp.mpf(x)-y)/scale for x,y,scale in zip(native,truth,scales)));assert error<1e-10;records.append(dict(input=row,maximum_scaled_error=error))
 paths=[Path(__file__),binary,ROOT/'src/contour_pressure_storage.cpp',ROOT/'include/ducklips/contour_pressure_storage.hpp',OUT/'linear-storage-input.txt',OUT/'linear-storage-output.jsonl'];report=dict(digits=80,cases=records,maximum_scaled_error=max(r['maximum_scaled_error'] for r in records),hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});(ROOT/'reports/linear_contour_storage_reference.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS 11 independent linear-span storage cases')
if __name__=='__main__':main()
