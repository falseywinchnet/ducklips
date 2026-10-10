#!/usr/bin/env python3
"""Independent matrix exponential for the aperture's two inertial states."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'out/radiation-boundary'
RHO=1.14;C=350.;COEFF=1.4*RHO/(4*np.pi*C)
def inverse(a):return np.sqrt(np.pi*a)/(RHO*.7)
def main():
 OUT.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(251019);cases=[];lines=[]
 for i in range(100):
  a,b=10**rng.uniform(-10,-3,2);h0,h1=inverse(a),inverse(b);g=h1*10**rng.uniform(-3,3);r=1/(COEFF*h1*h1)
  dt=10**rng.uniform(-8,3)/(r*(g+h1));p,q=rng.normal(size=2)*.03
  lines.append(' '.join(map(str,[a,b,g,dt,p,q])));cases.append((h0,h1,g,r,dt,p,q))
 fixture='\n'.join(lines)+'\n';(OUT/'matrix.input').write_text(fixture)
 binary=ROOT/'build/ducklips_radiation_boundary_probe';result=subprocess.run([str(binary),'--fixture'],input=fixture,text=True,capture_output=True,check=True)
 (OUT/'matrix.txt').write_text(result.stdout);native=np.loadtxt(OUT/'matrix.txt');errors=[]
 for row,(h0,h1,g,r,dt,p,q) in zip(native,cases):
  scale=np.sqrt([g,h1]);initial=scale*[p,q];direction=np.array([np.sqrt(g),-np.sqrt(h1)])
  expected=expm(-dt*r*np.outer(direction,direction))@initial
  e0=.5*initial@initial;heat=e0-.5*expected@expected;work=.5*(h1-h0)*q*q
  state_error=float(np.max(abs(row[:2]*scale-expected))/(np.linalg.norm(initial)+1e-300))
  heat_error=abs(row[2]-heat)/(e0+1e-300);work_error=abs(row[3]-work)/(abs(work)+e0+1e-300)
  assert max(state_error,heat_error,work_error)<2e-11,(state_error,heat_error,work_error)
  errors.append([state_error,heat_error,work_error])
 paths=[Path(__file__),ROOT/'src/radiation_boundary.cpp',ROOT/'include/ducklips/radiation_boundary.hpp',binary,OUT/'matrix.input',OUT/'matrix.txt']
 report=dict(scope='100 seeded positive-area fixed-coefficient radiation exchanges against SciPy expm in energy-scaled coordinates; includes separately prescribed aperture-work checks',cases=100,seed=251019,maximum_errors=np.max(errors,axis=0).tolist(),hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
 (ROOT/'reports/radiation_boundary_reference.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['maximum_errors']))
if __name__=='__main__':main()
