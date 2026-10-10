#!/usr/bin/env python3
"""Finite-slit heat diffusion, independent matrix evolution and static storage."""
import hashlib,json,re,subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.special import polygamma
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/slit-thermal'
def spectrum():
 starts=list(range(9))
 while len(starts)<24:starts.append(int(np.ceil(starts[-1]*1.5)))
 weights=[];poles=[]
 for a,b in zip(starts,starts[1:]+[np.inf]):
  s2=polygamma(1,a+.5)-(polygamma(1,b+.5) if np.isfinite(b) else 0.)
  s4=(polygamma(3,a+.5)-(polygamma(3,b+.5) if np.isfinite(b) else 0.))/6
  weights.append(2/np.pi**2*s2);poles.append(np.pi**2*s2/s4)
 return np.array(weights),np.array(poles)
def main():
 OUT.mkdir(parents=True,exist_ok=True);w,p=spectrum();header=(ROOT/'include/ducklips/thermal_spectrum.hpp').read_text()
 for key,values in [('weight',w),('pole',p)]:
  declared=np.array([float(v) for v in re.search(key+r'\{([^}]+)',header).group(1).split(',')])
  assert np.max(abs(declared/values-1))<1e-14
 omega=np.geomspace(1e-6,1e8,10001);s=1j*omega;exact=np.tanh(np.sqrt(s))/np.sqrt(s)
 approx=np.sum(w*p/(s[:,None]+p),axis=1)
 storage_error=float(np.max(abs((1+.4*approx)/(1+.4*exact)-1)))
 assert storage_error<8e-5 and abs(sum(w)-1)<1e-14
 np.savez_compressed(OUT/'frequency.npz',dimensionless_omega=omega,exact=exact,grouped=approx)
 # Independent full 25-coordinate generator in energy-scaled coordinates.
 # No elimination of the conserved mode and no native eigensolver.
 k=.4*w;A=np.zeros((25,25));A[0,0]=-np.sum(k*p)
 A[0,1:]=A[1:,0]=p*np.sqrt(k);A[1:,1:]=-np.diag(p)
 rng=np.random.default_rng(271019);fixtures=[];references=[];initial=[]
 for i in range(100):
  c=10**rng.uniform(-13,-8);gap=10**rng.uniform(-5,-2)
  dimensionless=10**rng.uniform(-13,-5);dt=dimensionless*gap*gap/2.3e-5
  x=rng.normal(size=25);charges=x*np.r_[1,np.sqrt(k)]*c*700
  expected=expm(A*dimensionless)@x
  heat=.5*c*700**2*(x@x-expected@expected)
  fixtures.append(' '.join(map(str,np.r_[c,gap,dt,charges])))
  references.append((expected*np.r_[1,np.sqrt(k)]*c*700,heat));initial.append((c,x))
 path=OUT/'matrix.input';path.write_text('\n'.join(fixtures)+'\n')
 result=subprocess.run([str(ROOT/'build/ducklips_slit_thermal_probe'),'--fixture'],input=path.read_text(),capture_output=True,text=True)
 assert result.returncode==0,result.stderr;(OUT/'matrix.txt').write_text(result.stdout)
 values=np.loadtxt(OUT/'matrix.txt');state_error=heat_error=0.
 for row,(expected,heat),(c,x) in zip(values,references,initial):
  state_error=max(state_error,float(np.max(abs(row[:25]-expected)/(c*700))))
  heat_error=max(heat_error,abs(row[25]-heat)/(.5*c*700**2*(x@x)))
 assert state_error<1e-10 and heat_error<1e-10,(state_error,heat_error)
 files=[Path(__file__),ROOT/'include/ducklips/thermal_spectrum.hpp',ROOT/'include/ducklips/slit_thermal.hpp',ROOT/'src/slit_thermal.cpp',ROOT/'tools/slit_thermal_probe.cpp',ROOT/'build/ducklips_slit_thermal_probe',*sorted(OUT.iterdir())]
 report=dict(scope='Isothermal-wall slit thermal diffusion; grouped positive spectrum compared with exact tanh transfer and independent matrix exponential. Linear local model, no moving-wall advection or full speech.',matrix_cases=100,maximum_scaled_state_error=state_error,maximum_relative_energy_error=heat_error,maximum_relative_complex_compliance_error=storage_error,dimensionless_frequency_range=[1e-6,1e8],static_thermal_weight_sum=float(sum(w)),hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
 (ROOT/'reports/slit_thermal_reference.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='hashes'},indent=2))
if __name__=='__main__':main()
