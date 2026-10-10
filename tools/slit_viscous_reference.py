#!/usr/bin/env python3
"""Independent slit Stokes impedance and full-matrix decay controls."""
import hashlib,json,re,subprocess
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/slit-viscous'
def spectrum():
 mp.mp.dps=40;starts=list(range(9))
 while len(starts)<24:starts.append(int(np.ceil(starts[-1]*1.5)))
 roots=[]
 for n in range(1,starts[-1]+1):
  b=(n+mp.mpf('.5'))*mp.pi;x=mp.findroot(lambda x:mp.sin(x)-x*mp.cos(x),(b-1/b,b-2/b))
  assert n*mp.pi<x<b;roots.append(x)
 weights=[];poles=[];used2=mp.mpf(0);used4=mp.mpf(0)
 for a,b in zip(starts[:-1],starts[1:]):
  s2=mp.fsum(2/x**2 for x in roots[a:b]);s4=mp.fsum(2/x**4 for x in roots[a:b]);used2+=s2;used4+=s4
  poles.append(float(s2/s4));weights.append(float(s2*s2/s4))
 s2=mp.mpf(1)/5-used2;s4=mp.mpf(1)/175-used4
 poles.append(float(s2/s4));weights.append(float(s2*s2/s4))
 return np.array(weights),np.array(poles)
def main():
 OUT.mkdir(parents=True,exist_ok=True);w,p=spectrum();header=(ROOT/'include/ducklips/viscous_spectrum.hpp').read_text()
 for key,values in [('weight',w),('pole',p)]:
  declared=np.array([float(v) for v in re.search(key+r'\{([^}]+)',header).group(1).split(',')]);assert np.max(abs(declared/values-1))<1e-14
 omega=np.geomspace(1e-6,1e8,10001);s=1j*omega;exact=s/(1-np.tanh(np.sqrt(s))/np.sqrt(s));approx=s+3+np.sum(w*s[:,None]/(s[:,None]+p),axis=1)
 impedance=float(max(abs(approx/exact-1)));drag=float(max(abs(approx.real/exact.real-1)));correction=float(max(abs((approx-s)/(exact-s)-1)))
 assert impedance<1.51e-4 and drag<.02 and correction<.02
 assert abs(np.sum(w/p)-.2)<1e-15 and abs(np.sum(w/p**2)-1/175)<1e-16
 np.savez_compressed(OUT/'frequency.npz',dimensionless_omega=omega,exact=exact,grouped=approx)
 rng=np.random.default_rng(291019);fixtures=[];references=[]
 for i in range(100):
  l=10**rng.uniform(0,5);fraction=10**rng.uniform(-4,0);g=fraction/l;gap=10**rng.uniform(-6,-2)
  time=10**rng.uniform(-13,-5);dt=time*gap*gap/(1.86e-5/1.14)
  # Independently form the energy-scaled full generator, including DC drag.
  a=np.zeros((25,25));a[0,0]=-fraction*(3+w.sum());a[0,1:]=a[1:,0]=np.sqrt(fraction*w*p);a[1:,1:]=-np.diag(p)
  metric=np.sqrt(np.r_[g,p/(l*w)]);x=rng.normal(size=25)*.001;momentum=x/metric;expected=expm(a*time)@x
  fixtures.append(' '.join(map(str,np.r_[l,g,gap,dt,momentum])))
  references.append((expected/metric,.5*(x@x-expected@expected),metric,.5*(x@x)))
 path=OUT/'matrix.input';path.write_text('\n'.join(fixtures)+'\n')
 result=subprocess.run([str(ROOT/'build/ducklips_slit_viscous_probe'),'--fixture'],input=path.read_text(),text=True,capture_output=True)
 assert result.returncode==0,result.stderr;(OUT/'matrix.txt').write_text(result.stdout);data=np.loadtxt(OUT/'matrix.txt')
 state_error=heat_error=0.
 for row,(expected,heat,metric,energy) in zip(data,references):
  state_error=max(state_error,float(np.max(abs(row[:25]-expected)*metric)/.001));heat_error=max(heat_error,abs(row[25]-heat)/energy)
 assert state_error<1e-10 and heat_error<1e-10,(state_error,heat_error)
 paths=[Path(__file__),ROOT/'src/slit_viscous.cpp',ROOT/'include/ducklips/slit_viscous.hpp',ROOT/'include/ducklips/viscous_spectrum.hpp',ROOT/'tools/slit_viscous_probe.cpp',ROOT/'build/ducklips_slit_viscous_probe',*sorted(OUT.iterdir())]
 report=dict(scope='No-slip stationary-wall slit impedance, with full branch inertia shared by a segment. Exact Stokes transfer and independent matrix decay; no anatomical calibration or full speech.',matrix_cases=100,maximum_scaled_state_error=state_error,maximum_relative_energy_error=heat_error,maximum_relative_impedance_error=impedance,maximum_relative_resistance_error=drag,maximum_relative_viscous_correction_error=correction,dimensionless_frequency_range=[1e-6,1e8],profile_inertia_fraction=float(np.sum(w/p)),hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
 (ROOT/'reports/slit_viscous_reference.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='hashes'},indent=2))
if __name__=='__main__':main()
