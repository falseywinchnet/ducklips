#!/usr/bin/env python3
"""Flat-strip squeeze singularity versus a point-closing affine span."""
import json,subprocess,hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from contour_pressure_loss_reference import conductance
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-contact'
RHO=1.14;GAMMA=1.4;BULK=139650.;WIDTH=.014;ETA=1.86e-5;H=.003;C0=.00125;B=3e-6;U=700.;V=-.0003
r=(1+(GAMMA-1)*U/BULK)**(1/(GAMMA-1));thermal=350*np.sqrt(2/GAMMA)*r**.2

def main():
 areas=[1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,0.,-1e-12,-1e-10,-1e-8];inp=''.join(f'{a} {B} {U} {U}\n' for a in areas);(OUT/'profile-input.txt').write_text(inp);binary=ROOT/'build/ducklips_contour_profile_limit_probe';records={}
 for name in ['flat','linear']:
  p=subprocess.run([str(binary),name],input=inp,text=True,capture_output=True);assert p.returncode==0,p.stderr;(OUT/f'{name}-profile.jsonl').write_text(p.stdout);runs=[]
  for a,line in zip(areas,p.stdout.splitlines()):
   s=json.loads(line);assert s['accepted'],(name,a);D=np.array(s['mass_area']).reshape(2,2);adot=np.array([V,0.]);mdot=D@adot;w=np.r_[0.,-mdot[0],-sum(mdot),adot];L=np.array(s['inertia']).reshape(5,5);R=np.array(s['resistance']).reshape(5,5);kinetic=float(.5*w@L@w);loss=float(w@R@w);row=dict(area=a,kinetic_energy_J=kinetic,loss_W=loss)
   if name=='flat' and a>0:
    delta=a/WIDTH*RHO*r*thermal/(2*ETA);row.update(exact_flat_cap_kinetic_J=RHO*r*V*V*C0**3/(6*a),exact_flat_cap_loss_W=RHO*r*thermal*WIDTH*V*V*C0**3/(3*a*a*conductance(delta)));assert kinetic>=row['exact_flat_cap_kinetic_J']*(1-1e-10) and loss>=row['exact_flat_cap_loss_W']*(1-1e-10)
   runs.append(row)
  records[name]=runs
 exact_kinetic=11*RHO*r*V*V*H**3/(96*B)
 def integrand(t):
  delta=B*t/WIDTH*RHO*r*thermal/(2*ETA);return RHO*r*thermal*WIDTH*V*V*H**3/B**2*(1-t/2)**2/conductance(delta)
 exact_loss=quad(integrand,0.,1.,epsabs=1e-14,epsrel=1e-10)[0];closed=records['linear'][7];assert abs(closed['kinetic_energy_J']/exact_kinetic-1)<1e-11 and abs(closed['loss_W']/exact_loss-1)<1e-7
 approach=records['linear'][6];assert abs(approach['kinetic_energy_J']/exact_kinetic-1)<1e-5 and abs(approach['loss_W']/exact_loss-1)<.002
 # Reversal is an exact coordinate change of the full-span equations. Check
 # all metric/loss entries down to openings far below onset, without flooring.
 reversal=[];perm=[2,1,0,4,3];sign=np.array([-1.,-1.,-1.,1.,1.])
 for a in [1e-8,1e-12,1e-16,1e-20,0.,-1e-12,-1e-8]:
  pair=f'{a} {B} 700 -100\n{B} {a} -100 700\n';p=subprocess.run([str(binary),'linear'],input=pair,text=True,capture_output=True,check=True);left,right=map(json.loads,p.stdout.splitlines());assert left['accepted'] and right['accepted'],a
  errors={}
  for key in ['inertia','resistance']:
   L=np.array(left[key]).reshape(5,5);R=np.array(right[key]).reshape(5,5);expected=L[np.ix_(perm,perm)]*sign[:,None]*sign[None,:];scale=np.sqrt(np.maximum(np.abs(np.diag(expected)),1e-280));error=float(np.max(np.abs(R-expected)/(scale[:,None]*scale[None,:]+1e-280)));assert error<1e-10,(a,key,error);errors[key]=error
  reversal.append(dict(area=a,maximum_scaled_error=errors))
 report=dict(reversal_checks=reversal,scope='prescribed constant-enthalpy squeeze, Qin=0; conditional constitutive limit, not proof of every unconstrained trajectory',profile_meaning=dict(flat='area controls at tissue centres with constant cap strips',linear='area controls at axial boundaries with an affine full span; different explicit shape assumption'),area_rate_m2_s=V,effort_Pa=U,records=records,independent_zero_gap_linear_limit=dict(kinetic_energy_J=exact_kinetic,loss_W=exact_loss),flat_cap_asymptote='Kcap=rho0*r*V^2*c0^3/(6a); Pcap=rho0*r*vth*width*V^2*c0^3/(3 a^2 G(delta)); G grows logarithmically as a tends to zero',hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'tools/contour_profile_limit_probe.cpp',binary,ROOT/'src/contour_pressure_storage.cpp',ROOT/'src/contour_pressure_flow.cpp',ROOT/'src/contour_pressure_loss.cpp',ROOT/'src/contour_pressure_intervals.hpp',OUT/'profile-input.txt',OUT/'flat-profile.jsonl',OUT/'linear-profile.jsonl']})
 (ROOT/'reports/contour_profile_limits.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS native open/closed profile limits, exact kinetic limit, independent BGK loss integral')
if __name__=='__main__':main()
