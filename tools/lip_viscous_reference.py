#!/usr/bin/env python3
"""Spatial lip-impedance quadrature against direct integration of slit Stokes flow."""
import csv,hashlib,json,subprocess
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad_vec,quad
from lip_profile_reference import oracle,vector,ROOT
from slit_viscous_reference import spectrum
OUT=ROOT/'out/lip-viscous-reference';RHO=1.14;NU=1.86e-5/RHO
def area_function(p):
 length,base,activation,place,core,width,target=p;activation=np.clip(activation,0,1);target=min(base,target)
 minimum=target+(base-target)*(1-activation)**2*(1+2*activation);tail=np.sqrt(-np.log(1e-4))
 def area(x):
  u=np.maximum(abs(x-place)-core,0)/width
  return np.where(u>tail,base,minimum+(base-minimum)*(-np.expm1(-u*u)))
 return area,[place-core,place+core,place-core-tail*width,place+core+tail*width]
def lip_segments(profile,order,width):
 metrics=oracle(profile);area,_=area_function(profile);x,w=leggauss(order)
 # Native identities alternate the negative/positive abscissae by radius.
 positive=x[order//2:];weights=w[order//2:];x=np.column_stack([-positive,positive]).ravel();w=np.repeat(weights,2)
 records=[]
 for path in range(1,4):
  a,b=metrics['centroid'][path],metrics['centroid'][path-1];position=(a+b)/2+(b-a)/2*x;areas=area(position)
  inertia=RHO*(b-a)/2*w/areas;inertia/=sum(inertia)*metrics['conductance'][path]
  for j in range(order):records.append(dict(branch=2+path,position=float(position[j]),segment_inertance=float(inertia[j]),branch_conductance=float(metrics['conductance'][path]),half_gap=float(areas[j]/(2*width))))
 return records
def main():
 OUT.mkdir(parents=True,exist_ok=True);rows=list(csv.DictReader((ROOT/'out/lip-profile/nato.csv').open()));profiles=[]
 for start in [1284,1374]:
  for k in range(9):
   t=start+k/4;i=int(t);fraction=t-i;profiles.append((t,(1-fraction)*vector(rows[i])+fraction*vector(rows[i+1])))
 frequencies=np.geomspace(20,200000,257);s=2j*np.pi*frequencies;width=.03;w,poles=spectrum()
 fixtures=[];cases=[]
 for time,profile in profiles:
  metrics=oracle(profile);area,cuts=area_function(profile)
  for path in range(1,4):
   a,b=metrics['centroid'][path],metrics['centroid'][path-1];points=[c for c in cuts if a<c<b]
   def integrand(x):
    value=area(x);u=np.sqrt(s*(value/(2*width))**2/NU)
    return RHO/value*s*np.tanh(u)/(u-np.tanh(u))
   exact=quad_vec(integrand,a,b,epsabs=0,epsrel=2e-11,points=points)[0]
   static=12*1.86e-5*width**2*quad(lambda x:1/area(x)**3,a,b,epsabs=0,epsrel=2e-11,points=points)[0]
   inertia=1/metrics['conductance'][path];effective_gap=np.sqrt(3*inertia*NU/static);u=np.sqrt(s*effective_gap**2/NU)
   uniform=inertia*s*np.tanh(u)/(u-np.tanh(u))
   for order in [4,8,16,32]:
    fixtures.append(str(order)+' '+' '.join(map(str,np.r_[profile,a,b,width])))
    cases.append((time,path,order,profile,exact,inertia,float(max(abs(uniform/exact-1)))))
 source=OUT/'geometry.input';source.write_text('\n'.join(fixtures)+'\n')
 result=subprocess.run([str(ROOT/'build/ducklips_lip_viscous_probe'),'--fixture'],input=source.read_text(),text=True,capture_output=True)
 assert result.returncode==0,result.stderr;(OUT/'geometry.txt').write_text(result.stdout)
 errors={n:dict(spatial_correction=0.,combined_correction=0.,combined_total_impedance=0.,inertance_normalization=0.,geometry=0.) for n in [4,8,16,32]};uniform_error=0
 for line,(time,path,n,profile,exact,inertia,uniform) in zip(result.stdout.splitlines(),cases):
  values=np.fromstring(line,sep=' ');segments=values[3:].reshape(n,3);reference=lip_segments(profile,n,width)[(path-1)*n:path*n]
  expected=np.array([[r['position'],r['segment_inertance'],r['half_gap']] for r in reference]);geometry=float(max(abs((segments/expected-1)).ravel()))
  l=segments[:,1];gap=segments[:,2];u=np.sqrt(s[:,None]*gap**2/NU)
  spatial=np.sum(l*s[:,None]*np.tanh(u)/(u-np.tanh(u)),axis=1)
  rates=NU/gap**2;local=3+np.sum(w*s[:,None,None]/(s[:,None,None]+rates[None,:,None]*poles),axis=2)
  grouped=np.sum(l*rates*local,axis=1);e=errors[n]
  e['spatial_correction']=max(e['spatial_correction'],float(max(abs(spatial/exact-1))))
  e['combined_correction']=max(e['combined_correction'],float(max(abs(grouped/exact-1))))
  e['combined_total_impedance']=max(e['combined_total_impedance'],float(max(abs((grouped-exact)/(s*inertia+exact)))))
  e['inertance_normalization']=max(e['inertance_normalization'],abs(values[2]));e['geometry']=max(e['geometry'],geometry);uniform_error=max(uniform_error,uniform)
  assert abs(sum(l)*values[0]-1)<1e-13 and geometry<1e-10
 assert errors[8]['spatial_correction']<3e-4 and errors[16]['spatial_correction']<1e-7 and errors[32]['spatial_correction']<1e-10,errors
 assert all(e['combined_correction']<.08 for e in errors.values()),errors
 paths=[Path(__file__),ROOT/'tools/lip_profile_reference.py',ROOT/'tools/slit_viscous_reference.py',ROOT/'src/lip_viscous.cpp',ROOT/'include/ducklips/lip_viscous.hpp',ROOT/'src/lip_profile.cpp',ROOT/'tools/lip_viscous_probe.cpp',ROOT/'build/ducklips_lip_viscous_probe',ROOT/'out/lip-profile/nato.csv',source,OUT/'geometry.txt']
 report=dict(scope='Positive quadrature of local stationary-wall slit impedances along the three open internal lip paths, retaining exact integrated plug inertia. Direct adaptive spatial integration, no full moving-wall PDE or speech claim.',profiles=len(profiles),native_fixtures=len(cases),frequency_range_Hz=[20,200000],frequencies=len(frequencies),slit_width_m=width,maximum_relative_errors_by_order=errors,maximum_correction_error_single_DC_matched_slit=uniform_error,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
 (ROOT/'reports/lip_viscous_reference.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='hashes'},indent=2))
if __name__=='__main__':main()
