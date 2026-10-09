#!/usr/bin/env python3
"""Check local seam geometry against independent 80-digit clipped integrals."""
import hashlib,json
from pathlib import Path
import mpmath as mp
import numpy as np
from contour_pocket_reference import metric,density_factor
from contour_sealed_reference import metrics
ROOT=Path(__file__).resolve().parents[1]
mp.mp.dps=80
maxima=dict(volume=0.,inverse_inertia=0.,inverse_area_cubed=0.,energy_gradient=0.,volume_gradient=0.)
rows=[]
for seam in np.geomspace(1e-7,1e-16,40):
 for density in [1.,1e6,1e12]:
    area=-1.8e-6;v=metric([1e-9,1e-9,area,seam])[0];c=np.r_[v*np.array([density,1.003]),area,seam];p=np.array([.001,.002]);native=metric(c)
    c0=list(map(mp.mpf,c));p0=list(map(mp.mpf,p))
    def exact(x):return metrics([x[0],x[1],x[2],(6*x[3]-x[2])/5],mp.log)
    ev,eg,ei=exact(c0);e=sum(p0[j]**2*eg[j+1]/2 for j in range(2))
    for k in range(2):maxima['volume']=max(maxima['volume'],float(abs(native[0][k]/ev[k]-1)));maxima['inverse_inertia']=max(maxima['inverse_inertia'],float(abs(native[1][k]/eg[k+1]-1)))
    for j in range(2):
     for k in range(2):
      if ei[j+1][k]:maxima['inverse_area_cubed']=max(maxima['inverse_area_cubed'],float(abs(native[2][j,k]/ei[j+1][k]-1)))
    for k in range(4):
     h=abs(c[k])*1e-20;z=c.astype(complex);z[k]+=h*1j;vz,gz,_,_=metric(z);derivative=np.imag(np.sum(p*p*gz/2))/h
     expected=mp.diff(lambda value:sum(p0[j]**2*exact([value if n==k else c0[n] for n in range(4)])[1][j+1]/2 for j in range(2)),c0[k]);maxima['energy_gradient']=max(maxima['energy_gradient'],float(abs(derivative-expected)*abs(c0[k])/e))
     if k>=2:
      for j in range(2):
       expected=mp.diff(lambda value:exact([value if n==k else c0[n] for n in range(4)])[0][j],c0[k]);maxima['volume_gradient']=max(maxima['volume_gradient'],float(abs(vz[j].imag/h-expected)*abs(c0[k])/ev[j]))
    rows.append(dict(seam=float(seam),relative_density=density))
assert max(maxima.values())<1e-10,maxima
loss_error=0.;loss_cases=0
for r0 in np.geomspace(.5,1e12,25):
 for ratio in [1.,1+1e-12,1.01,2.,1e6]:
  r1=r0/ratio;a,b=map(mp.mpf,[r0,r1]);gamma=mp.mpf('1.4')
  expected=1/(a*a) if a==b else (gamma+1)/(gamma-1)*(a**(gamma-1)-b**(gamma-1))/(a**(gamma+1)-b**(gamma+1))
  actual=density_factor(r0,r1);loss_error=max(loss_error,float(abs(actual/expected-1)));assert actual>0 and actual==density_factor(r1,r0);loss_cases+=1
assert loss_error<1e-12,loss_error
reports={n:json.loads((ROOT/f'out/contour-pocket/{n}.json').read_text()) for n in ['continuous','continuous-fine','continuous-finer-start']}
for r in reports.values():
 for name,h in r['hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
pressure_runs={n:json.loads((ROOT/f'out/contour-pocket/{n}.json').read_text()) for n in ['pressure-loss','pressure-loss-small']}
for r in pressure_runs.values():
 assert r['success'] and r['message']=='A termination event occurred.'
 for name,h in r['hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
 assert abs(r['events'][-1]['work_residual_J'])<1e-12
large=pressure_runs['pressure-loss']['events'][-1];small=pressure_runs['pressure-loss-small']['events'][-1]
assert small['mass'][0]<large['mass'][0] and small['area_rate'][1]<0

a,b=reports['continuous'],reports['continuous-fine'];endpoint=max(abs(np.array(a['final_scaled_state'])-b['final_scaled_state']));assert endpoint<1e-6,endpoint
rebound_a=next(r for r in a['events'] if r['tag']=='rebound');rebound_b=next(r for r in b['events'] if r['tag']=='rebound')
rebound_relative=max(abs(rebound_a[k][0]/rebound_b[k][0]-1) for k in ['volume','mass','relative_density','pressure']);assert rebound_relative<1e-5,rebound_relative
max_work=max(abs(e['work_residual_J']) for r in reports.values() for e in r['events']);assert max_work<1e-12,max_work
report=dict(scope='diagnosis of the existing two-chamber model, not a physical validation of its extreme-compression rebound',kernel=dict(cases=len(rows),digits=80,maxima=maxima),pressure_loss=dict(cases=loss_cases,maximum_relative_error=loss_error,runs=pressure_runs),continuous=dict(endpoint_scaled_refinement_difference=float(endpoint),rebound_relative_refinement_difference=float(rebound_relative),maximum_work_residual_J=max_work,runs=reports),hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['tools/contour_pocket_audit.py','tools/contour_pocket_reference.py','tools/contour_sealed_reference.py','tools/contour_pocket_probe.cpp','build/ducklips_contour_pocket_probe']})
(ROOT/'reports/contour_pocket_reference.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(kernel=report['kernel'],endpoint_difference=endpoint,rebound_difference=rebound_relative,work=max_work),indent=2))
