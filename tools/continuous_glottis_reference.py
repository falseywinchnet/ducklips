#!/usr/bin/env python3
"""Analytic linear-limit and time-refinement audit of coupled contour motion."""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import expm,eigvalsh
from contour_pressure_flow_reference import Reference
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-glottis'
S=np.array([1e-8]*3+[1e-4]*4);TIME=.001;ENERGY=1e-6;WIDTH=.014

def vector(r):return np.r_[np.array(r['mass'])/1e-8,np.array(r['position'])/1e-4,np.array(r['flow'])/1e-5,np.array(r['velocity'])/.1]
def main():
 ref=Reference([3e-6,3e-6,0.,0.]);fluid=ref.matrix(1e-12);mapping=np.zeros((5,7));mapping[:3,:3]=np.eye(3);mapping[3,3:5]=WIDTH;mapping[4,5:]=WIDTH
 tissue_mass=np.array([1.25e-4/.9,1.25e-4/(.9*1.17),2.5e-5/.9,2.5e-5/(.9*1.17)]);stiffness=np.array([80*.9,80*.9*1.17,8*.9,8*.9*1.17]);coupling=np.array([25*.9,25*.9*1.17])
 metric=mapping.T@fluid@mapping;metric[3:,3:]+=np.diag(tissue_mass)
 mass_map=np.array([[1.,-1.,0.],[0.,1.,-1.]]);area_map=mapping[3:,3:]
 compress=np.c_[mass_map,-ref.D@area_map];potential=compress.T@np.linalg.solve(ref.C,compress);potential[3:,3:]+=np.diag(stiffness)
 for side,k in enumerate(coupling):
  i,j=3+side,5+side;potential[i,i]+=k;potential[j,j]+=k;potential[i,j]-=k;potential[j,i]-=k
 M=metric*np.outer(S,S)/(TIME**2*ENERGY);K=potential*np.outer(S,S)/ENERGY
 generator=np.block([[np.zeros((7,7)),np.eye(7)],[-np.linalg.solve(M,K),np.zeros((7,7))]])
 base_x=np.array([1e-6,-2e-6,-1e-6,3e-6]);base_u=np.array([100.,-30.]);base_v=np.array([.001,-.002,.003,-.001]);base_flow=np.array([1e-6,-2e-6,3e-6]);records=[]
 for amplitude in [.001,.0001]:
  x=base_x*amplitude;mass=ref.C@(base_u*amplitude)+ref.D@(area_map@x);q=np.r_[mass[0],0.,-mass[1],x]/S;v=np.r_[base_flow*amplitude,base_v*amplitude]*TIME/S
  exact=expm(generator*(1e-5/TIME))@np.r_[q,v];qe,ve=exact[:7]*S,exact[7:]*S/TIME
  exact_vector=np.r_[(ref.C@np.zeros(2)+np.array([6e-9,3e-9])+mass_map@qe[:3])/1e-8,qe[3:]/1e-4,ve[:3]/1e-5,ve[3:]/.1]
  for n in [100,200,400,800,1600]:
   path=OUT/f'linear-{n}-{amplitude}.json';r=json.loads(path.read_text());error=float(np.max(abs(vector(r)-exact_vector))/amplitude)
   records.append(dict(amplitude=amplitude,steps=n,maximum_normalized_error_per_amplitude=error,maximum_energy_residual_J=r['maximum_energy_residual_J']));print('linear',records[-1])
 assert records[-1]['maximum_normalized_error_per_amplitude']<1e-5
 temporal=[]
 for mode in ['free','driven','lower','upper']:
  steps=[1600,3200,6400] if mode=='lower' else [800,1600,3200]
  runs=[json.loads((OUT/f'{mode}-{n}.json').read_text()) for n in steps];differences=[float(np.linalg.norm(vector(a)-vector(b))) for a,b in zip(runs,runs[1:])];ratio=differences[0]/differences[1]
  temporal.append(dict(mode=mode,steps=steps,normalized_endpoint_differences=differences,ratio=ratio,runs=runs));print('temporal',mode,ratio)
  assert 3.3<ratio<4.7,(mode,ratio)
  assert max(r['maximum_energy_residual_J'] for r in runs)<1e-12
 long=json.loads((OUT/'driven-1000-long.json').read_text())
 report=dict(scope='continuous fluid/tissue dynamics within fixed port topology; analytic small-amplitude limit and temporal refinement; not phonation or new speech',linear_limit=dict(reference='independent continuous-field quadrature, analytic ambient compliance/stiffness, exact matrix exponential of 14-state linear ODE',metric=M.tolist(),stiffness=K.tolist(),frequencies_Hz=(np.sqrt(np.maximum(0,eigvalsh(K,M)))/(2*np.pi*TIME)).tolist(),cases=records),temporal=temporal,driven_1ms=long,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'tools/contour_pressure_flow_reference.py',ROOT/'src/continuous_glottis.cpp',ROOT/'include/ducklips/continuous_glottis.hpp',ROOT/'src/gas_storage.cpp',ROOT/'src/contour_pressure_storage.cpp',ROOT/'build/ducklips_continuous_glottis_probe',*sorted(OUT.glob('linear-*.json')),*[OUT/f'{m}-{n}.json' for m in ['free','driven','lower','upper'] for n in [100,200,400,800,1600,3200]],OUT/'lower-6400.json',OUT/'driven-1000-long.json']})
 (ROOT/'reports/continuous_glottis_reference.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS analytic linear limit and four motion refinements')
if __name__=='__main__':main()
