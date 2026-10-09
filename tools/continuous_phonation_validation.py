#!/usr/bin/env python3
"""Cross-check native onset against a separate smooth continuous-field ODE."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
from continuous_phonation_reference import Reference,YSCALE,VSCALE,REST,MASS,DAMP
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-phonation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def native(rows,name,binary,args=()):
 inp='\n'.join(' '.join(map(str,row)) for row in rows)+'\n';(OUT/(name+'.input')).write_text(inp);p=subprocess.run([str(ROOT/binary),*args],input=inp,text=True,capture_output=True,check=True);(OUT/(name+'.jsonl')).write_text(p.stdout);assert not p.stderr,p.stderr;result=[json.loads(s) for s in p.stdout.splitlines()];assert len(result)==len(rows) and all(s['accepted'] for s in result);return result

def main():
 ref=Reference();paths=[];static=[]
 configurations=[[3e-6,4e-6,300.,-100.],[3.7e-6,4.2e-6,120.,-3.5],[4.5e-6,2e-6,900.,-200.]]
 results=native(configurations,'constitutive','build/ducklips_contour_profile_limit_probe',['linear'])
 for row,s in zip(configurations,results):
  r=ref.configuration(np.array(row[:2]),np.array(row[2:]));G=np.array(s['inertia']).reshape(5,5);R=np.array(s['resistance']).reshape(5,5);reference_R=r['R'].copy();reference_R[3:,3:]-=np.diag(DAMP);errors={}
  for name,a,b in [('inertia',G,r['L']),('resistance',R,reference_R)]:
   scale=np.sqrt(np.diag(b));error=float(np.max(np.abs(a-b)/(scale[:,None]*scale[None,:])));assert error<1e-7,(row,name,error);errors[name]=error
  static.append(dict(input=row,maximum_scaled_errors=errors))
 y=np.array([.3,-.1,3.,4.,1.5,1.0,.7,.2,-.1]);target=ref.rhs(y,700);fine=Reference(48).rhs(y,700);reference_error=float(np.max(np.abs(fine-target)/np.maximum(np.abs(fine),100.)));assert reference_error<1e-9
 steps=[1e-9,5e-10,2.5e-10];runs=native([[1,dt,700,REST,*(y*YSCALE)] for dt in steps],'tangent','build/ducklips_continuous_phonation_step_probe');rates=[];tangent=[]
 for dt,s in zip(steps,runs):
  rate=(np.array(s['final'])-s['initial'])/YSCALE/dt;rates.append(rate);tangent.append(dict(dt=dt,maximum_scaled_error=float(np.max(np.abs(rate-fine)/np.maximum(np.abs(fine),100.))),**s))
 extrapolated=2*rates[-1]-rates[-2];error=float(np.max(np.abs(extrapolated-fine)/np.maximum(np.abs(fine),100.)));assert error<1e-5,error
 threshold=[]
 for order,h in [(24,.001),(48,.0005)]:
  r=Reference(order)
  def growth(p):
   _,state,_=r.equilibrium(p,[.12,-.003,3.7,4.2,1.3]);return max(np.linalg.eigvals(r.jacobian(state,p,h)).real)
  pressure=brentq(growth,700,1050,xtol=1e-5);_,state,_=r.equilibrium(pressure,[.12,-.003,3.7,4.2,1.3]);eigen=np.linalg.eigvals(r.jacobian(state,pressure,h));mode=eigen[np.argmax(eigen.real)];threshold.append(dict(quadrature_order=order,jacobian_step=h,pressure_Pa=pressure,frequency_Hz=float(mode.imag/(2*np.pi))))
 assert abs(threshold[0]['pressure_Pa']-threshold[1]['pressure_Pa'])<.01
 modes=[]
 for m in json.loads((OUT/'modal-inputs.json').read_text()):
  A=np.array(m['jacobian']);eigen,V=np.linalg.eig(A);j=max((j for j in range(9) if eigen[j].imag>0),key=lambda j:eigen[j].real);left=np.linalg.inv(V)[j];records=[];growth=[]
  for n in [200,400,800]:
   path=OUT/f"mode-{m['pressure']}-{n}.jsonl";minus,plus=map(json.loads,path.read_text().splitlines());assert minus['accepted'] and plus['accepted'];d0=(np.array(plus['initial'])-minus['initial'])/YSCALE;d1=(np.array(plus['final'])-minus['final'])/YSCALE;rate=np.log((left@d1)/(left@d0))/.002;expected=expm(A*.002)@d0;endpoint_error=float(np.max(np.abs(d1-expected))/np.max(np.abs(d0)));growth.append(rate);records.append(dict(steps=n,dt_s=.002/n,growth_s=float(rate.real),frequency_Hz=float(rate.imag/(2*np.pi)),relative_endpoint_error=endpoint_error,roundoff_limited_steps=[minus['roundoff_limited_steps'],plus['roundoff_limited_steps']],maximum_step_work_residual=max(minus['maximum_step_work_residual'],plus['maximum_step_work_residual'])))
  ratios=[a['relative_endpoint_error']/b['relative_endpoint_error'] for a,b in zip(records,records[1:])];assert all(1.8<r<2.2 for r in ratios);extrapolated=2*growth[-1]-growth[-2];assert abs(extrapolated-eigen[j])<.02;assert all(r['roundoff_limited_steps'][0]>0 and r['roundoff_limited_steps'][1]>0 for r in records)
  modes.append(dict(pressure_Pa=m['pressure'],reference_growth_s=float(eigen[j].real),reference_frequency_Hz=float(eigen[j].imag/(2*np.pi)),native=records,endpoint_refinement_ratios=ratios,richardson_growth_s=float(extrapolated.real),richardson_frequency_Hz=float(extrapolated.imag/(2*np.pi))))
 report=dict(scope='open linear-span model in the symmetric tissue subspace; model onset, not measured physiological threshold',base_commit='dc04d2abbf49ba076774c0602d54144cd2bfc2bc',constitutive=static,independent_ode_quadrature_refinement=reference_error,tangent=dict(steps=tangent,richardson_maximum_scaled_error=error),threshold=threshold,modes=modes)
 paths=[Path(__file__),ROOT/'tools/continuous_phonation_reference.py',ROOT/'tools/continuous_phonation_modal_experiment.py',ROOT/'tools/continuous_phonation_step_probe.cpp',ROOT/'tools/continuous_phonation_controls.cpp',ROOT/'build/ducklips_continuous_phonation_step_probe',ROOT/'build/ducklips_contour_profile_limit_probe',ROOT/'src/continuous_glottis.cpp',ROOT/'include/ducklips/continuous_glottis.hpp',*OUT.glob('mode-*.jsonl'),*OUT.glob('mode-*.input'),OUT/'modal-inputs.json',OUT/'constitutive.input',OUT/'constitutive.jsonl',OUT/'tangent.input',OUT/'tangent.jsonl']
 report['hashes']={str(p.relative_to(ROOT)):sha(p) for p in paths};(ROOT/'reports/continuous_phonation_onset.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS independent field matrices, nonlinear tangent, onset refinement, and native modal growth/decay')
if __name__=='__main__':main()
