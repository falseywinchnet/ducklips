#!/usr/bin/env python3
"""Independent quadrature and regression audit for continuous BGK slit loss."""
import json,subprocess,hashlib
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.integrate import quad_vec,quad
from scipy.interpolate import CubicSpline
from contour_pressure_flow_reference import Reference,SCALE,H,C0,C1
from contour_pressure_reference import oracle,unpack
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/slit-kinetic';BINARY=ROOT/'build/ducklips_contour_pressure_loss_probe'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(binary,args,rows,name):
 text=''.join(' '.join(map(str,row))+'\n' for row in rows);(OUT/(name+'-input.txt')).write_text(text)
 p=subprocess.run([str(binary),*args],input=text,text=True,capture_output=True);(OUT/(name+'-output.txt')).write_text(p.stdout);assert p.returncode==0,p.stderr;return p.stdout
reference=json.loads((ROOT/'reports/slit_kinetic_reference.json').read_text());table=reference['table'];ds=np.array([r['delta'] for r in table]);gs=np.array([r['conductance'] for r in table]);step=np.log(ds[1]/ds[0])
def low(d):return (-np.log(d)+1.5*(1-np.euler_gamma))/np.sqrt(np.pi)
def high(d):return d/6+1.0162+1.0653/d-(2.1354/d)/d
gs[0]=low(ds[0]);gs[-1]=high(ds[-1]);d=ds[-1]
spline=CubicSpline(np.log(ds),np.log(gs),bc_type=((1,-1/(np.sqrt(np.pi)*gs[0])),(1,(d/6-1.0653/d+4.2708/d**2)/gs[-1])))
def conductance(d):return low(d) if d<=1e-8 else high(d) if d>=100 else np.exp(spline(np.log(d)))
def matrix(ref,tol):
 def f(z):
  a,r,_,_=ref.field(z);thermal=350*np.sqrt(2/1.4)*r**.2;delta=a/.014*1.14*r*thermal/(2*1.86e-5);q=ref.flow(z)/a/SCALE
  return (1.14*thermal*.014/(r*conductance(delta))*np.outer(q,q)).ravel()
 parts=[]
 for a,b in zip(ref.nodes,ref.nodes[1:]):
  # Direct axial Gauss-Kronrod, independent of the native cubic mapping.
  value,error=quad_vec(f,a,b,epsabs=1e-5,epsrel=tol)
  parts.append(value)
 return sum(parts).reshape(5,5)*np.outer(SCALE,SCALE)
def main():
 OUT.mkdir(exist_ok=True,parents=True)
 rows=[r['input'] for r in json.loads((ROOT/'reports/contour_pressure_flow_reference.json').read_text())['cases']]
 native=[json.loads(s)['resistance'] for s in run(BINARY,['--batch'],rows,'loss').splitlines()];records=[]
 for row,values in zip(rows,native):
  ref=Reference(row);fine=matrix(ref,1e-9);coarse=matrix(ref,1e-7);actual=np.array(values).reshape(5,5)
  dimensional=np.outer(SCALE,SCALE);scale=max(np.diag(fine/dimensional));error=float(np.max(abs((actual-fine)/dimensional))/scale);change=float(np.max(abs((fine-coarse)/dimensional))/scale);eigen=float(min(np.linalg.eigvalsh(actual/dimensional))/scale)
  record=dict(input=row,matrix_scaled_error=error,reference_refinement=change,minimum_scaled_eigenvalue=eigen);records.append(record);print('matrix',json.dumps(record),flush=True)
  assert error<1e-7 and change<1e-8 and eigen>-1e-12,record
 checks=reference['interpolation_checks'];deltas=[r['delta'] for r in checks]+[1e-300,1e-12,1e-8,100.,200.,500.,1e6,1e300]
 native_g=np.array([float(s) for s in run(BINARY,['--conductance'],[[d] for d in deltas],'conductance').splitlines()]);kinetic_error=max(abs(native_g[:len(checks)]/np.array([r['reference'] for r in checks])-1));port_error=max(abs(native_g/np.array([conductance(d) for d in deltas])-1));assert kinetic_error<2e-5 and port_error<1e-12
 # Constant-density moving tip, analytic flux: Q=r*integral(-A_t) dz.
 contact_row=[-1e-6,3e-6,700.,700.];ref=Reference(contact_row);adot=np.array([-.0003,-.0003]);mdot=ref.D@adot;rate=np.r_[0.,-mdot[0],-sum(mdot),adot];r=ref.field(ref.hi)[1];thermal=350*np.sqrt(2/1.4)*r**.2
 def loss(z):
  a=ref.field(z)[0];q=r*.0003*(z-ref.lo);d=a/.014*1.14*r*thermal/(2*1.86e-5);return 1.14*thermal*.014/(r*conductance(d))*(q/a)**2
 cutoffs=[1e-5,1e-6,1e-7,1e-8,1e-9,0.];powers=[]
 for eps in cutoffs:
  nodes=[ref.lo+eps,*[z for z in ref.nodes if ref.lo+eps<z<ref.hi],ref.hi];powers.append(sum(quad(loss,a,b,epsabs=1e-15,epsrel=1e-9)[0] for a,b in zip(nodes,nodes[1:])))
 m=np.array(json.loads(run(BINARY,['--batch'],[contact_row],'contact').strip())['resistance']).reshape(5,5);power=float(rate@m@rate);assert abs(power/powers[-1]-1)<1e-7
 # Fresh 80-digit storage regression; historical reports remain immutable.
 storage_rows=[r['input'] for r in json.loads((ROOT/'reports/contour_pressure_reference.json').read_text())['cases_detail']]
 storage_binary=ROOT/'build/ducklips_contour_pressure_probe';states=[json.loads(s) for s in run(storage_binary,['--batch'],storage_rows,'storage').splitlines()];mp.mp.dps=80;errors=[]
 for i,(row,state) in enumerate(zip(storage_rows,states)):
  expected=oracle(row);v=float(expected[0]);scale=[v,v,v,139650*v,139650*v,139650*.003,139650*.003,v/139650,v/139650,v/139650,.003,.003,.003]
  errors.append(float(max(abs(mp.mpf(a)-b)/c for a,b,c in zip(unpack(state),expected,scale))))
  if i%12==11:print('storage',i+1,max(errors),flush=True)
 assert max(errors)<1e-10
 report=dict(scope='native local BGK slit loss; independent spatial quadrature and fresh storage regression, no coupled time integration or new speech',cases=records,native_conductance_max_relative_kinetic_error=float(kinetic_error),native_conductance_max_relative_port_error=float(port_error),contact=dict(input=contact_row,rate=rate.tolist(),cutoffs_m=cutoffs,truncated_loss_W=powers,native_full_loss_W=power),storage=dict(cases=len(errors),maximum_scaled_error=max(errors),precision_digits=80),hashes={str(p.relative_to(ROOT)):sha(p) for p in [BINARY,storage_binary,Path(__file__),ROOT/'tools/slit_kinetic_reference.py',ROOT/'tools/contour_pressure_flow_reference.py',ROOT/'tools/contour_pressure_reference.py',ROOT/'src/slit_bgk_table.inc',ROOT/'src/contour_pressure_loss.cpp',ROOT/'src/contour_pressure_flow.cpp',ROOT/'src/contour_pressure_storage.cpp',*sorted(OUT.glob('*-input.txt')),*sorted(OUT.glob('*-output.txt'))]})
 (ROOT/'reports/contour_pressure_loss_reference.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS independent loss, contact convergence, native conductance, 80-digit storage',flush=True)
if __name__=='__main__':main()
