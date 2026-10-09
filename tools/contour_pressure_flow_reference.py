#!/usr/bin/env python3
"""Independent local-continuity and axial-inertia audit."""
import hashlib,json,math,subprocess
from pathlib import Path
import numpy as np
from scipy.integrate import quad,quad_vec
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/contour-flow'
BULK=139650.;GAMMA=1.4;H=.003;C0=.00125;C1=.00275;ELL=C1-C0
SCALE=np.array([1.,1.,1.,H,H])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

class Reference:
 def __init__(self,row):
  self.a,self.b,self.u0,self.u1=row;self.lo=0. if self.a>0 else C0-ELL*self.a/(self.b-self.a);self.hi=H if self.b>0 else C0+ELL*self.a/(self.a-self.b)
  self.nodes=sorted(set([self.lo,self.hi,*[z for z in [C0,C1] if self.lo<z<self.hi]]))
  def matrices(z):
   area,r,slope,phi=self.field(z);outer=np.outer(phi,phi);return np.r_[(area*slope*outer).ravel()/1e-14,(r*outer).ravel()/H]
  value=self.integral(matrices,self.lo,self.hi);self.C=value[:4].reshape(2,2)*1e-14;self.D=value[4:].reshape(2,2)*H
  right=np.c_[np.array([[1,-1,0],[0,1,-1]]),-self.D];self.E=np.linalg.solve(self.C,right)
  self.closed=[j for j,v in [(0,self.a),(2,self.b)] if v<=0];self.E[:,self.closed]=0
  self.base=1.14*H*H/self.integral(lambda z:np.array([self.field(z)[0]*self.field(z)[1]]),self.lo,self.hi)[0]
 def field(self,z):
  t=min(1.,max(0.,(z-C0)/ELL));phi=np.array([1-t,t]);a=self.a*(1-t)+self.b*t;mu=self.u0*(1-t)+self.u1*t;r=(1+(GAMMA-1)*mu/BULK)**(1/(GAMMA-1));return a,r,r**(2-GAMMA)/BULK,phi
 def integral(self,f,lo,hi,tol=1e-12):
  nodes=[lo]+[z for z in self.nodes if lo<z<hi]+[hi];return sum((quad_vec(f,a,b,epsabs=tol,epsrel=tol)[0] for a,b in zip(nodes,nodes[1:])),np.zeros_like(f((lo+hi)/2))) if hi>lo else np.zeros_like(f(lo))
 def flow(self,z):
  right=z>(self.lo+self.hi)/2
  def terms(x):
   a,r,slope,phi=self.field(x);return np.r_[a*slope*phi/1e-14,r*phi/H]
  value=self.integral(terms,z if right else self.lo,self.hi if right else z)
  result=(value[:2]*1e-14)@self.E;result[3:]+=value[2:]*H;result*=1 if right else -1;result[2 if right else 0]+=1;result[self.closed]=0;return result
 def matrix(self,tol):
  def f(z):
   a,r,_,_=self.field(z);q=self.flow(z)/SCALE;return (1.14*np.outer(q,q)/(r*a)/self.base).ravel()
  return self.integral(f,self.lo,self.hi,tol).reshape(5,5)*self.base*np.outer(SCALE,SCALE)

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 rows=[[3e-6,3e-6,700,700],[1e-6,4e-6,-100,900],[5e-6,1e-6,1500,-500],[-4.8e-6,1e-6,-100,900],[-5e-6,1e-6,-100,900],[-5.2e-6,1e-6,-100,900],[3e-6,-.5e-6,700,-100],[3e-6,-2e-6,900,-500],[0,3e-6,700,200],[3e-6,0,900,-100],[2e-6,5e-6,-20000,10000],[-1e-6,3e-6,-20000,10000]]
 inp=OUT/'input.txt';inp.write_text(''.join(' '.join(map(str,row))+'\n' for row in rows));binary=ROOT/'build/ducklips_contour_pressure_flow_probe';run=subprocess.run([str(binary),'--batch'],input=inp.read_text(),capture_output=True,text=True);assert run.returncode==0,run.stderr
 output=OUT/'output.jsonl';output.write_text(run.stdout);native=[json.loads(line) for line in run.stdout.splitlines()];assert len(native)==len(rows);records=[]
 for row,s in zip(rows,native):
  ref=Reference(row);matrix=ref.matrix(1e-11);coarse=ref.matrix(1e-9);scale=ref.base*np.outer(SCALE,SCALE)
  err=float(np.max(abs(np.array(s['inertia']).reshape(5,5)-matrix)/scale));refinement=float(np.max(abs(matrix-coarse)/scale));field_error=0
  for p in s['samples']:field_error=max(field_error,float(max(abs(np.array(p['flow'])-ref.flow(p['z']))/SCALE)))
  start=max(C0,ref.lo);end=min(C1,ref.hi);weak=ref.integral(lambda z:ref.flow(z)/SCALE,start,end)/ELL;target=np.array([0.,1.,0.,0.,0.]);weak_error=float(max(abs(weak-target)))
  normalized=np.array(s['inertia']).reshape(5,5)/scale;eigen=float(min(np.linalg.eigvalsh(normalized)));assert err<2e-8 and refinement<1e-9 and field_error<1e-9 and weak_error<1e-9 and eigen>-1e-11,(row,err,refinement,field_error,weak_error,eigen)
  records.append(dict(input=row,matrix_scaled_error=err,reference_refinement=refinement,flow_scaled_error=field_error,weak_continuity_error=weak_error,minimum_scaled_eigenvalue=eigen));print(json.dumps(records[-1]),flush=True)
 # A moving contact tip has finite axial kinetic energy but logarithmic
 # divergence of the no-slip Poiseuille dissipation. Evaluate cutoffs only
 # as an asymptotic audit; they are never inserted into the source model.
 row=[-1e-6,3e-6,700.,700.];ref=Reference(row);area_rate=np.array([-.0003,-.0003]);mass_rate=ref.D@area_rate;rate=np.r_[0.,-mass_rate[0],-sum(mass_rate),area_rate];assert max(abs(ref.E@rate))<1e-7
 slope=(ref.b-ref.a)/ELL;tip_speed=-sum(ref.field(ref.lo)[3]*area_rate)/slope;coefficient=12*1.86e-5*.014**2*tip_speed**2/abs(slope)
 cutoffs=np.array([1e-5,1e-6,1e-7,1e-8,1e-9]);loss=[];speeds=[]
 for eps in cutoffs:
  def integrand(z):
   a,r,_,_=ref.field(z);q=ref.flow(z)@rate;return 12*1.86e-5*.014**2*q*q/(r*r*a**3)
  value=0.
  for a,b in zip([ref.lo+eps,*[z for z in ref.nodes if ref.lo+eps<z<ref.hi]],[*[z for z in ref.nodes if ref.lo+eps<z<ref.hi],ref.hi]):value+=quad(integrand,a,b,epsabs=1e-15,epsrel=1e-10)[0]
  loss.append(value);z=ref.lo+eps;a,r,_,_=ref.field(z);speeds.append(float((ref.flow(z)@rate)/(r*a)))
 native_coefficient=float(area_rate@np.array(native[-1]['contact_loss']).reshape(2,2)@area_rate);assert abs(native_coefficient/coefficient-1)<1e-13
 slopes=np.diff(loss)/np.log(10);assert abs(slopes[-1]/coefficient-1)<1e-5 and abs(speeds[-1]/tip_speed-1)<1e-6
 report=dict(scope='strong local continuity and axial fluid inertia; contact loss divergence diagnosed, no momentum advance or speech',cases=records,contact=dict(input=row,rate=rate.tolist(),tip_speed=tip_speed,near_tip_speed=speeds,cutoffs_m=cutoffs.tolist(),truncated_loss_W=loss,logarithmic_slopes=slopes.tolist(),expected_logarithmic_coefficient_W=coefficient,native_logarithmic_coefficient_W=native_coefficient),hashes={str(p.relative_to(ROOT)):sha(p) for p in [inp,output,binary,ROOT/'include/ducklips/contour_pressure_flow.hpp',ROOT/'src/contour_pressure_flow.cpp',ROOT/'include/ducklips/contour_pressure_storage.hpp',ROOT/'src/contour_pressure_storage.cpp',ROOT/'tools/contour_pressure_flow_probe.cpp',Path(__file__)]})
 (ROOT/'reports/contour_pressure_flow_reference.json').write_text(json.dumps(report,indent=2)+'\n');print('contact',json.dumps(report['contact']),flush=True)
if __name__=='__main__':main()
