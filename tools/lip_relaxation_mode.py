#!/usr/bin/env python3
"""Find a damped lip mode without subtracting large stiff-memory matrices."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import root
from planar_source_reference import B,EDGES
from slit_viscous_reference import spectrum
from slit_thermal_reference import spectrum as thermal_spectrum
from radiation_boundary_reference import inverse,COEFF
from relaxation_graph_reference import fixture
from radau5_reference import A as A5
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/lip-relaxation-mode'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 OUT.mkdir(parents=True,exist_ok=True);snapshot=ROOT/'out/planar-joint-fine/opening-1-2048.json';record=json.loads(snapshot.read_text())['steps'][-1]
 c=np.array(record['compliance']);g=np.array(record['conductance']);assert np.min(g)>0;sqrtc=np.sqrt(c)
 weights,poles=spectrum();tw,tp=thermal_spectrum();no=[];bo=[];nm=[];bm=[];resistance=np.zeros(15)
 for site in record['thermal_sites']:
  owner=site['node'];no.extend([owner]*24);nm.extend(zip(.4*tw*c[owner],tp*2.3e-5/site['half_gap']**2))
 for site in record['viscous_sites']:
  owner=site['branch'];inertia=site['segment_inertance'];gap=site['half_gap'];assert inertia>0;rate=(1.86e-5/1.14)/gap**2
  bo.extend([owner]*24);bm.extend(zip(inertia*weights/poles,poles*rate));resistance[owner]+=3*inertia*rate
 h=inverse(record['radiation_state']['area']);bo.append(6);bm.append((1/h,1/(COEFF*h)))
 no=np.array(no);bo=np.array(bo);nm=np.array(nm);bm=np.array(bm);other=np.array([j for j in range(14) if j!=5])
 def operator(s):
  cap=c.astype(complex);np.add.at(cap,no,nm[:,0]*nm[:,1]/(s+nm[:,1]))
  impedance=s/g+resistance+0j;np.add.at(impedance,bo,bm[:,0]*bm[:,1]*s/(s+bm[:,1]))
  return (np.diag(s*cap)+(B/impedance)@B.T)/sqrtc[:,None]/sqrtc[None,:],impedance
 def secular(x):
  matrix,_=operator(x[0]+1j*x[1]);value=matrix[5,5]-matrix[5,other]@np.linalg.solve(matrix[np.ix_(other,other)],matrix[other,5]);return [value.real,value.imag]
 lam,v=eigh((B*g)@B.T/sqrtc[:,None]/sqrtc[None,:]);initial_frequency=np.sqrt(lam[np.argmax(v[5]**2)])
 solution=root(secular,[-1000,initial_frequency],tol=1e-11);assert solution.success and max(abs(np.array(secular(solution.x))))<1e-6
 pole=solution.x[0]+1j*solution.x[1];matrix,impedance=operator(pole);shape=np.zeros(14,dtype=complex);shape[5]=1;shape[other]=-np.linalg.solve(matrix[np.ix_(other,other)],matrix[other,5]);pressure=shape/sqrtc;pressure/=pressure[5]
 flow=(B.T@pressure)/impedance;thermal=nm[:,0]*nm[:,1]*pressure[no]/(pole+nm[:,1]);viscous=bm[:,0]*bm[:,1]*flow[bo]/(pole+bm[:,1])
 mode=np.r_[c*pressure,flow/g,thermal,viscous];metric=np.sqrt(np.r_[1/c,g,1/nm[:,0],1/bm[:,0]])
 charge_balance=pole*(c*pressure)+B@flow;np.add.at(charge_balance,no,pole*thermal)
 pressure_balance=pole*flow/g+resistance*flow-B.T@pressure;np.add.at(pressure_balance,bo,pole*viscous)
 assert max(abs(charge_balance))<1e-14 and max(abs(pressure_balance))<1e-11
 energy=.5*abs(mode*metric)**2;partition=dict(node=float(sum(energy[:14])),branch=float(sum(energy[14:29])),thermal=float(sum(energy[29:29+len(no)])),viscous_and_radiation=float(sum(energy[29+len(no):])))
 edges=[(u,14 if v<0 else v) for u,v in EDGES];binary=ROOT/'build/ducklips_relaxation_graph_probe';cases=[]
 for order,a,flag in [(3,np.array([[5/12,-1/12],[.75,.25]]),'--fixture'),(5,A5,'--fixture-fifth')]:
  b=a[-1]
  for n in [256,512,1024,2048]:
   dt=.001/n;steps=n//4;duration=steps*dt;factor=1+dt*pole*(b@np.linalg.solve(np.eye(len(b))-dt*pole*a,np.ones(len(b))));effective=np.log(factor)/dt
   path=OUT/f'order-{order}-{n}.input';path.write_text(fixture(c,g,resistance,edges,no,bo,nm,bm,mode.real,np.zeros(14),np.zeros(15),dt,steps)+'\n')
   result=subprocess.run([str(binary),flag],input=path.read_text(),capture_output=True,text=True);assert result.returncode==0,result.stderr;target=path.with_suffix('.txt');target.write_text(result.stdout);data=np.fromstring(result.stdout,sep=' ')
   expected=(mode*factor**steps).real;exact=(mode*np.exp(pole*duration)).real
   discrete=float(max(abs((data[:len(mode)]-expected)*metric)));assert discrete<1e-14,discrete
   cases.append(dict(order=order,subdivisions_per_ms=n,duration_s=duration,native_to_discrete_mode_energy_coordinate_error=discrete,native_to_continuous_mode_energy_coordinate_error=float(max(abs((data[:len(mode)]-exact)*metric))),extra_numerical_decay_s=float(pole.real-effective.real),frequency_error_Hz=float((effective.imag-pole.imag)/(2*np.pi)),relative_amplitude_error=float(abs(factor**steps)/abs(np.exp(pole*duration))-1),last_lip_pressure_error_Pa=float((data[5]-exact[5])/c[5]),hashes={str(p.relative_to(ROOT)):sha(p) for p in [path,target]}))
 paths=[Path(__file__),snapshot,binary,ROOT/'tools/relaxation_graph_reference.py',ROOT/'tools/radau5_reference.py',ROOT/'tools/slit_viscous_reference.py',ROOT/'tools/slit_thermal_reference.py',ROOT/'tools/radiation_boundary_reference.py',ROOT/'tools/planar_source_reference.py',ROOT/'src/relaxation_graph.cpp',ROOT/'include/ducklips/relaxation_graph.hpp']
 report=dict(scope='Frozen actual opening endpoint, all tract viscous/thermal/radiation memories, source port unforced. Rational pressure operator avoids cancellation of fast-memory rates. A mode and native time-response diagnostic, not a complete moving-source spectrum or a proof of audible effect.',frequency_Hz=float(pole.imag/(2*np.pi)),physical_decay_s=float(-pole.real),secular_residual=float(max(abs(np.array(secular(solution.x))))),charge_balance_residual_m3_s=float(max(abs(charge_balance))),pressure_balance_residual_Pa=float(max(abs(pressure_balance))),complex_pressure_shape=[[float(x.real),float(x.imag)] for x in pressure],complex_mode_energy_J=float(sum(energy)),energy_partition_J=partition,cases=cases,hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})
 (ROOT/'reports/lip_relaxation_mode.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['hashes','complex_pressure_shape']},indent=2))
if __name__=='__main__':main()
