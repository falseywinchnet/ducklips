#!/usr/bin/env python3
"""Independent ODE for continuous source, cyclic tract, and aperture radiation."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from continuous_outlet_column_reference import ColumnReference
from continuous_phonation_reference import YSCALE,VSCALE,STIFF,REST
from graph_continuous_glottis_reference import EDGES,B,C,G,NODAL,SERIES,PORT,XSCALE,INITIAL,XINITIAL
from radiation_boundary_reference import inverse,COEFF
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-radiation';DURATION=.002

def area(t,moving):
 if not moving:return 1e-4,0.
 return 1e-4*(1+.3*np.sin(2*np.pi*300*t)),3e-5*2*np.pi*300*np.cos(2*np.pi*300*t)

def continuous(moving,order,tolerance):
 ref=ColumnReference(order)
 def rhs(t,y):
  physical=y[:9]*YSCALE;x=y[9:16]*XSCALE;p=y[16]*.1
  a,adot=area(t,moving);h=inverse(a);hdot=.5*h*adot/a;r=1/(COEFF*h*h)
  pressure=r*(G[-1]*x[-1]-h*p)
  derivative=ref.rhs(y[:9],1400.);state=ref.configuration(physical[2:4],physical[:2])
  load=np.array([0.,0.,-x[PORT]/C[PORT],0.,0.])
  derivative[4:]+=np.linalg.solve(state['G']*VSCALE[:,None]*VSCALE[None,:],VSCALE*load)
  qdot=-B@(G*x[3:])+NODAL;qdot[PORT]+=physical[6]
  pdot=B.T@(x[:3]/C)+SERIES;pdot[-1]-=pressure
  return np.r_[derivative,np.r_[qdot,pdot]/XSCALE,pressure/.1,.5*hdot*p*p/1e-6,pressure*pressure/r/1e-6]
 sol=solve_ivp(rhs,(0.,DURATION),np.r_[INITIAL/YSCALE,XINITIAL/XSCALE,.004/.1,0.,0.],method='Radau',rtol=tolerance,atol=tolerance*.01,max_step=1e-5,dense_output=True)
 assert sol.success,sol.message
 return sol

def run(moving,n,oracle):
 dt=DURATION/n;lines=[f'3 4 {n} {PORT}',*[f'{i} {j}' for i,j in EDGES], ' '.join(map(str,INITIAL)),'1','80','1e-4 .004',' '.join(map(str,XINITIAL))]
 for k in range(n):lines.append(' '.join(map(str,np.r_[dt,1400.,2.,80.,area((k+.5)*dt,moving)[0],area((k+1)*dt,moving)[0],G[-1],G[-1],C,C,G,G,NODAL,SERIES])))
 name=f'{"moving" if moving else "fixed"}-{n}';fixture='\n'.join(lines)+'\n';(OUT/(name+'.input')).write_text(fixture)
 result=subprocess.run([str(ROOT/'build/ducklips_graph_continuous_glottis_probe'),'--radiation-fixture'],input=fixture,capture_output=True,text=True)
 (OUT/(name+'.txt')).write_text(result.stdout);(OUT/(name+'.stderr')).write_text(result.stderr)
 assert result.returncode==0,result.stderr
 data=np.loadtxt(OUT/(name+'.txt'));assert data.shape==(n,49),data.shape
 ledger=data[:,23:42];times=np.arange(1,n+1)*dt
 native=np.c_[data[:,:9]/YSCALE,data[:,9:16]/XSCALE,data[:,44]/.1]
 error=abs(native-oracle.sol(times).T[:,:17])
 ref=ColumnReference(48);s=ref.configuration(INITIAL[2:4],INITIAL[:2]);p=ref.fields(ref.t,INITIAL[2:4],INITIAL[:2])[3]
 energy=INITIAL[:2]@s['mass']-np.sum(s['w']*s['a']*p)+.5*INITIAL[4:]@s['G']@INITIAL[4:]+.5*(INITIAL[2:4]-REST)@STIFF@(INITIAL[2:4]-REST)
 energy+=.5*np.sum(XINITIAL[:3]**2/C)+.5*np.sum(G*XINITIAL[3:]**2)+.5*inverse(1e-4)*.004**2
 graph=data[:,9:16];energies=ledger[:,0]+.5*np.sum(graph[:,:3]**2/C,axis=1)+.5*np.sum(G*graph[:,3:]**2,axis=1)+data[:,45]
 work=np.cumsum(ledger[:,3:5].sum(axis=1)+ledger[:,11]+data[:,42]+data[:,47]);loss=np.cumsum(ledger[:,5:9].sum(axis=1)+data[:,46])
 residual=float(np.max(abs(energies-energy+loss-work)));assert residual<1e-13,residual
 mass0=s['mass'].sum()+XINITIAL[:3].sum();mass=ledger[:,1:3].sum(axis=1)+graph[:,:3].sum(axis=1)
 net=np.cumsum(dt*(ledger[:,17]+NODAL.sum()-data[:,22]));mass_error=float(np.max(abs(mass-mass0-net)));assert mass_error<1e-19,mass_error
 assert np.min(data[:,46])>=0 and np.max(abs(data[:,48]))<1e-18
 return dict(steps=n,endpoint_scaled_error=float(error[-1].max()),maximum_sampled_scaled_error=float(error.max()),
  joint_work_residual_J=residual,mass_residual_m3=mass_error,radiated_energy_J=float(data[:,46].sum()),
  independent_radiated_energy_J=float(oracle.y[18,-1]*1e-6),radiated_energy_error_J=float(abs(data[:,46].sum()-oracle.y[18,-1]*1e-6)),
  radiation_geometry_work_J=float(data[:,47].sum()),independent_geometry_work_J=float(oracle.y[17,-1]*1e-6),
  geometry_work_error_J=float(abs(data[:,47].sum()-oracle.y[17,-1]*1e-6)))

def main():
 OUT.mkdir(parents=True,exist_ok=True);cases=[]
 for moving in [False,True]:
  fine=continuous(moving,48,2e-10);coarse=continuous(moving,32,2e-9);times=np.linspace(0,DURATION,6401)
  delta=float(np.max(abs(fine.sol(times)-coarse.sol(times))));assert delta<2e-7,delta
  np.savez_compressed(OUT/f'continuous-{moving}.npz',times=times,fine=fine.sol(times),coarse=coarse.sol(times))
  runs=[]
  for n in [200,400,800,1600]:
   record=run(moving,n,fine);runs.append(record);print(json.dumps(dict(moving=moving,**record)),flush=True)
  ratios=[a['endpoint_scaled_error']/b['endpoint_scaled_error'] for a,b in zip(runs,runs[1:])]
  radiated_ratios=[a['radiated_energy_error_J']/b['radiated_energy_error_J'] for a,b in zip(runs,runs[1:])]
  assert 3.5<ratios[-1]<4.5,ratios
  assert 3.5<radiated_ratios[-1]<4.5,radiated_ratios
  cases.append(dict(moving_aperture=moving,reference_refinement_scaled_error=delta,runs=runs,endpoint_refinement_ratios=ratios,radiated_energy_refinement_ratios=radiated_ratios))
 paths=[Path(__file__),ROOT/'tools/continuous_outlet_column_reference.py',ROOT/'tools/continuous_phonation_reference.py',ROOT/'tools/radiation_boundary_reference.py',
 ROOT/'tools/graph_continuous_glottis_probe.cpp',ROOT/'src/graph_continuous_glottis.cpp',ROOT/'include/ducklips/graph_continuous_glottis.hpp',ROOT/'src/radiation_boundary.cpp',ROOT/'include/ducklips/radiation_boundary.hpp',
 ROOT/'build/ducklips_graph_continuous_glottis_probe',*sorted(OUT.glob('*.input')),*sorted(OUT.glob('*.txt')),*sorted(OUT.glob('*.stderr')),*sorted(OUT.glob('*.npz'))]
 report=dict(scope='Open continuous gas/tissue source plus cyclic fixed tract and passive one-pole radiation, checked against independent time-continuous equations. Changing aperture work and radiated energy are separate. No speech or high-frequency aperture calibration.',cases=cases,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
 (ROOT/'reports/continuous_radiation_reference.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
