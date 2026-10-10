#!/usr/bin/env python3
"""Independent continuous source/tract/radiation/thermal-storage/viscous ODE."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix
from continuous_outlet_column_reference import ColumnReference
from continuous_phonation_reference import YSCALE,VSCALE,STIFF,REST
from graph_continuous_glottis_reference import EDGES,B,C,G,NODAL,SERIES,PORT,XSCALE,INITIAL,XINITIAL
from radiation_boundary_reference import inverse,COEFF
from slit_thermal_reference import spectrum
from slit_viscous_reference import spectrum as viscous_spectrum
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-viscous';DURATION=.0005
W,POLE=spectrum();K=.4*W;GAPS=np.array([.001,.002,.0005]);RATES=2.3e-5/GAPS[:,None]**2*POLE
ZSCALE=C[:,None]*700*np.sqrt(K);H=inverse(1e-4);R=1/(COEFF*H*H)

VW,VP=viscous_spectrum();VBRANCH=np.array([1,1,3]);VL=np.array([.4/G[1],.6/G[1],.8/G[3]]);VGAP=np.array([.001,.002,.0005]);VRATE=(1.86e-5/1.14)/VGAP[:,None]**2*VP
VRESIST=VL[:,None]*(1.86e-5/1.14)/VGAP[:,None]**2*VW;VSTATIC=3*VL*(1.86e-5/1.14)/VGAP**2
VZSCALE=.001/np.sqrt(VP/(VL[:,None]*VW))

def continuous(order,tolerance):
 ref=ColumnReference(order)
 def rhs(t,y):
  physical=y[:9]*YSCALE;x=y[9:16]*XSCALE;p=y[16]*.1;z=y[17:89].reshape(3,24)*ZSCALE
  pressure=R*(G[-1]*x[-1]-H*p)
  derivative=ref.rhs(y[:9],1400.);state=ref.configuration(physical[2:4],physical[:2])
  load=np.array([0.,0.,-x[PORT]/C[PORT],0.,0.])
  derivative[4:]+=np.linalg.solve(state['G']*VSCALE[:,None]*VSCALE[None,:],VSCALE*load)
  difference=K*x[:3,None]-z;zdot=RATES*difference
  qdot=-B@(G*x[3:])+NODAL-zdot.sum(axis=1);qdot[PORT]+=physical[6]
  pdot=B.T@(x[:3]/C)+SERIES;pdot[-1]-=pressure
  heat=np.sum(RATES*difference**2/(C[:,None]*K))
  vm=y[91:163].reshape(3,24)*VZSCALE;vf=G[VBRANCH]*x[3+VBRANCH];vdot=VRESIST*vf[:,None]-VRATE*vm
  for i,branch in enumerate(VBRANCH):pdot[branch]-=VSTATIC[i]*vf[i]+vdot[i].sum()
  vheat=np.sum(VSTATIC*vf**2)+np.sum(vdot**2/VRESIST)
  return np.r_[derivative,np.r_[qdot,pdot]/XSCALE,pressure/.1,(zdot/ZSCALE).ravel(),heat/1e-6,pressure*pressure/R/1e-6,(vdot/VZSCALE).ravel(),vheat/1e-6]
 sparsity=np.zeros((164,164),dtype=bool);sparsity[:9,:9]=True;sparsity[:9,9+PORT]=True
 sparsity[9:16,:17]=True;sparsity[16,15:17]=True
 for i in range(3):
  ids=np.arange(17+24*i,17+24*(i+1));sparsity[9+i,ids]=True;sparsity[ids,9+i]=True;sparsity[ids,ids]=True
 sparsity[89,9:12]=True;sparsity[89,17:89]=True;sparsity[90,15:17]=True
 for i,branch in enumerate(VBRANCH):
  ids=np.arange(91+24*i,91+24*(i+1));sparsity[12+branch,ids]=True;sparsity[ids,12+branch]=True;sparsity[ids,ids]=True
  sparsity[163,12+branch]=True;sparsity[163,ids]=True
 sol=solve_ivp(rhs,(0.,DURATION),np.r_[INITIAL/YSCALE,XINITIAL/XSCALE,.004/.1,np.zeros(147)],method='Radau',rtol=tolerance,atol=tolerance*.01,max_step=5e-6,dense_output=True,jac_sparsity=csr_matrix(sparsity))
 assert sol.success,sol.message
 return sol

def run(n,oracle):
 dt=DURATION/n;lines=[f'3 4 {n} {PORT}',*[f'{i} {j}' for i,j in EDGES],' '.join(map(str,INITIAL)),'1','80','1e-4 .004',' '.join(map(str,XINITIAL)),
 ' '.join(map(str,np.c_[C,GAPS].ravel())),' '.join(['0']*72),'3',*[' '.join(map(str,[VBRANCH[i],VL[i],G[VBRANCH[i]],VGAP[i]])) for i in range(3)],' '.join(['0']*72)]
 for k in range(n):lines.append(' '.join(map(str,np.r_[dt,1400.,2.,80.,1e-4,1e-4,G[-1],G[-1],np.c_[C,GAPS].ravel(),np.c_[VL,G[VBRANCH],VGAP].ravel(),C,C,G,G,NODAL,SERIES])))
 path=OUT/f'fixed-{n}.input';path.write_text('\n'.join(lines)+'\n')
 result=subprocess.run([str(ROOT/'build/ducklips_graph_continuous_glottis_probe'),'--viscous-fixture'],input=path.read_text(),capture_output=True,text=True)
 path.with_suffix('.txt').write_text(result.stdout);path.with_suffix('.stderr').write_text(result.stderr);assert result.returncode==0,result.stderr
 data=np.loadtxt(path.with_suffix('.txt'));assert data.shape==(n,202),data.shape
 ledger=data[:,23:42];times=np.arange(1,n+1)*dt
 native=np.c_[data[:,:9]/YSCALE,data[:,9:16]/XSCALE,data[:,44]/.1,data[:,49:121]/ZSCALE.ravel(),data[:,126:198]/VZSCALE.ravel()]
 reference=oracle.sol(times).T;error=abs(native-np.c_[reference[:,:89],reference[:,91:163]])
 ref=ColumnReference(48);s=ref.configuration(INITIAL[2:4],INITIAL[:2]);p=ref.fields(ref.t,INITIAL[2:4],INITIAL[:2])[3]
 energy=INITIAL[:2]@s['mass']-np.sum(s['w']*s['a']*p)+.5*INITIAL[4:]@s['G']@INITIAL[4:]+.5*(INITIAL[2:4]-REST)@STIFF@(INITIAL[2:4]-REST)
 energy+=.5*np.sum(XINITIAL[:3]**2/C)+.5*np.sum(G*XINITIAL[3:]**2)+.5*H*.004**2
 graph=data[:,9:16];energies=ledger[:,0]+.5*np.sum(graph[:,:3]**2/C,axis=1)+.5*np.sum(G*graph[:,3:]**2,axis=1)+data[:,45]+data[:,121]+data[:,198]
 work=np.cumsum(ledger[:,3:5].sum(axis=1)+ledger[:,11]+data[:,42]+data[:,47]+data[:,124]+data[:,200]);loss=np.cumsum(ledger[:,5:9].sum(axis=1)+data[:,46]+data[:,123]+data[:,199])
 residual=float(np.max(abs(energies-energy+loss-work)));assert residual<1e-13,residual
 mass0=s['mass'].sum()+XINITIAL[:3].sum();mass=ledger[:,1:3].sum(axis=1)+graph[:,:3].sum(axis=1)+data[:,122]
 net=np.cumsum(dt*(ledger[:,17]+NODAL.sum()-data[:,22]));mass_error=float(np.max(abs(mass-mass0-net)));assert mass_error<1e-19,mass_error
 assert np.min(data[:,123])>=0 and np.max(abs(data[:,125]))<1e-18
 return dict(viscous_dissipation_J=float(data[:,199].sum()),independent_viscous_dissipation_J=float(oracle.y[163,-1]*1e-6),viscous_dissipation_error_J=float(abs(data[:,199].sum()-oracle.y[163,-1]*1e-6)),steps=n,endpoint_scaled_error=float(error[-1].max()),maximum_sampled_scaled_error=float(error.max()),joint_work_residual_J=residual,mass_residual_m3=mass_error,
  thermal_dissipation_J=float(data[:,123].sum()),independent_thermal_dissipation_J=float(oracle.y[89,-1]*1e-6),thermal_dissipation_error_J=float(abs(data[:,123].sum()-oracle.y[89,-1]*1e-6)))
def main():
 OUT.mkdir(parents=True,exist_ok=True);print('Independent ODE, fine reference',flush=True);fine=continuous(48,2e-10)
 print('Independent ODE, quadrature/tolerance control',flush=True);coarse=continuous(32,2e-9);times=np.linspace(0,DURATION,3201)
 delta=float(np.max(abs(fine.sol(times)-coarse.sol(times))));assert delta<2e-7,delta
 np.savez_compressed(OUT/'continuous.npz',times=times,fine=fine.sol(times),coarse=coarse.sol(times))
 runs=[]
 for n in [200,400,800,1600]:
  record=run(n,fine);runs.append(record);print(json.dumps(record),flush=True)
 ratios=[a['endpoint_scaled_error']/b['endpoint_scaled_error'] for a,b in zip(runs,runs[1:])]
 assert all(a['endpoint_scaled_error']>b['endpoint_scaled_error'] for a,b in zip(runs,runs[1:])),ratios
 paths=[Path(__file__),ROOT/'tools/slit_viscous_reference.py',ROOT/'src/slit_viscous.cpp',ROOT/'include/ducklips/slit_viscous.hpp',ROOT/'include/ducklips/viscous_spectrum.hpp',ROOT/'tools/slit_thermal_reference.py',ROOT/'tools/continuous_outlet_column_reference.py',ROOT/'tools/continuous_phonation_reference.py',ROOT/'tools/graph_continuous_glottis_reference.py',ROOT/'tools/radiation_boundary_reference.py',ROOT/'tools/graph_continuous_glottis_probe.cpp',ROOT/'src/graph_continuous_glottis.cpp',ROOT/'include/ducklips/graph_continuous_glottis.hpp',ROOT/'src/slit_thermal.cpp',ROOT/'include/ducklips/slit_thermal.hpp',ROOT/'include/ducklips/thermal_spectrum.hpp',ROOT/'build/ducklips_graph_continuous_glottis_probe',*sorted(OUT.iterdir())]
 report=dict(scope='Open continuous source with fixed cyclic tract, aperture radiation, grouped finite-slit thermal storage and three viscous segments (two sharing one flow), compared with independent time-continuous ODE. Tests trajectory accuracy of the split coupling, not speech or anatomical calibration.',reference_refinement_scaled_error=delta,duration_s=DURATION,runs=runs,endpoint_refinement_ratios=ratios,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
 (ROOT/'reports/continuous_viscous_reference.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
