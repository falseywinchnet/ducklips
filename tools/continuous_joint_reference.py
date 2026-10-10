#!/usr/bin/env python3
"""Fresh joint source/tract trajectories against a hash-verified saved ODE oracle."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from continuous_viscous_reference import (ColumnReference,INITIAL,YSCALE,XSCALE,XINITIAL,C,G,NODAL,SERIES,EDGES,PORT,ZSCALE,VZSCALE,VBRANCH,VL,VGAP,GAPS,H,REST,STIFF,DURATION)
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-joint'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

class SavedOracle:
 def __init__(self,path):
  saved=np.load(path);self.times=saved['times'];self.y=saved['fine'];self.coarse=saved['coarse']
 def sol(self,times):
  indices=np.rint(times/DURATION*(len(self.times)-1)).astype(int)
  assert np.max(abs(self.times[indices]-times))<1e-18
  return self.y[:,indices]

def run(n,oracle,*,fixture='--joint-fixture',out=OUT):
 dt=DURATION/n;lines=[f'3 4 {n} {PORT}',*[f'{i} {j}' for i,j in EDGES],' '.join(map(str,INITIAL)),'1','80','1e-4 .004',' '.join(map(str,XINITIAL)),
 ' '.join(map(str,np.c_[C,GAPS].ravel())),' '.join(['0']*72),'3',*[' '.join(map(str,[VBRANCH[i],VL[i],G[VBRANCH[i]],VGAP[i]])) for i in range(3)],' '.join(['0']*72)]
 for k in range(n):lines.append(' '.join(map(str,np.r_[dt,1400.,2.,80.,1e-4,1e-4,G[-1],G[-1],np.c_[C,GAPS].ravel(),np.c_[VL,G[VBRANCH],VGAP].ravel(),C,C,G,G,NODAL,SERIES])))
 path=out/f'fixed-{n}.input';path.write_text('\n'.join(lines)+'\n')
 result=subprocess.run([str(ROOT/'build/ducklips_graph_continuous_glottis_probe'),fixture],input=path.read_text(),capture_output=True,text=True)
 path.with_suffix('.txt').write_text(result.stdout);path.with_suffix('.stderr').write_text(result.stderr);assert result.returncode==0,result.stderr
 data=np.loadtxt(path.with_suffix('.txt'));assert data.shape==(n,204),data.shape
 ledger=data[:,23:42];times=np.arange(1,n+1)*dt
 native=np.c_[data[:,:9]/YSCALE,data[:,9:16]/XSCALE,data[:,44]/.1,data[:,49:121]/ZSCALE.ravel(),data[:,126:198]/VZSCALE.ravel()]
 reference=oracle.sol(times).T;error=abs(native-np.c_[reference[:,:89],reference[:,91:163]])
 ref=ColumnReference(48);s=ref.configuration(INITIAL[2:4],INITIAL[:2]);p=ref.fields(ref.t,INITIAL[2:4],INITIAL[:2])[3]
 energy=INITIAL[:2]@s['mass']-np.sum(s['w']*s['a']*p)+.5*INITIAL[4:]@s['G']@INITIAL[4:]+.5*(INITIAL[2:4]-REST)@STIFF@(INITIAL[2:4]-REST)
 energy+=.5*np.sum(XINITIAL[:3]**2/C)+.5*np.sum(G*XINITIAL[3:]**2)+.5*H*.004**2
 graph=data[:,9:16];energies=ledger[:,0]+.5*np.sum(graph[:,:3]**2/C,axis=1)+.5*np.sum(G*graph[:,3:]**2,axis=1)+data[:,45]+data[:,121]+data[:,198]
 work=np.cumsum(ledger[:,3:5].sum(axis=1)+ledger[:,11]+data[:,42]+data[:,47]+data[:,124]+data[:,200]);loss=np.cumsum(ledger[:,5:9].sum(axis=1)+data[:,46]+data[:,123]+data[:,199]+data[:,202])
 residual=float(np.max(abs(energies-energy+loss-work)));assert residual<1e-13,residual
 mass0=s['mass'].sum()+XINITIAL[:3].sum();mass=ledger[:,1:3].sum(axis=1)+graph[:,:3].sum(axis=1)+data[:,122]
 net=np.cumsum(dt*(ledger[:,17]+NODAL.sum()-data[:,22]));mass_error=float(np.max(abs(mass-mass0-net)));assert mass_error<1e-19,mass_error
 assert min(np.min(data[:,123]),np.min(data[:,202]))>=0 and np.max(abs(data[:,203]))<1e-18
 return dict(passive_numerical_dissipation_J=float(data[:,202].sum()),maximum_passive_work_residual_J=float(np.max(abs(data[:,203]))),viscous_dissipation_J=float(data[:,199].sum()),independent_viscous_dissipation_J=float(oracle.y[163,-1]*1e-6),viscous_dissipation_error_J=float(abs(data[:,199].sum()-oracle.y[163,-1]*1e-6)),steps=n,endpoint_scaled_error=float(error[-1].max()),maximum_sampled_scaled_error=float(error.max()),joint_work_residual_J=residual,mass_residual_m3=mass_error,
  thermal_dissipation_J=float(data[:,123].sum()),independent_thermal_dissipation_J=float(oracle.y[89,-1]*1e-6),thermal_dissipation_error_J=float(abs(data[:,123].sum()-oracle.y[89,-1]*1e-6)))
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 receipt_path=ROOT/'reports/continuous_viscous_reference.json';saved_receipt=json.loads(receipt_path.read_text())
 oracle_path=ROOT/'out/continuous-viscous/continuous.npz'
 reused=['tools/continuous_viscous_reference.py','tools/continuous_outlet_column_reference.py','tools/continuous_phonation_reference.py','tools/graph_continuous_glottis_reference.py','tools/radiation_boundary_reference.py','tools/slit_thermal_reference.py','tools/slit_viscous_reference.py','out/continuous-viscous/continuous.npz']
 for name in reused:assert sha(ROOT/name)==saved_receipt['hashes'][name],name
 oracle=SavedOracle(oracle_path);delta=float(np.max(abs(oracle.y-oracle.coarse)));assert delta==saved_receipt['reference_refinement_scaled_error']
 runs=[]
 for n in [200,400,800,1600]:
  record=run(n,oracle);runs.append(record);print(json.dumps(record),flush=True)
 ratios=[a['endpoint_scaled_error']/b['endpoint_scaled_error'] for a,b in zip(runs,runs[1:])]
 assert all(a['endpoint_scaled_error']>b['endpoint_scaled_error'] for a,b in zip(runs,runs[1:])),ratios
 assert all(a['maximum_sampled_scaled_error']>b['maximum_sampled_scaled_error'] for a,b in zip(runs,runs[1:]))
 paths=[Path(__file__),receipt_path,*[ROOT/name for name in reused],ROOT/'src/relaxation_graph.cpp',ROOT/'include/ducklips/relaxation_graph.hpp',ROOT/'src/graph_joint_relaxation.cpp',ROOT/'src/graph_continuous_glottis.cpp',ROOT/'include/ducklips/graph_continuous_glottis.hpp',ROOT/'tools/graph_continuous_glottis_probe.cpp',ROOT/'build/ducklips_graph_continuous_glottis_probe',*sorted(OUT.iterdir())]
 report=dict(scope='Fresh native common passive Radau update coupled to the continuous gas/tissue source, fixed cyclic tract, aperture radiation, thermal storage and three viscous segments. Compared to a reused, hash-verified independent continuous ODE trajectory. No new ODE integration was run; all independent generator files and the saved trajectory are verified against their original receipt. Tests trajectory accuracy, not speech.',reference_refinement_scaled_error=delta,duration_s=DURATION,runs=runs,endpoint_refinement_ratios=ratios,hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})
 (ROOT/'reports/continuous_joint_reference.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
