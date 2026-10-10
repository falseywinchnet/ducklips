#!/usr/bin/env python3
"""Check coupled step-doubling against the independently integrated source ODE."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from continuous_joint_reference import SavedOracle
from continuous_viscous_reference import ROOT,ColumnReference,INITIAL,YSCALE,XSCALE,XINITIAL,C,G,H,REST,STIFF,ZSCALE,VZSCALE,NODAL,DURATION
OUT=ROOT/'out/joint-accuracy-reference'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(parents=True,exist_ok=True);original=json.loads((ROOT/'reports/continuous_viscous_reference.json').read_text());reused=[name for name in original['hashes'] if name.startswith('tools/') and name.endswith('.py')]+['out/continuous-viscous/continuous.npz']
 for name in reused:assert sha(ROOT/name)==original['hashes'][name],name
 source=ROOT/'out/continuous-joint/fixed-200.input';old=json.loads((ROOT/'reports/continuous_joint_reference.json').read_text());assert sha(source)==old['hashes'][str(source.relative_to(ROOT))]
 oracle=SavedOracle(ROOT/'out/continuous-viscous/continuous.npz');reference=oracle.sol(np.arange(1,201)*DURATION/200).T;reference=np.c_[reference[:,:89],reference[:,91:163]]
 ref=ColumnReference(48);config=ref.configuration(INITIAL[2:4],INITIAL[:2]);pressure=ref.fields(ref.t,INITIAL[2:4],INITIAL[:2])[3]
 energy=INITIAL[:2]@config['mass']-np.sum(config['w']*config['a']*pressure)+.5*INITIAL[4:]@config['G']@INITIAL[4:]+.5*(INITIAL[2:4]-REST)@STIFF@(INITIAL[2:4]-REST)+.5*np.sum(XINITIAL[:3]**2/C)+.5*np.sum(G*XINITIAL[3:]**2)+.5*H*.004**2
 mass=config['mass'].sum()+XINITIAL[:3].sum();runs=[];binary=ROOT/'build/ducklips_graph_continuous_glottis_probe'
 for tolerance in [1e-4,2.5e-5,6.25e-6]:
  command=[str(binary),'--joint-accuracy-fixture',str(tolerance),'5e-6'];result=subprocess.run(command,input=source.read_text(),capture_output=True,text=True)
  output=OUT/f'tol-{tolerance}.txt';output.write_text(result.stdout);output.with_suffix('.stderr').write_text(result.stderr);assert result.returncode==0,result.stderr
  data=np.loadtxt(output);assert data.shape==(200,175),data.shape
  native=np.c_[data[:,:9]/YSCALE,data[:,9:16]/XSCALE,data[:,16]/.1,data[:,17:89]/ZSCALE.ravel(),data[:,89:161]/VZSCALE.ravel()];error=abs(native-reference)
  work=np.cumsum(data[:,163]);loss=np.cumsum(data[:,164]+data[:,165]);work_error=float(max(abs(data[:,161]-energy+loss-work)));assert work_error<1e-13,work_error
  transfer=np.cumsum(data[:,169]-data[:,170]+DURATION/200*NODAL.sum());mass_error=float(max(abs(data[:,162]-mass-transfer)));assert mass_error<1e-19,mass_error
  assert np.min(data[:,164:166])>=0 and np.max(data[:,173])<=tolerance
  assert np.max(abs(data[:,168]))<1e-17
  record=dict(tolerance=tolerance,macro_steps=200,accepted_steps=int(data[:,171].sum()),rejected_trials=int(data[:,172].sum()),minimum_step_s=float(data[:,174].min()),maximum_accepted_local_difference=float(data[:,173].max()),maximum_sampled_scaled_error=float(error.max()),endpoint_scaled_error=float(error[-1].max()),maximum_source_effort_error_Pa=float(error[:,:2].max()*1000),maximum_tract_pressure_error_Pa=float(np.max(error[:,9:12]*XSCALE[:3]/C)),work_residual_J=work_error,mass_residual_m3=mass_error,command=command,hashes={str(p.relative_to(ROOT)):sha(p) for p in [source,output,output.with_suffix('.stderr')]})
  runs.append(record);print(json.dumps(record),flush=True)
 assert all(a['maximum_sampled_scaled_error']>b['maximum_sampled_scaled_error'] for a,b in zip(runs,runs[1:]))
 paths=[Path(__file__),ROOT/'tools/continuous_joint_reference.py',ROOT/'reports/continuous_viscous_reference.json',ROOT/'reports/continuous_joint_reference.json',*[ROOT/name for name in reused],binary,ROOT/'src/joint_accuracy.cpp',ROOT/'include/ducklips/joint_accuracy.hpp',ROOT/'tools/graph_continuous_glottis_probe.cpp']
 report=dict(scope='Fresh native coupled adaptive trajectories compared with a reused, hash-verified independent 164-state continuous ODE. Fixed tract, source, thermal/viscous/radiation memories. Local tolerance is not a global error bound. No new speech.',runs=runs,hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})
 (ROOT/'reports/joint_accuracy_reference.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
