#!/usr/bin/env python3
"""Actual geometry with coupled local-error control and complete macro receipts."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from planar_continuous_reference import ROOT,setup,vector,run,sha

def energy(row):return row['continuous_source']['energy']+.5*np.sum(np.array(row['charge'])**2/row['compliance'])+.5*np.sum(np.array(row['conductance'])*np.array(row['momentum'])**2)+row['radiation_state']['energy']+sum(s['energy'] for kind in ['thermal_sites','viscous_sites'] for s in row[kind])
def air(row):return sum(row['volume'])+sum(row['charge'])+sum(row['continuous_source']['mass'])+sum(s['stored_air'] for s in row['thermal_sites'])
def audit(data,tolerance):
 initial=data['initial'];old=initial;work=loss=net=0.;residual=mass_residual=0.;minimum=1.;accepted=rejected=0;times=[];clock=0.
 for row in data['steps']:
  dt=row['dt'];r=row['accuracy'];assert r['maximum_accepted_error']<=tolerance
  assert r['physical_loss']>=0 and r['numerical_loss']>=0
  assert abs(r['source_port_work']-r['graph_port_work'])<1e-17
  assert abs(energy(row)-energy(old)-r['energy_change'])<1e-13
  assert len(r['accepted_end_times'])==r['accepted_steps'] and all(t>clock for t in r['accepted_end_times'])
  assert all(b>a for a,b in zip(r['accepted_end_times'],r['accepted_end_times'][1:]))
  clock+=dt;assert abs(r['accepted_end_times'][-1]-clock)<1e-17
  work+=r['external_work'];loss+=r['physical_loss']+r['numerical_loss'];net+=r['reservoir_volume']-dt*r['average_flow'][6]
  residual=max(residual,abs(energy(row)-energy(initial)+loss-work));mass_residual=max(mass_residual,abs(air(row)-air(initial)-net));minimum=min(minimum,r['minimum_step']);accepted+=r['accepted_steps'];rejected+=r['rejected'];times+=r['accepted_end_times'];old=row
 assert residual<1e-12 and mass_residual<1e-17,(residual,mass_residual)
 return dict(accepted_steps=accepted,rejected_trials=rejected,minimum_step_s=minimum,work_residual_J=residual,mass_residual_m3=mass_residual,initial_air_m3=air(initial),final_air_m3=air(old),net_external_air_m3=net)
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--tolerances',nargs='+',type=float,default=[1e-4,2.5e-5,6.25e-6]);parser.add_argument('--subdivisions',type=int,default=64);parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_planar_inertive_motion_probe');parser.add_argument('--output',type=Path,default=ROOT/'out/joint-accuracy-moving');parser.add_argument('--report',type=Path,default=ROOT/'reports/joint_accuracy_moving_reference.json');args=parser.parse_args()
 out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
 source_paths=[Path(__file__),args.binary.resolve(),ROOT/'tools/planar_continuous_reference.py',ROOT/'tools/planar_inertive_motion_probe.cpp',ROOT/'src/joint_accuracy.cpp',ROOT/'include/ducklips/joint_accuracy.hpp',ROOT/'src/continuous_glottis.cpp',ROOT/'src/graph_joint_relaxation.cpp',ROOT/'src/relaxation_graph.cpp',ROOT/'include/ducklips/graph_continuous_glottis.hpp',ROOT/'bodies/reference.body']
 starting_hashes={str(p.relative_to(ROOT)):sha(p) for p in source_paths}
 rows=list(csv.DictReader((ROOT/'out/lip-profile/nato.csv').open()));airway=list(csv.DictReader((ROOT/'out/nato-spatial-airway.csv').open()));runs=[];comparisons=[];inputs={}
 for name,start,pressure in [('closing',1284,0.),('opening',1374,700.)]:
  prior=None
  for tolerance in args.tolerances:
   folder=out/f'tol-{tolerance}';folder.mkdir(exist_ok=True);previous=None;history=[];first_air=None;net=0.
   for interval in [0,1]:
    index=start+interval;p,q,paths=setup(index/1000);profiles=np.array([vector(rows[i]) for i in [index,index+1]])
    epi=np.array([[[float(airway[t][f'{field}_{i}']) for i in reversed(range(8))] for field in ['area','length']] for t in [index,index+1]])
    data,receipt=run(f'{name}-{interval}',p,q,profiles,epi,args.subdivisions,pressure,previous,args.binary.resolve(),folder,True,.03,.03,.03,8,True,True,True,tolerance,1e-5)
    receipt.update(audit(data,tolerance));receipt['tolerance']=tolerance;runs.append(receipt);history+=data['steps'];previous=data['steps'][-1]
    if first_air is None:first_air=receipt['initial_air_m3']
    net+=receipt['net_external_air_m3'];inputs.update({str(path.relative_to(ROOT)):sha(path) for path in paths})
    print(json.dumps(dict(case=name,interval=interval,tolerance=tolerance,accepted_steps=receipt['accepted_steps'],minimum_step_s=receipt['minimum_step_s'])),flush=True)
   assert abs(air(previous)-first_air-net)<1e-17
   now=np.array([r['pressure'] for r in history]);record=dict(case=name,tolerance=tolerance,endpoint_pressure=now[-1].tolist())
   if prior is not None:
    delta=abs(now-prior);record.update(maximum_sampled_pressure_difference_Pa=float(delta.max()),endpoint_pressure_difference_Pa=float(delta[-1].max()),maximum_throat_difference_Pa=float(delta[:,6:].max()),maximum_lip_difference_Pa=float(delta[:,2:6].max()))
   comparisons.append(record);prior=now;print(json.dumps(record),flush=True)
 assert all(sha(ROOT/name)==value for name,value in starting_hashes.items()),'source or binary changed during trajectory run'
 report=dict(scope='Local step-doubling of the complete source/passive system; geometry recomputed at actual trial endpoints and midpoint. Direct whole-request energy and air checks. Adjacent tolerances are not errors against an independent whole-moving-source oracle; fixed-source ODE accuracy is checked separately. No source contact transition rule or new speech.',runs=runs,comparisons=comparisons,input_poses=inputs,hashes=starting_hashes)
 args.report.resolve().write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
