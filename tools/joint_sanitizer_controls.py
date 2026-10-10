#!/usr/bin/env python3
"""Replay the independently audited moving inputs through AddressSanitizer."""
import argparse,hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/planar-joint-asan'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def difference(a,b,path=""):
 if isinstance(a,dict):
  assert a.keys()==b.keys();return max([difference(a[k],b[k],path+"."+k) for k in a]+[0.])
 if isinstance(a,list):
  assert len(a)==len(b);return max([difference(x,y,path) for x,y in zip(a,b)]+[0.])
 if isinstance(a,(float,int)):
  tolerance=1e-6 if path.endswith(".effort") else 2e-9*max(1.,abs(a),abs(b))
  assert abs(a-b)<=tolerance,(path,a,b)
  return abs(a-b)/max(1.,abs(a),abs(b))
 assert a==b;return 0.
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,default=ROOT/'out/planar-joint');parser.add_argument('--output',type=Path,default=OUT);parser.add_argument('--subdivisions',type=int,default=256);parser.add_argument('--fifth',action='store_true');parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_joint_asan_reference.json');args=parser.parse_args();out=args.output.resolve()
 out.mkdir(parents=True,exist_ok=True);binary=ROOT/'.build/asan/ducklips_planar_inertive_motion_probe';runs=[]
 for name in ['closing','opening']:
  for interval in [0,1]:
   path=args.input.resolve()/f'{name}-{interval}-{args.subdivisions}.txt';reference=path.with_suffix('.json')
   assert path.exists() and reference.exists(),path
   command=[str(binary),'--lip-profile','--epilarynx','--continuous-folds',str(ROOT/'bodies/reference.body'),'--compact-geometry','--radiation','--thermal-slit-width','.03','--viscous-slit-width','.03','--lip-viscous-width','.03','--outlet-viscous','--joint-relaxation']
   if interval:command+=['--initial-state']
   if args.fifth:command+=['--joint-fifth']
   result=subprocess.run(command,input=path.read_text(),capture_output=True,text=True)
   target=out/path.with_suffix('.json').name;target.write_text(result.stdout);error=target.with_suffix('.stderr');error.write_text(result.stderr)
   assert result.returncode==0,(result.returncode,result.stderr)
   data=json.loads(result.stdout)
   for step in data['steps']:
    assert abs(step['joint_passive']['energy_residual'])<1e-15 and abs(step['joint_passive']['air_residual'])<1e-18
    assert abs(step['source_work_residual'])<1e-12 and abs(step['source_port_residual'])<1e-17
   assert result.stderr==f"steps={len(data['steps'])}; inviscid canonical moving-volume reference; no speech\n",result.stderr
   maximum=difference(data,json.loads(reference.read_text()))
   runs.append(dict(case=name,interval=interval,steps=len(data['steps']),maximum_relative_difference=maximum,command=command,hashes={str(p.relative_to(ROOT)):sha(p) for p in [path,reference,target,error]}))
 report=dict(passive_order=5 if args.fifth else 3,scope='Fresh AddressSanitizer runs using the release audit inputs. Full JSON compared against Release (2e-9 scaled tolerance, source effort 1e-6 absolute because its reconstruction cancels large terms); per-step passive/source work and air receipts checked; independent full-stage audit is in the corresponding moving reference receipt and is not rerun here.',runs=runs,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),binary,ROOT/'bodies/reference.body']})
 args.report.resolve().write_text(json.dumps(report,indent=2)+'\n');print('PASS',len(runs),'sanitized moving intervals')
if __name__=='__main__':main()
