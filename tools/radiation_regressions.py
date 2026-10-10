#!/usr/bin/env python3
"""Exact legacy and non-radiating continuous-path output regressions."""
import argparse,hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'out/radiation-regressions');parser.add_argument('--include-radiation',action='store_true');options=parser.parse_args()
OUT=options.output.resolve();OUT.mkdir(parents=True,exist_ok=True)
cases=[('out/planar-lip-profile/audit/closing-0-16.txt',['--lip-profile']),('out/planar-lip-profile/audit/opening-0-16.txt',['--lip-profile']),
 ('out/articulated-lip/closure-128.txt',['--lip-tube']),('out/articulated-lip/release-128.txt',['--lip-tube']),('out/articulated-lip/lip-narrowing-128.txt',['--lip-tube']),
 ('out/articulated-lip/original-closure-128.txt',[]),('out/articulated-lip/original-release-128.txt',[])]
for name in ['closing','opening']:
 for interval in [0,1]:
  args=['--lip-profile','--epilarynx','--continuous-folds',str(ROOT/'bodies/reference.body'),'--compact-geometry']
  if interval:args+=['--initial-state']
  cases.append((f'out/planar-continuous/{name}-{interval}-256.txt',args))
if options.include_radiation:
 for name in ['closing','opening']:
  for interval in [0,1]:
   command=['--lip-profile','--epilarynx','--continuous-folds',str(ROOT/'bodies/reference.body'),'--compact-geometry','--radiation']
   if interval:command+=['--initial-state']
   cases.append((f'out/planar-radiation/{name}-{interval}-256.txt',command))
records=[];binary=ROOT/'build/ducklips_planar_inertive_motion_probe'
for i,(path,args) in enumerate(cases):
 p=ROOT/path
 result=subprocess.run([str(binary),*args],input=p.read_text(),capture_output=True,text=True)
 assert result.returncode==0,(path,result.stderr)
 assert json.loads(result.stdout)==json.loads(p.with_suffix('.json').read_text()),path
 target=OUT/f'{i}.json';target.write_text(result.stdout)
 records.append(dict(input=path,arguments=args,exact=True,hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [p,p.with_suffix('.json'),target]}))
report=dict(cases=records,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),binary]})
(OUT/'receipt.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS',len(records),'exact historical JSON regressions')
