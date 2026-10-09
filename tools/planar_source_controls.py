#!/usr/bin/env python3
"""Closed-network controls and exact historical output regressions for the source adapter."""
import sys,json,csv,copy,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
from planar_source_reference import setup,run,audit,ROOT,sha
from lip_profile_reference import vector
out=ROOT/'out/planar-source-controls';out.mkdir(exist_ok=True)
rows=list(csv.DictReader((ROOT/'out/lip-profile/nato.csv').open()));airway=list(csv.DictReader((ROOT/'out/nato-spatial-airway.csv').open()))
p,_q,_=setup(1.285);q=copy.deepcopy(p);q['time']+=.001;profiles=[vector(rows[1285])]*2
epi=np.array([[[float(airway[1285][f'{field}_{i}']) for i in reversed(range(8))] for field in ['area','length']]]*2)
receipts=[]
for name,area,reservoir in [('equilibrium',3e-6,700.),('sealed-source',0.,900.)]:
 native,receipt=run(name,p,q,profiles,epi,16,700.,area,reservoir,None,ROOT/'build/ducklips_planar_inertive_motion_probe',out,False)
 checked,_=audit(native,p,q,profiles,epi,700.,area,reservoir,None,True);receipt.update(checked)
 receipt['maximum_pressure_drift_Pa']=max(abs(v-700) for s in native['steps'] for v in s['pressure'])
 receipt['maximum_source_flow_m3_s']=max(abs(s['source_flow']) for s in native['steps'])
 assert receipt['maximum_pressure_drift_Pa']<1e-7
 assert receipt['maximum_source_flow_m3_s']<1e-14
 receipts.append(receipt)
# Historical full JSON equivalence verifies changed indexing did not redirect lip outlet.
regressions=[]
for path,args in [(ROOT/'out/planar-lip-profile/audit/closing-0-16.txt',['--lip-profile']),
                  (ROOT/'out/planar-lip-profile/audit/opening-0-16.txt',['--lip-profile']),
                  (ROOT/'out/articulated-lip/closure-128.txt',['--lip-tube']),
                  (ROOT/'out/articulated-lip/release-128.txt',['--lip-tube']),
                  (ROOT/'out/articulated-lip/lip-narrowing-128.txt',['--lip-tube']),
                  (ROOT/'out/articulated-lip/original-closure-128.txt',[]),
                  (ROOT/'out/articulated-lip/original-release-128.txt',[])]:
 with path.open() as f:r=subprocess.run([str(ROOT/'build/ducklips_planar_inertive_motion_probe'),*args],stdin=f,capture_output=True,text=True)
 assert r.returncode==0,(path,r.stderr)
 assert json.loads(r.stdout)==json.loads(path.with_suffix('.json').read_text()),path
 regressions.append(dict(input=str(path.relative_to(ROOT)),input_sha256=sha(path),historical_output_sha256=sha(path.with_suffix('.json')),exact=True))
(out/'receipt.json').write_text(json.dumps(dict(controls=receipts,regressions=regressions,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__).resolve(),ROOT/'tools/planar_source_reference.py',ROOT/'build/ducklips_planar_inertive_motion_probe']}),indent=2)+'\n')
print(json.dumps(dict(controls=[dict(name=r['name'],pressure_drift=r['maximum_pressure_drift_Pa'],flow=r['maximum_source_flow_m3_s']) for r in receipts],regressions=len(regressions))))
