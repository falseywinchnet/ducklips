#!/usr/bin/env python3
"""Fresh default-profile motion and constitutive checks against saved receipts."""
import concurrent.futures,hashlib,json,subprocess
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-contact'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(mode):
 old=ROOT/f'out/continuous-glottis/{mode}-100.json';previous=json.loads(old.read_text());p=subprocess.run([str(ROOT/'build/ducklips_continuous_glottis_probe'),mode,'100','.0001'],capture_output=True,text=True,check=True);(OUT/f'default-{mode}.json').write_text(p.stdout);current=json.loads(p.stdout);errors={}
 for key,scale in [('mass',1e-8),('effort',1000.),('area',1e-6),('flow',1e-5),('position',1e-4),('velocity',.1),('momentum',.001)]:
  error=float(np.max(np.abs(np.array(previous[key])-current[key]))/scale);assert error<1e-7,(mode,key,error);errors[key]=error
 return dict(mode=mode,maximum_scaled_changes=errors,reference_hash=sha(old))
def main():
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:motion=list(pool.map(run,['free','driven','lower','upper']))
 previous=json.loads((ROOT/'reports/contour_pressure_loss_validation.json').read_text());inp=ROOT/'out/slit-kinetic/loss-input.txt';expected=ROOT/'out/slit-kinetic/loss-output.txt'
 for p in [inp,expected]:assert sha(p)==previous['hashes'][str(p.relative_to(ROOT))]
 p=subprocess.run([str(ROOT/'build/ducklips_contour_pressure_loss_probe'),'--batch'],input=inp.read_text(),text=True,capture_output=True,check=True);(OUT/'default-loss.jsonl').write_text(p.stdout);old=[json.loads(s)['resistance'] for s in expected.read_text().splitlines()];new=[json.loads(s)['resistance'] for s in p.stdout.splitlines()];assert len(old)==len(new);error=max(float(np.max(np.abs(np.array(a)-b))/np.max(np.abs(a))) for a,b in zip(old,new));assert error<1e-10,error
 report=dict(scope='fresh native runs against historical default-profile fixtures, not fresh independent oracles',motion=motion,constitutive_maximum_scaled_change=error,constitutive_fixture_hashes={str(p.relative_to(ROOT)):sha(p) for p in [inp,expected]});(OUT/'regression.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS four default motion regressions and saved constitutive fixture')
if __name__=='__main__':main()
