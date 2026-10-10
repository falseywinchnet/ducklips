#!/usr/bin/env python3
"""Verify the selectable third-order path against its preserved fixtures."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/radau5-compatibility'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(parents=True,exist_ok=True);binary=ROOT/'build/ducklips_relaxation_graph_probe';runs=[]
 old=json.loads((ROOT/'reports/relaxation_graph_reference.json').read_text())
 for source in sorted((ROOT/'out/relaxation-graph-reference').glob('*.input')):
  expected=source.with_suffix('.txt')
  for p in [source,expected]:assert sha(p)==old['hashes'][str(p.relative_to(ROOT))]
  result=subprocess.run([str(binary),'--fixture'],input=source.read_text(),capture_output=True,text=True);assert result.returncode==0,result.stderr
  target=OUT/expected.name;target.write_text(result.stdout);a=np.loadtxt(target);b=np.loadtxt(expected)
  assert a.shape==b.shape;error=float(np.max(abs(a-b)/np.maximum(1.,abs(b))));assert error<1e-15,(source,error)
  runs.append(dict(input=str(source.relative_to(ROOT)),bit_identical=result.stdout==expected.read_text(),maximum_scaled_difference=error,hashes={str(p.relative_to(ROOT)):sha(p) for p in [source,expected,target]}))
 report=dict(scope='Preserved third-order random, transient and stiff driven-flow fixtures rerun after generalizing the implementation. Whole numerical output compared at 1e-15 scaled tolerance; bit equality reported separately.',runs=runs,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),binary,ROOT/'reports/relaxation_graph_reference.json']})
 (ROOT/'reports/radau5_compatibility.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS',len(runs),'fixtures;',sum(x['bit_identical'] for x in runs),'bit identical')
if __name__=='__main__':main()
