#!/usr/bin/env python3
"""Sanitizer checks at motion endpoints, hidden splits, contact and worst errors."""
import hashlib
import json
from pathlib import Path
import subprocess

from planar_regions_reference import check_geometry

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'out/planar-motion-asan'
BINARY=ROOT/'.build/asan/ducklips_planar_regions_probe'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    main_receipt=ROOT/'reports/planar_motion_reference.json'
    report=json.loads(main_receipt.read_text())
    selected=[5.020,5.028625,5.030125,5.030250,5.031250,5.031375,5.062375,5.062500,5.067500,5.070250,5.085]
    results=[]
    for time in selected:
        case=min(report['samples'],key=lambda c:abs(c['time']-time))
        assert abs(case['time']-time)<1e-9
        source=ROOT/case['input'];destination=OUT/f'pose-{time:.9f}.json'
        command=[str(BINARY),'--between',format(case['fraction'],'.17g')]
        with source.open() as f,destination.open('w') as g:
            subprocess.run(command,stdin=f,stdout=g,check=True)
        native=json.loads(destination.read_text())
        p=json.loads((ROOT/case['metrics']).read_text())
        p.update(domain=native['domain'],solids=native['solids'])
        result=check_geometry(p,native)
        result.update(command=command,input=case['input'],input_sha256=sha(source),
                      output=str(destination.relative_to(ROOT)),output_sha256=sha(destination))
        results.append(result)
        print(time,'passed',flush=True)
    hashes={name:sha(ROOT/name) for name in ['include/ducklips/planar_airspace.hpp','src/planar_airspace.cpp',
            'tools/planar_regions_probe.cpp','tools/planar_regions_reference.py','tools/planar_motion_sanitizer.py',
            '.build/asan/ducklips_planar_regions_probe']}
    (OUT/'receipt.json').write_text(json.dumps(dict(source_receipt_sha256=sha(main_receipt),cases=results,hashes=hashes),indent=2)+'\n')


if __name__=='__main__':
    main()
