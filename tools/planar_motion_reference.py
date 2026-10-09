#!/usr/bin/env python3
"""Audit intermediate contour geometry and quantify endpoint-metric blending.

A declared piecewise linear vertex path, not a recovered muscle trajectory.
All region/port results are independently audited at the sampled instants.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import shapely
from planar_regions_reference import check_geometry

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'out/planar-motion'
BINARY = ROOT / 'build/ducklips_planar_regions_probe'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def flatten(p,native):
    nx,ny=len(p['x_edges'])-1,len(p['y_edges'])-1
    v=np.zeros(nx*ny);a=np.zeros((nx+1)*ny);b=np.zeros(nx*(ny+1))
    for cell,volume,_,_ in native['regions']:
        v[cell]+=volume
    for _,_,face,vertical,area in native['links']:
        if vertical:a[face]+=area
        else:b[face]+=area
    return v,a,b


def pair_input(first,last,path):
    nx,ny=len(first['x_edges'])-1,len(first['y_edges'])-1
    assert first['x_edges']==last['x_edges'] and first['y_edges']==last['y_edges'] and first['depth']==last['depth']
    lines=[f'{nx} {ny} {len(first["solids"])} {first["depth"]}', ' '.join(map(str,first['x_edges'])), ' '.join(map(str,first['y_edges']))]
    for pose in (first,last):
        for poly in [pose['domain'],*pose['solids']]:
            lines.append(str(len(poly)));lines.extend(' '.join(map(str,point)) for point in poly)
    path.write_text('\n'.join(lines)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--subdivisions',type=int,default=8)
    args=parser.parse_args();division=args.subdivisions
    if division<2 or division>32:raise ValueError('bounded subdivisions required')
    OUT.mkdir(parents=True,exist_ok=True)
    paths=sorted((ROOT/'out/planar-airspace/motion').glob('pose-*.json'))
    poses=[json.loads(path.read_text()) for path in paths]
    assert len(poses)==66
    results=[];states=[];max_vertex_error=0.;worst_volume=(0.,None);worst_aperture=(0.,None)
    for step,(first,last) in enumerate(zip(poses,poses[1:])):
        source=OUT/f'interval-{step:03}.txt';pair_input(first,last,source)
        for sub in range(0 if step==0 else 1,division+1):
            fraction=sub/division;time=first['time']+fraction*(last['time']-first['time'])
            destination=OUT/f'pose-{time:.9f}.json'
            command=[str(BINARY),'--between',format(fraction,'.17g')]
            with source.open() as stream,destination.open('w') as output:
                subprocess.run(command,stdin=stream,stdout=output,check=True)
            native=json.loads(destination.read_text())
            sample=dict(time=time,depth=first['depth'],x_edges=first['x_edges'],y_edges=first['y_edges'],domain=native['domain'],solids=native['solids'])
            for p0,p1,p in zip([first['domain'],*first['solids']],[last['domain'],*last['solids']],[sample['domain'],*sample['solids']]):
                expected=(1-fraction)*np.array(p0)+fraction*np.array(p1)
                error=float(np.max(np.abs(expected-np.array(p))))
                assert error<1e-15;max_vertex_error=max(max_vertex_error,error)
            result=check_geometry(sample,native)
            v,a,b=flatten(sample,native)
            if 0<fraction<1:
                blended=(1-fraction)*np.array(first['volume'])+fraction*np.array(last['volume'])
                k=int(np.argmax(abs(v-blended)));error=float(abs(v[k]-blended[k]))
                if error>worst_volume[0]:worst_volume=(error,dict(time=time,cell=k,actual_m3=float(v[k]),blended_m3=float(blended[k])))
                for axis,actual,key in [('vertical',a,'vertical_face'),('horizontal',b,'horizontal_face')]:
                    blended=(1-fraction)*np.array(first[key])+fraction*np.array(last[key])
                    k=int(np.argmax(abs(actual-blended)));error=float(abs(actual[k]-blended[k]))
                    if error>worst_aperture[0]:worst_aperture=(error,dict(time=time,axis=axis,face=k,actual_m2=float(actual[k]),blended_m2=float(blended[k])))
            sample.update(volume=v.tolist(),vertical_face=a.tolist(),horizontal_face=b.tolist())
            # Metrics plus actual native sampled contours allow independent
            # pressure experiments without reconstructing them by interpolation.
            metrics=OUT/f'metrics-{time:.9f}.json';metrics.write_text(json.dumps(sample,separators=(',',':'))+'\n')
            result.update(command=command,input=str(source.relative_to(ROOT)),input_sha256=sha(source),
                          output=str(destination.relative_to(ROOT)),output_sha256=sha(destination),
                          metrics=str(metrics.relative_to(ROOT)),metrics_sha256=sha(metrics),interval=step,fraction=fraction)
            results.append(result);states.append(sample)
        print(json.dumps(dict(interval=step,time=last['time'],audited_samples=len(results),latest_split_cells=results[-1]['split_cells'])),flush=True)
    transitions=[]
    for old,new in zip(results,results[1:]):
        if old['components']!=new['components']:
            transitions.append(dict(after=old['time'],by=new['time'],before=old['components'],after_count=new['components']))
    report=dict(use='sampled piecewise linear contour motion; no moving pressure-state identity or speech integration',
                subdivisions=division,maximum_vertex_error_m=max_vertex_error,
                maximum_volume_blending_error_m3=worst_volume[0],worst_volume_blending=worst_volume[1],
                maximum_aperture_blending_error_m2=worst_aperture[0],worst_aperture_blending=worst_aperture[1],
                component_transition_brackets=transitions,samples=results,
                versions=dict(numpy=np.__version__,shapely=shapely.__version__,geos=shapely.geos_version_string),
                audit_command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],
                input_poses={str(path.relative_to(ROOT)):sha(path) for path in paths})
    sources=['include/ducklips/planar_airspace.hpp','src/planar_airspace.cpp','tests/planar_airspace_tests.cpp',
             'tools/planar_regions_probe.cpp','tools/planar_regions_reference.py','tools/planar_motion_reference.py','build/ducklips_planar_regions_probe']
    report['hashes']={name:sha(ROOT/name) for name in sources}
    (ROOT/'reports/planar_motion_reference.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['samples','input_poses','hashes']},indent=2),flush=True)


if __name__=='__main__':
    main()
