#!/usr/bin/env python3
"""Independent time integration of moving-face openings; no speech rendering."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import scipy
from scipy.integrate import quad
import shapely
from shapely.geometry import Polygon, LineString
from shapely.ops import unary_union
from planar_motion_reference import pair_input

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometry(first,last,fraction,axis,face):
    polygons=[]
    for a,b in zip([first['domain'],*first['solids']],[last['domain'],*last['solids']]):
        polygon=Polygon((1-fraction)*np.asarray(a)+fraction*np.asarray(b))
        assert polygon.is_valid and polygon.area>0
        polygons.append(polygon)
    nx=len(first['x_edges'])-1
    if axis=='v':
        row,col=divmod(face,nx+1);x=first['x_edges'][col]
        line=LineString([(x,first['y_edges'][row]),(x,first['y_edges'][row+1])])
    else:
        row,col=divmod(face,nx);y=first['y_edges'][row]
        line=LineString([(first['x_edges'][col],y),(first['x_edges'][col+1],y)])
    # Selected audit faces are internal to the domain, not coincident outer
    # boundary faces. Solid boundaries are removed by difference, including
    # a solid edge coincident with the line.
    return line.intersection(polygons[0]).difference(unary_union(polygons[1:])).length*first['depth']


def closed_spans(windows):
    spans=[]
    for a,b,area,_ in windows:
        if area!=0:continue
        if spans and spans[-1][1]==a:spans[-1][1]=b
        else:spans.append([a,b])
    return spans


def audit(name,first,last,axis,face,binary,out):
    source=out/f'{name}.txt';pair_input(first,last,source)
    destination=out/f'{name}.json';command=[str(binary),f'--face-{axis}',str(face)]
    with source.open() as f,destination.open('w') as g:subprocess.run(command,stdin=f,stdout=g,check=True)
    native=json.loads(destination.read_text());windows=native['windows']
    assert windows[0][0]==0 and windows[-1][1]==1
    assert all(a<b and v>=0 and e>=0 for a,b,v,e in windows)
    assert all(a[1]==b[0] for a,b in zip(windows,windows[1:]))
    f=lambda t:geometry(first,last,t,axis,face)
    independent=[]
    # These quadratures receive uniform partitions, never native event roots.
    # They are an independent convergence check, not a certified error bound.
    for divisions in [32,128,512,2048]:
        result=quad(f,0,1,points=np.linspace(0,1,divisions+1),epsabs=1e-12,epsrel=1e-12,
                    limit=2*divisions+1024,full_output=True)
        value,error=result[:2];warning=result[3] if len(result)>3 else None
        independent.append(dict(divisions=divisions,average_area_m2=value,estimated_error_m2=error,
                                evaluations=result[2]['neval'],warning=warning))
        if divisions>=128 and warning is None:break
    assert independent[-1]['warning'] is None,(name,independent)
    error=abs(native['average_area']-independent[-1]['average_area_m2'])
    refinement=abs(independent[-2]['average_area_m2']-independent[-1]['average_area_m2'])
    assert error<3e-12 and refinement<3e-12,(name,error,refinement)
    spans=closed_spans(windows);closed_max=0
    for a,b in spans:
        # Independently inspect each closed run at three interior times.
        for fraction in (.211324865405187,.5,.788675134594813):
            closed_max=max(closed_max,f(a+(b-a)*fraction))
    assert closed_max<1e-16,(name,closed_max)
    result=dict(name=name,axis=axis,face=face,time_first=first['time'],time_last=last['time'],
                average_area_m2=native['average_area'],native_quadrature_error_estimate_m2=native['error_estimate'],
                independent=independent,absolute_error_m2=error,independent_refinement_difference_m2=refinement,
                event_windows=len(windows),closed_fraction_intervals=spans,
                closed_duration_fraction=sum(b-a for a,b in spans),maximum_closed_probe_area_m2=closed_max,
                endpoint_trapezoid_area_m2=.5*(f(0)+f(1)),command=command,
                input=str(source.relative_to(ROOT)),input_sha256=sha(source),output=str(destination.relative_to(ROOT)),output_sha256=sha(destination))
    print(json.dumps({k:result[k] for k in ['name','average_area_m2','absolute_error_m2','event_windows','closed_duration_fraction']}),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_planar_regions_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/planar-face-motion')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_face_motion_reference.json')
    parser.add_argument('--sanitizer',action='store_true',help='small selected subset, full independent integral checks')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    paths=sorted((ROOT/'out/planar-airspace/motion').glob('pose-*.json'))
    poses=[json.loads(p.read_text()) for p in paths];assert len(poses)==66
    selected=[(50,'v',3743)]
    for interval in ([11,42] if args.sanitizer else [0,8,10,11,12,29,40,42,47,51,55,64]):
        for axis,key in [('v','vertical_face'),('h','horizontal_face')]:
            face=int(np.argmax(abs(np.asarray(poses[interval+1][key])-np.asarray(poses[interval][key]))))
            selected.append((interval,axis,face))
    results=[]
    for interval,axis,face in selected:
        results.append(audit(f'foxtrot-{interval:03}-{axis}-{face}',poses[interval],poses[interval+1],axis,face,args.binary,args.output))
    # Convex triangles remain convex under the positive diagonal affine paths;
    # two independently moving solids exercise changes in union ownership.
    rng=np.random.default_rng(91374)
    for index in range(2 if args.sanitizer else 24):
        base=np.array([[-1.2,-.9],[1.3,-.6],[.1,1.4]])
        p=dict(time=0.,depth=.022,x_edges=[-3.,0.,3.],y_edges=[-2.,2.],domain=[[-3,-2],[3,-2],[3,2],[-3,2]],solids=[])
        q=dict(p,time=1.,solids=[])
        for solid in range(2):
            p['solids'].append((base*rng.uniform(.4,1.2,2)+rng.uniform(-.7,.7,2)).tolist())
            q['solids'].append((base*rng.uniform(.4,1.2,2)+rng.uniform(-.7,.7,2)).tolist())
        results.append(audit(f'affine-triangles-{index:02}',p,q,'v',1,args.binary,args.output))
    sources=['include/ducklips/planar_face_motion.hpp','src/planar_face_motion.cpp','tests/planar_face_motion_tests.cpp',
             'tools/planar_regions_probe.cpp','tools/planar_face_motion_reference.py','tools/planar_motion_reference.py']
    report=dict(use='geometric temporal face openings; no moving region identity or production speech integration',
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],cases=results,
                maximum_absolute_integral_error_m2=max(r['absolute_error_m2'] for r in results),
                maximum_independent_refinement_difference_m2=max(r['independent_refinement_difference_m2'] for r in results),
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),args.binary]},
                input_poses={str(p.relative_to(ROOT)):sha(p) for p in paths},
                versions=dict(numpy=np.__version__,scipy=scipy.__version__,shapely=shapely.__version__,geos=shapely.geos_version_string))
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',len(results),'cases; max integral error',report['maximum_absolute_integral_error_m2'],flush=True)


if __name__=='__main__':main()
