#!/usr/bin/env python3
"""Independent domain-cap area, region, centroid and normal audit."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import shapely
from shapely.geometry import LineString,Polygon
from shapely.ops import unary_union
from planar_regions_reference import check_geometry
from planar_region_motion_reference import exact_regions

ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def check_boundaries(p,native,selected,regions=None):
    if regions is None:regions,_,_=exact_regions(p,native)
    tissue=unary_union([Polygon(c) for c in p['solids']])
    points=np.array(p['domain']);relative=points-points[0]
    signed=sum(relative[i,0]*relative[i+1,1]-relative[i,1]*relative[i+1,0] for i in range(1,len(points)-1))
    expected={}
    for edge in selected:
        a,b=points[edge],points[(edge+1)%len(points)];u=b-a
        normal=np.sign(signed)*np.array([u[1],-u[0]])/np.linalg.norm(u)
        aperture=LineString([a,b]).difference(tissue)
        for i,region in enumerate(regions):
            opening=aperture.intersection(region)
            if opening.length>0:
                centre=opening.centroid;expected[(edge,i)]=(opening.length*p['depth'],np.array([centre.x,centre.y]),normal)
    actual={(int(e),int(i)):(area,np.array([x,y]),np.array([nx,ny])) for e,i,area,x,y,nx,ny in native['boundary_ports']}
    assert len(actual)==len(native['boundary_ports']) and actual.keys()==expected.keys(),(p['time'],set(actual)-set(expected),set(expected)-set(actual))
    errors=dict(area_m2=0.,centre_m=0.,normal=0.)
    for key,(area,centre,normal) in expected.items():
        a,c,n=actual[key];errors['area_m2']=max(errors['area_m2'],abs(a-area))
        errors['centre_m']=max(errors['centre_m'],float(max(abs(c-centre))))
        errors['normal']=max(errors['normal'],float(max(abs(n-normal))))
    assert errors['area_m2']<1e-15 and errors['centre_m']<1e-12 and errors['normal']<1e-12,errors
    return expected,errors


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_planar_regions_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/planar-boundary')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_boundary_reference.json')
    parser.add_argument('--sanitizer',action='store_true')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True);runs=[]
    paths=sorted((ROOT/'out/planar-boundary/poses').glob('pose-*.json'))
    assert len(paths)==66
    if args.sanitizer:paths=[paths[i] for i in [0,11,12,20,42,43,65]]
    for path in paths:
        p=json.loads(path.read_text());old_path=ROOT/'out/planar-airspace/motion'/path.name;old=json.loads(old_path.read_text())
        assert {k:v for k,v in p.items() if k!='boundary_edges'}==old,'adding markers changed historical geometry'
        selected=list(p['boundary_edges'].values());source=args.output/(path.stem+'.txt');destination=source.with_suffix('.ports.json')
        lines=[f'{len(p["x_edges"])-1} {len(p["y_edges"])-1} {len(p["solids"])} {p["depth"]}',
               ' '.join(map(str,p['x_edges'])),' '.join(map(str,p['y_edges']))]
        for poly in [p['domain'],*p['solids']]:lines.append(str(len(poly)));lines.extend(' '.join(map(str,q)) for q in poly)
        lines.extend([str(len(selected)),' '.join(map(str,selected))]);source.write_text('\n'.join(lines)+'\n')
        with source.open() as f,destination.open('w') as g:subprocess.run([str(args.binary),'--ports'],stdin=f,stdout=g,check=True)
        native=json.loads(destination.read_text());interior=check_geometry(p,native);expected,errors=check_boundaries(p,native,selected)
        totals={name:sum(a for (e,i),(a,c,n) in expected.items() if e==edge) for name,edge in p['boundary_edges'].items()}
        runs.append(dict(time=p['time'],ports=len(expected),boundary_edges=p['boundary_edges'],areas_m2=totals,errors=errors,interior=interior,
                         **{k:str(v.relative_to(ROOT)) for k,v in [('pose',path),('previous_pose',old_path),('input',source),('output',destination)]},
                         **{k+'_sha256':sha(v) for k,v in [('pose',path),('previous_pose',old_path),('input',source),('output',destination)]}))
        print(json.dumps(dict(time=p['time'],ports=len(expected),areas_m2=totals,errors=errors)),flush=True)
    sources=['tools/planar_boundary_reference.py','tools/planar_regions_reference.py','tools/planar_region_motion_reference.py',
             'tools/planar_regions_probe.cpp','tools/airspace_geometry_probe.cpp','include/ducklips/planar_airspace.hpp','src/planar_airspace.cpp']
    report=dict(use='end-cap attachments of the tongue/floor domain; not complete lip or glottal anatomy',runs=runs,
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),args.binary]},
                versions=dict(numpy=np.__version__,shapely=shapely.__version__,geos=shapely.geos_version_string))
    args.report.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
