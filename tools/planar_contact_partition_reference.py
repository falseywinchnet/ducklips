#!/usr/bin/env python3
"""Locate two Foxtrot contacts and retain their neck in a pressure partition.

Selected-event geometry and frozen-pressure controls, not a full contact
scheduler or a moving speech trajectory. The linear mobility is a declared
numerical fixture, not calibrated vocal-tract friction.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import shapely
from shapely.geometry import Polygon,box
from shapely.ops import unary_union
from planar_regions_reference import check_geometry
import cavity_graph_reference as pressure_reference

ROOT=Path(__file__).resolve().parents[1]
BULK=1.14*350**2
MOBILITY=1e-7
DT=.0001


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def cross(a,b):return a[0]*b[1]-a[1]*b[0]


def bounds_for(p,cell):
    nx=len(p['x_edges'])-1;row,col=divmod(cell,nx)
    return np.array([p['x_edges'][col],p['y_edges'][row],p['x_edges'][col+1],p['y_edges'][row+1]])


def candidates(p,q,cell):
    """Vertex-edge contact candidates in this cell, not all arrangement events."""
    bounds=bounds_for(p,cell);grid=np.array(box(*bounds).exterior.coords)[:-1]
    first=[np.array(c) for c in [p['domain'],*p['solids']]]+[grid]
    last=[np.array(c) for c in [q['domain'],*q['solids']]]+[grid]
    result=[]
    for ci,(poly,end) in enumerate(zip(first,last)):
        for vi,(v,z) in enumerate(zip(poly,end)):
            if np.any(np.maximum(v,z)<bounds[:2]) or np.any(np.minimum(v,z)>bounds[2:]):continue
            dv=z-v
            for cj,(other,new) in enumerate(zip(first,last)):
                if ci==cj:continue
                for ei in range(len(other)):
                    a,b=other[ei],other[(ei+1)%len(other)];aa,bb=new[ei],new[(ei+1)%len(other)]
                    hull=np.array([a,b,aa,bb])
                    if np.any(hull.max(axis=0)<bounds[:2]) or np.any(hull.min(axis=0)>bounds[2:]):continue
                    da,db=aa-a,bb-b;w,dw=v-a,dv-da;e,de=b-a,db-da
                    polynomial=np.array([cross(w,e),cross(w,de)+cross(dw,e),cross(dw,de)])
                    if not np.any(polynomial):continue
                    polynomial/=max(abs(polynomial))
                    for root in np.polynomial.polynomial.polyroots(polynomial):
                        if abs(root.imag)>1e-9 or not 0<root.real<1:continue
                        t=float(root.real);point=v+t*dv;edge_a=a+t*da;direction=b+t*db-edge_a
                        if np.dot(direction,direction)==0:continue
                        position=np.dot(point-edge_a,direction)/np.dot(direction,direction)
                        if not -1e-10<=position<=1+1e-10 or np.any(point<bounds[:2]-1e-12) or np.any(point>bounds[2:]+1e-12):continue
                        result.append(dict(fraction=t,point_m=point.tolist(),vertex=[ci,vi],edge=[cj,ei],normalized_polynomial=polynomial.tolist()))
    return sorted(result,key=lambda r:r['fraction'])


def contours(p,q,t):
    return [((1-t)*np.array(a)+t*np.array(b)).tolist() for a,b in zip([p['domain'],*p['solids']],[q['domain'],*q['solids']])]


def fluid(p,q,t):
    polygons=[Polygon(c) for c in contours(p,q,t)];assert all(g.is_valid for g in polygons)
    return polygons[0].difference(unary_union(polygons[1:]))


def components(g):return len([part for part in shapely.get_parts(g) if part.area>0])


def locate(p,q,cell,bracket,before,after):
    all_candidates=candidates(p,q,cell)
    selected=[]
    for c in all_candidates:
        t=c['fraction']
        if bracket[0]<t<bracket[1] and components(fluid(p,q,t-1e-7))==before and components(fluid(p,q,t+1e-7))==after:
            selected.append(c)
    assert len(selected)==1,selected
    event=selected[0];lo,hi=bracket
    assert components(fluid(p,q,lo))==before and components(fluid(p,q,hi))==after
    for _ in range(60):
        mid=(lo+hi)/2
        if mid==lo or mid==hi:break
        count=components(fluid(p,q,mid));assert count in [before,after]
        if count==before:lo=mid
        else:hi=mid
    error=abs(event['fraction']-(lo+hi)/2);assert error<1e-10
    event.update(cell=cell,time=p['time']+event['fraction']*(q['time']-p['time']),
                 component_counts=[before,after],independent_bisection_fraction_bracket=[lo,hi],root_difference_fraction=error,
                 candidate_count=len(all_candidates))
    return event


def local_pose(p,q,t,cell,cut):
    bounds=bounds_for(p,cell);c=contours(p,q,t)
    domain=Polygon(c[0]).intersection(box(*bounds));assert domain.geom_type=='Polygon'
    return dict(time=p['time']+t*(q['time']-p['time']),depth=p['depth'],domain=list(domain.exterior.coords)[:-1],solids=c[1:],
                x_edges=[bounds[0],cut,bounds[2]] if cut is not None else [bounds[0],bounds[2]],y_edges=[bounds[1],bounds[3]])


def native_geometry(name,p,binary,out):
    lines=[f'{len(p["x_edges"])-1} 1 {len(p["solids"])} {p["depth"]}',
           ' '.join(map(str,p['x_edges'])),' '.join(map(str,p['y_edges']))]
    for c in [p['domain'],*p['solids']]:
        lines.append(str(len(c)));lines.extend(' '.join(map(str,x)) for x in c)
    source=out/f'{name}.txt';destination=source.with_suffix('.json');source.write_text('\n'.join(lines)+'\n')
    with source.open() as f,destination.open('w') as g:subprocess.run([str(binary)],stdin=f,stdout=g,check=True)
    native=json.loads(destination.read_text());audit=check_geometry(p,native)
    audit.update(command=[str(binary)],input=str(source.relative_to(ROOT)),input_sha256=sha(source),output=str(destination.relative_to(ROOT)),output_sha256=sha(destination))
    return native,audit


def partition_control(name,coarse,fine,binary,out):
    assert len(coarse['regions'])==1 and len(fine['regions'])==2
    parent=coarse['regions'][0][1];v=np.array([r[1] for r in fine['regions']]);pressure=600.
    lines=['1 1',f'{parent:.17g} {parent/BULK:.17g} {pressure}',f'2 0 {DT}']
    lines.extend(f'0 {x:.17g} {x/BULK:.17g} {x:.17g} {x/BULK:.17g} 0' for x in v)
    source=out/f'{name}-partition.txt';destination=source.with_suffix('.csv');source.write_text('\n'.join(lines)+'\n')
    command=[str(binary),'--partition']
    with source.open() as f,destination.open('w') as g:
        run=subprocess.run(command,stdin=f,stdout=g,stderr=subprocess.PIPE,text=True,check=True)
    data=np.loadtxt(destination,delimiter=',',skiprows=1);actual=data[data[:,0]==1,2];assert len(actual)==2
    mass_error=abs(sum(v*(1+actual/BULK))-parent*(1+pressure/BULK))
    energy_error=abs(.5*np.dot(v/BULK,actual**2)-.5*parent/BULK*pressure**2)
    assert max(abs(actual-pressure))<1e-8 and mass_error<1e-21 and energy_error<1e-20
    return dict(maximum_pressure_change_Pa=float(max(abs(actual-pressure))),mass_error_m3=float(mass_error),energy_error_J=float(energy_error),
                command=command,input=str(source.relative_to(ROOT)),input_sha256=sha(source),output=str(destination.relative_to(ROOT)),output_sha256=sha(destination),native_report=run.stderr.strip())


def pressure_control(name,p,native):
    assert len(native['regions'])==2 and len(native['links'])==1
    assert [r[0] for r in native['regions']]==[0,1]
    a=native['links'][0][4];distance=.5*(p['x_edges'][-1]-p['x_edges'][0]);g=MOBILITY*a/distance
    v=np.array([r[1] for r in native['regions']]);c=v/BULK;initial=np.array([600.,0.])
    links=[(0,1,a,1/g,0,0,0)]
    pressure,flow,receipt=pressure_reference.run(name,v,c,initial,np.zeros(2),links,DT)
    expected_flow=g*600/(1+DT*g*(1/c[0]+1/c[1]))
    expected=np.array([600-DT*expected_flow/c[0],DT*expected_flow/c[1]])
    error=float(max(abs(pressure-expected)));mass=float(abs(np.dot(c,pressure-initial)))
    local_mass_error=float(max(abs(c*(pressure-initial)+DT*flow[0]*np.array([1.,-1.]))))
    assert error<1e-7 and mass<1e-19 and local_mass_error<1e-21
    mixed=float(np.dot(c,initial)/sum(c))
    return dict(opening_m2=a,conductance_m3_per_Pa_s=g,pressure_Pa=pressure.tolist(),flow_m3_per_s=float(flow[0]),
                analytic_pressure_error_Pa=error,mass_error_m3=mass,local_mass_error_m3=local_mass_error,instantaneously_merged_pressure_Pa=mixed,**receipt)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry-binary',type=Path,default=ROOT/'build/ducklips_planar_regions_probe')
    parser.add_argument('--pressure-binary',type=Path,default=ROOT/'build/ducklips_cavity_graph_probe')
    parser.add_argument('--partition-binary',type=Path,default=ROOT/'build/ducklips_cavity_region_motion_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/contact-partition')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/contact_partition_reference.json')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    pressure_reference.OUT=args.output;pressure_reference.BINARY=args.pressure_binary
    results=[];input_paths=[]
    for label,start,cell,bracket,before,after,direction in [
        ('closure',5.031,3760,[.25,.375],1,2,-1),('release',5.062,3758,[.375,.5],2,1,1)]:
        paths=[ROOT/f'out/planar-airspace/motion/pose-{time:.3f}.json' for time in [start,start+.001]]
        input_paths.extend(paths);p,q=[json.loads(path.read_text()) for path in paths]
        event=locate(p,q,cell,bracket,before,after);cut=event['point_m'][0];samples=[]
        for exponent in range(3,9):
            offset=10.**-exponent;t=event['fraction']+direction*offset;name=f'{label}-{exponent}'
            refined=local_pose(p,q,t,cell,cut);coarse=local_pose(p,q,t,cell,None)
            fine,fine_audit=native_geometry(name+'-refined',refined,args.geometry_binary,args.output)
            unrefined,coarse_audit=native_geometry(name+'-coarse',coarse,args.geometry_binary,args.output)
            split=partition_control(name,unrefined,fine,args.partition_binary,args.output)
            response=pressure_control(name+'-pressure',refined,fine)
            samples.append(dict(offset_fraction=offset,time=refined['time'],refined_geometry=fine_audit,coarse_geometry=coarse_audit,partition=split,pressure=response))
            print(json.dumps(dict(event=label,offset=offset,area=response['opening_m2'],right_pressure=response['pressure_Pa'][1],merged_pressure=response['instantaneously_merged_pressure_Pa'])),flush=True)
        areas=[s['pressure']['opening_m2'] for s in samples];front=[s['pressure']['pressure_Pa'][1] for s in samples]
        assert all(a>b>0 for a,b in zip(areas,areas[1:])) and all(a>b>0 for a,b in zip(front,front[1:]))
        assert areas[-1]<areas[0]*1.01e-5 and front[-1]<front[0]*1.1e-5
        closed=local_pose(p,q,event['fraction']-direction*1e-6,cell,cut)
        native,audit=native_geometry(label+'-closed',closed,args.geometry_binary,args.output)
        assert len(native['regions'])==2 and native['links']==[]
        v=np.array([r[1] for r in native['regions']]);initial=np.array([600.,0.])
        pressure,flows,receipt=pressure_reference.run(label+'-closed-pressure',v,v/BULK,initial,np.zeros(2),[],DT)
        assert np.array_equal(pressure,initial) and len(flows)==0
        results.append(dict(event=event,samples=samples,closed_geometry=audit,closed_pressure=receipt))
    sources=['include/ducklips/cavity_network.hpp','src/cavity_network.cpp','tests/cavity_network_tests.cpp',
             'tools/cavity_region_motion_probe.cpp','tools/planar_contact_partition_reference.py','tools/planar_regions_reference.py',
             'tools/cavity_graph_reference.py','include/ducklips/planar_airspace.hpp','src/planar_airspace.cpp']
    report=dict(use='two located Foxtrot contacts, conservative pressure partition and frozen local-neck controls; no moving speech integration',
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],bulk_modulus_Pa=BULK,
                mobility_m2_per_Pa_s=MOBILITY,frozen_pressure_step_s=DT,events=results,
                input_poses={str(path.relative_to(ROOT)):sha(path) for path in input_paths},
                hashes={str(path.relative_to(ROOT)):sha(path) for path in [*(ROOT/s for s in sources),args.geometry_binary,args.pressure_binary,args.partition_binary]},
                versions=dict(numpy=np.__version__,shapely=shapely.__version__,geos=shapely.geos_version_string))
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print('PASS two contacts, twelve open-neck samples, two sealed controls and twelve native partitions',flush=True)


if __name__=='__main__':main()
