#!/usr/bin/env python3
"""Audit moving-region correspondence using space clear throughout a step.

A positive witness establishes geometric continuation, not material parcel
identity. Ambiguous split/merge relations and missing witnesses are explicit;
this tool does not turn an overlap fraction into a mass-transfer weight.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import shapely
from shapely.geometry import Polygon, Point, box
from shapely.ops import unary_union

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fluid(p):
    polygons=[Polygon(c) for c in [p['domain'],*p['solids']]]
    assert all(g.is_valid and g.area>0 for g in polygons)
    return polygons[0].difference(unary_union(polygons[1:]))


def swept_boundaries(first,last):
    """Conservative enclosure of every boundary point along linear motion.

    A moving edge point is a convex combination of its four endpoint corners.
    Their hull encloses its ruled sweep, including rotating/deforming edges.
    A hull can overestimate that sweep; this makes witnesses conservative.
    It is exact for a translating straight edge, up to floating arithmetic.
    """
    quads=[]
    for a,b in zip([first['domain'],*first['solids']],[last['domain'],*last['solids']],strict=True):
        a,b=np.asarray(a),np.asarray(b);assert a.shape==b.shape
        quads.extend(np.stack([a,np.roll(a,-1,axis=0),np.roll(b,-1,axis=0),b],axis=1))
    hulls=shapely.convex_hull(shapely.multipoints(np.asarray(quads)))
    # Zero-area sweeps cannot remove a positive-area witness; excluding them
    # from the area subtraction also avoids mixed
    # dimensional Boolean overlays for edge-touching empty intersections.
    return unary_union(hulls[shapely.area(hulls)>0])


def exact_regions(p,native):
    """Match independently clipped fragments to native IDs by interior points."""
    nx,ny=len(p['x_edges'])-1,len(p['y_edges'])-1
    cells={}
    for i,(cell,volume,x,y) in enumerate(native['regions']):cells.setdefault(cell,[]).append(i)
    ordered=sorted(cells)
    squares=[box(p['x_edges'][cell%nx],p['y_edges'][cell//nx],p['x_edges'][cell%nx+1],p['y_edges'][cell//nx+1]) for cell in ordered]
    clipped=shapely.intersection(squares,fluid(p));geometries=[None]*len(native['regions']);error=0
    for cell,geometry in zip(ordered,clipped):
        parts=[g for g in shapely.get_parts(geometry) if g.area>0]
        ids=cells[cell];assert len(parts)==len(ids),(p['time'],cell,len(parts),len(ids))
        used=set()
        for i in ids:
            _,volume,x,y=native['regions'][i]
            matches=[k for k,g in enumerate(parts) if g.covers(Point(x,y))]
            assert len(matches)==1 and matches[0] not in used
            used.add(matches[0]);geometries[i]=parts[matches[0]]
            error=max(error,abs(geometries[i].area*p['depth']-volume))
    assert error<1e-18
    return geometries,cells,error


def correspondence(first,last,old_geometry,old_cells,new_geometry,new_cells):
    assert first['depth']==last['depth'] and first['x_edges']==last['x_edges'] and first['y_edges']==last['y_edges']
    pairs=[(i,j) for cell,ids in old_cells.items() for i in ids for j in new_cells.get(cell,[])]
    old=[old_geometry[i] for i,j in pairs];new=[new_geometry[j] for i,j in pairs]
    overlap=shapely.intersection(old,new)
    swept=swept_boundaries(first,last)
    witness=shapely.difference(overlap,swept)
    overlap_volume=shapely.area(overlap)*first['depth'];witness_volume=shapely.area(witness)*first['depth']
    links=[];removed=[];old_targets={};new_sources={};max_subset_error=0
    for (i,j),raw,clear,geometry in zip(pairs,overlap_volume,witness_volume,witness):
        max_subset_error=max(max_subset_error,float(clear-raw))
        if clear>0:
            point=geometry.representative_point()
            assert old_geometry[i].covers(point) and new_geometry[j].covers(point) and not swept.contains(point)
            links.append([i,j,float(raw),float(clear),point.x,point.y])
            old_targets.setdefault(i,[]).append(j);new_sources.setdefault(j,[]).append(i)
        elif raw>0:removed.append([i,j,float(raw)])
    # A thin fragment's area can be ill-conditioned relative to the absolute
    # coordinates used by Boolean intersections. Check an explicit coordinate
    # roundoff allowance (boundary length times positional roundoff), rather
    # than an absolute m^3 floor or a relative tolerance on the tiny area.
    # This is a consistency allowance, not a certified GEOS error bound, and
    # it changes no geometry, witness volume, or correspondence classification.
    subset_allowance=np.asarray([64*np.finfo(float).eps*first['depth']*
        max(abs(x) for x in (*a.bounds,*b.bounds))*(a.length+b.length) for a,b in zip(old,new)])
    assert np.all((overlap_volume>0)|(witness_volume==0))
    assert np.all(witness_volume-overlap_volume<=subset_allowance)
    # Independent point-in-fluid checks along the motion guard the constructive
    # sweep argument against implementation/rounding mistakes. They supplement
    # the hull enclosure; finite samples alone do not certify a clear path.
    checked=0
    if links:
        x=np.asarray([row[4] for row in links]);y=np.asarray([row[5] for row in links])
        for fraction in [.25,.5,.75]:
            contours=[]
            for a,b in zip([first['domain'],*first['solids']],[last['domain'],*last['solids']]):
                a,b=np.asarray(a),np.asarray(b);contours.append((a+fraction*(b-a)).tolist())
            sample=dict(first,domain=contours[0],solids=contours[1:])
            clear=shapely.contains_xy(fluid(sample),x,y)
            assert np.all(clear),(first['time'],last['time'],fraction,np.flatnonzero(~clear).tolist())
            checked+=len(links)
    # These are graph correspondences only. They must not be normalized into
    # mass weights: swept-out old air still needs an interval flow equation.
    continuation=[[i,ids[0]] for i,ids in old_targets.items() if len(ids)==1 and len(new_sources[ids[0]])==1]
    return dict(time_first=first['time'],time_last=last['time'],links=links,
                one_to_one=continuation,
                split_candidates=[[i,ids] for i,ids in old_targets.items() if len(ids)>1],
                merge_candidates=[[j,ids] for j,ids in new_sources.items() if len(ids)>1],
                old_without_witness=[i for i in range(len(old_geometry)) if i not in old_targets],
                new_without_witness=[i for i in range(len(new_geometry)) if i not in new_sources],
                endpoint_overlaps_removed=removed,maximum_subset_error_m3=max_subset_error,intermediate_witness_checks=checked)


def fixture(domain,solids,time):
    return dict(time=time,domain=list(domain.exterior.coords)[:-1],solids=[list(s.exterior.coords)[:-1] for s in solids],
                depth=1.,x_edges=[-2.,2.],y_edges=[-1.,1.])


def fixture_regions(p):
    parts=sorted([g for g in shapely.get_parts(fluid(p)) if g.area>0],key=lambda g:g.centroid.x)
    return parts,{0:list(range(len(parts)))}


def controls():
    mouth=box(-2,-1,2,1)
    first=fixture(mouth,[box(-.6,-2,-.4,2)],0)
    last=fixture(mouth,[box(.4,-2,.6,2)],1)
    a,ac=fixture_regions(first);b,bc=fixture_regions(last)
    moving=correspondence(first,last,a,ac,b,bc)
    assert moving['one_to_one']==[[0,0],[1,1]]
    assert len(moving['endpoint_overlaps_removed'])==1
    removed=moving['endpoint_overlaps_removed'][0]
    assert removed[:2]==[1,0] and abs(removed[2]-1.6)<1e-14
    # Assign different passive tracer labels. Endpoint overlap says some old
    # right-chamber space belongs to the new left chamber, despite the seal.
    # Witness correspondence never creates that edge.
    assert all(i==j for i,j,*_ in moving['links'])
    # If an entire pocket translates farther than its own width, a stationary
    # witness need not exist even though continuation is real. Do not invent it.
    small0=fixture(box(-1,-.5,-.5,.5),[],0)
    small1=fixture(box(.5,-.5,1,.5),[],1)
    c,cc=fixture_regions(small0);d,dc=fixture_regions(small1)
    missing=correspondence(small0,small1,c,cc,d,dc)
    assert missing['old_without_witness']==[0] and missing['new_without_witness']==[0] and not missing['links']
    # A wall enters from above and divides a cell. The relation is one-to-many,
    # not an instruction to choose one target or immediately divide its mass.
    open_pose=fixture(mouth,[box(-.1,.5,.1,2)],0)
    closed_pose=fixture(mouth,[box(-.1,-2,.1,2)],1)
    e,ec=fixture_regions(open_pose);f,fc=fixture_regions(closed_pose)
    split=correspondence(open_pose,closed_pose,e,ec,f,fc)
    assert len(e)==1 and len(f)==2 and split['split_candidates']==[[0,[0,1]]]
    reverse=correspondence(closed_pose,open_pose,f,fc,e,ec)
    assert reverse['merge_candidates']==[[0,[0,1]]]
    # Positive one-to-one endpoint witnesses do NOT certify unchanged topology
    # throughout the step. A shutter crosses a narrow joining corridor while
    # both large chambers retain clear witness points.
    neck=unary_union([box(-2,-1,-.25,1),box(.25,-1,2,1),box(-.25,-.1,.25,.1)])
    gate0=fixture(neck,[box(-.1,.2,.1,.4)],0)
    gate1=fixture(neck,[box(-.1,-.4,.1,-.2)],1)
    g,gc=fixture_regions(gate0);h,hc=fixture_regions(gate1)
    transient=correspondence(gate0,gate1,g,gc,h,hc)
    middle=fixture(neck,[box(-.1,-.1,.1,.1)],.5)
    mid_regions,_=fixture_regions(middle)
    assert len(g)==len(h)==1 and transient['one_to_one']==[[0,0]] and len(mid_regions)==2
    transient['middle_fluid_components']=len(mid_regions)
    transient['whole_interval_topology_certified']=False
    return dict(translating_sealed_wall=moving,translated_pocket_without_witness=missing,
                dividing_wall=split,retracting_wall=reverse,temporary_internal_closure=transient)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stride',type=int,default=1,help='stride through saved 125 microsecond poses')
    parser.add_argument('--output',type=Path,default=ROOT/'out/planar-region-motion')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_region_motion_reference.json')
    args=parser.parse_args();assert args.stride in [1,2,4,8]
    args.output.mkdir(parents=True,exist_ok=True)
    control=controls()
    receipt_path=ROOT/'reports/planar_motion_reference.json';receipt=json.loads(receipt_path.read_text())
    cases=receipt['samples'][::args.stride];assert len(cases)>1
    states=[];max_error=0;results=[]
    old=None
    for number,case in enumerate(cases):
        metric_path=ROOT/case['metrics'];native_path=ROOT/case['output']
        assert sha(metric_path)==case['metrics_sha256'] and sha(native_path)==case['output_sha256']
        p=json.loads(metric_path.read_text());native=json.loads(native_path.read_text())
        geometries,cells,error=exact_regions(p,native);max_error=max(max_error,error)
        states.append(dict(time=p['time'],metrics=case['metrics'],metrics_sha256=case['metrics_sha256'],native=case['output'],native_sha256=case['output_sha256']))
        if old is not None:
            first,a,ac=old
            result=correspondence(first,p,a,ac,geometries,cells)
            destination=args.output/f'interval-{number-1:04}.json';destination.write_text(json.dumps(result,separators=(',',':'))+'\n')
            summary={k:v for k,v in result.items() if k not in ['links','one_to_one']}
            summary.update(positive_witness_pairs=len(result['links']),one_to_one_pairs=len(result['one_to_one']),
                           output=str(destination.relative_to(ROOT)),output_sha256=sha(destination))
            results.append(summary)
            if number%16==0 or result['split_candidates'] or result['merge_candidates']:
                print(json.dumps({k:summary[k] for k in ['time_last','one_to_one_pairs','split_candidates','merge_candidates','old_without_witness','new_without_witness']}),flush=True)
        old=p,geometries,cells
    report=dict(use='conservative geometric continuation witnesses, not moving mass transfer',
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],stride=args.stride,
                controls=control,intervals=results,inputs=states,maximum_region_volume_error_m3=max_error,
                source_receipt=str(receipt_path.relative_to(ROOT)),source_receipt_sha256=sha(receipt_path),
                hashes={'tools/planar_region_motion_reference.py':sha(Path(__file__).resolve())},
                versions=dict(numpy=np.__version__,shapely=shapely.__version__,geos=shapely.geos_version_string))
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',len(results),'intervals; geometric correspondence only',flush=True)


if __name__=='__main__':main()
