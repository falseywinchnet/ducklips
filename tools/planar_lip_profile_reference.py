#!/usr/bin/env python3
"""Integrated axial lip storage coupled to moving planar mouth regions.

Two separately initialized, continuous two-interval Bravo controls; no audio.
"""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import shapely
from shapely.geometry import LineString
from lip_profile_reference import oracle as lip_metrics, vector, FIELDS
from planar_inertive_motion_reference import ROOT, RHO, BULK, write_input, dense_step, sha
from planar_pressure_motion_reference import sampled
from planar_regions_reference import check_geometry
from planar_region_motion_reference import exact_regions, correspondence, swept_boundaries
from planar_boundary_reference import check_boundaries
from planar_open_boundary_reference import clearance

PORTS=[(0,1,1,False)]


def run(name,p,q,profiles,subdivisions,pressure,previous,binary,out,compact):
    source=out/f'{name}-{subdivisions}.txt';destination=source.with_suffix('.json');log=source.with_suffix('.log')
    write_input(p,q,subdivisions,[pressure]*2,PORTS,source)
    mouth=p['boundary_edges']['mouth_end'];distance=.5*(p['x_edges'][-1]-p['x_edges'][0])
    attachments=[(mouth,i,distance) for i in range(2)]
    with source.open('a') as f:
        f.write('4 2\n');f.writelines(' '.join(map(str,v))+'\n' for v in profiles)
        f.write(' '.join(map(str,[pressure]*4))+'\n')
        f.writelines(' '.join(map(str,row))+'\n' for row in attachments)
        if previous is not None:
            f.writelines(' '.join(map(str,previous[k]))+'\n' for k in ['charge','momentum'])
    command=[str(binary),'--lip-profile']
    if previous is not None:command.append('--initial-state')
    if compact:command.append('--compact-geometry')
    with source.open() as f,destination.open('w') as g:result=subprocess.run(command,stdin=f,stdout=g,stderr=subprocess.PIPE,text=True)
    log.write_text(result.stderr);assert result.returncode==0,(name,result.stderr)
    native=json.loads(destination.read_text())
    if previous is not None:
        for k in ['volume','compliance','conductance','charge','momentum']:
            assert native['initial'][k]==previous[k],('continuation changed accepted state',k)
    receipt=dict(name=name,subdivisions=subdivisions,steps=len(native['steps']),command=command,attachments=attachments,
                 continued_from_previous=previous is not None,
                 **{k:str(v.relative_to(ROOT)) for k,v in [('input',source),('output',destination),('log',log)]},
                 **{k+'_sha256':sha(v) for k,v in [('input',source),('output',destination),('log',log)]})
    return native,receipt


def pose_at(p,q,record):
    pose=sampled(p,q,record['fraction']);contours=[]
    for a,b,c in zip([p['domain'],*p['solids']],[q['domain'],*q['solids']],[pose['domain'],*pose['solids']]):
        a=np.array(a);b=np.array(b);contours.append(np.where(a==b,a,c).tolist())
    pose.update(domain=contours[0],solids=contours[1:])
    if 'domain' in record:pose.update(domain=record['domain'],solids=record['solids'])
    return pose


def audit(native,p,q,profiles,pressure,previous_oracle,full):
    n=6;mouth=p['boundary_edges']['mouth_end'];distance=.5*(p['x_edges'][-1]-p['x_edges'][0])
    line=LineString([(p['x_edges'][0],.09),(p['x_edges'][-1],.09)])
    core_distance=.5*(p['y_edges'][-1]-p['y_edges'][0])
    edges=[(0,1),(0,2),(1,2),(2,3),(3,4),(4,5),(5,-1)];b=np.zeros((n,len(edges)))
    for k,(i,j) in enumerate(edges):
        b[i,k]=1
        if j>=0:b[j,k]=-1
    maxima=dict(volume_m3=0.,cap_area_m2=0.,internal_area_m2=0.,conductance=0.,pressure_Pa=0.,momentum_Pa_s=0.,
                average_flow_m3_s=0.,local_mass_m3=0.,global_mass_m3=0.,sealed_lip_flow_m3_s=0.)
    geometry_checks=0;witnesses=0
    def geometry(record):
        nonlocal geometry_checks
        f=record['fraction'];pose=pose_at(p,q,record)
        checked=check_geometry(pose,record);regions,cells,error=exact_regions(pose,record)
        boundaries,errors=check_boundaries(pose,record,[mouth],regions)
        maxima['volume_m3']=max(maxima['volume_m3'],error)
        maxima['internal_area_m2']=max(maxima['internal_area_m2'],checked['maximum_port_area_error_m2'])
        maxima['cap_area_m2']=max(maxima['cap_area_m2'],errors['area_m2'])
        lip=lip_metrics((1-f)*profiles[0]+f*profiles[1]);v=np.r_[[r.area*p['depth'] for r in regions],lip['volume']]
        aperture=regions[0].intersection(regions[1]).intersection(line).length*p['depth']
        cap=np.array([boundaries[(mouth,i)][0] for i in range(2)]);shares=cap/sum(cap)
        entry=np.zeros(2) if lip['conductance'][0]==0 else 1/(RHO*distance/cap+1/(lip['conductance'][0]*shares))
        g=np.r_[aperture/(RHO*core_distance),entry,lip['conductance'][1:]]
        geometry_checks+=1
        return pose,regions,cells,v,v/BULK,g
    first=native['initial'];old=first;outflow=0.;mass0=sum(np.array(first['volume'])+first['charge']);density=[]
    continuous_open=None
    if full:
        old_pose,old_regions,old_cells,old_v,c,g=geometry(first['geometry'])
        z=np.r_[c*pressure,np.zeros(len(edges))] if previous_oracle is None else previous_oracle.copy()
        final_regions,_,_=exact_regions(q,native['steps'][-1]['geometry'])
        opening=old_regions[0].intersection(old_regions[1]).intersection(final_regions[0]).intersection(final_regions[1]).intersection(line)
        continuous_open=opening.difference(swept_boundaries(p,q)).length*p['depth'];assert continuous_open>0
    else:z=None
    p0,p1=profiles;expected=[]
    if p0[-1]==p1[-1]==0 and (p0[2]<1)!=(p1[2]<1):
        t=(1-p0[2])/(p1[2]-p0[2])
        if 0<t<1:expected.append(t)
    assert native['lip_contact_cuts']==expected and not native['contact_cuts']
    for state in native['steps']:
        dt=state['dt']
        if full:
            _,_,_,_,c1,g1=geometry(state['stage1']);_,_,_,_,c2,g2=geometry(state['stage2'])
            pose,regions,cells,v,c,g=geometry(state['geometry'])
            match=correspondence(old_pose,pose,old_regions,old_cells,regions,cells)
            assert match['one_to_one']==[[0,0],[1,1]] and not any(match[k] for k in ['split_candidates','merge_candidates','old_without_witness','new_without_witness'])
            witnesses+=match['intermediate_witness_checks']
            z[n:][(g1==0)&(g2==0)]=0
            z,flow,_=dense_step(z,c1,c2,g1,g2,-(v-old_v)/dt,dt,b)
            z[n:][g==0]=0
            maxima['pressure_Pa']=max(maxima['pressure_Pa'],float(max(abs(np.array(state['pressure'])-z[:n]/c))))
            maxima['momentum_Pa_s']=max(maxima['momentum_Pa_s'],float(max(abs(np.array(state['momentum'])-z[n:]))))
            maxima['average_flow_m3_s']=max(maxima['average_flow_m3_s'],float(max(abs(np.array(state['average_flow'])-flow))))
            maxima['volume_m3']=max(maxima['volume_m3'],float(max(abs(np.array(state['volume'])-v))))
            maxima['conductance']=max(maxima['conductance'],float(max(abs(np.array(state['conductance'])-g))))
            old_pose,old_regions,old_cells,old_v=pose,regions,cells,v
        current=np.array(state['volume'])+state['charge'];prior=np.array(old['volume'])+old['charge']
        balance=current-prior+dt*b@np.array(state['average_flow'])
        maxima['local_mass_m3']=max(maxima['local_mass_m3'],float(max(abs(balance))))
        outflow+=dt*state['average_flow'][-1]
        maxima['global_mass_m3']=max(maxima['global_mass_m3'],abs(sum(current)+outflow-mass0))
        if state['conductance'][-1]==0 and old['conductance'][-1]==0:
            maxima['sealed_lip_flow_m3_s']=max(maxima['sealed_lip_flow_m3_s'],abs(state['average_flow'][-1]))
        density.extend((1+np.array(state['pressure'])/BULK).tolist());old=state
    assert maxima['local_mass_m3']<1e-17 and maxima['global_mass_m3']<1e-17 and maxima['sealed_lip_flow_m3_s']==0,maxima
    if full:assert maxima['pressure_Pa']<1e-3 and maxima['average_flow_m3_s']<1e-9 and maxima['conductance']<1e-12,maxima
    result=dict(independent_audit=full,geometry_checks=geometry_checks,intermediate_witness_checks=witnesses,maxima=maxima,
                internal_port_continuously_open_witness_m2=continuous_open,lip_contact_cuts=expected,
                final_pressure_Pa=native['steps'][-1]['pressure'],relative_density_range=[min(density),max(density)],
                minimum_lip_storage_m3=min(min(s['volume'][2:]) for s in [first,*native['steps']]),
                net_outlet_volume_m3=outflow)
    return result,z


def setup(start):
    paths=[ROOT/f'out/planar-lip-profile/poses/pose-{t:.3f}.json' for t in [start,start+.001]]
    p,q=[json.loads(path.read_text()) for path in paths]
    for pose in [p,q]:
        mouth=pose['boundary_edges']['mouth_end'];entry=pose['boundary_edges']['supraglottal_entry']
        pose.update(x_edges=[-.01,pose['domain'][mouth][0]],y_edges=[pose['domain'][entry][1],.09,.105])
    return p,q,paths


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_planar_inertive_motion_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/planar-lip-profile/audit')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_lip_profile_reference.json')
    parser.add_argument('--sanitizer',action='store_true');args=parser.parse_args()
    args.binary=args.binary.resolve();args.output=args.output.resolve();args.report=args.report.resolve();args.output.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((ROOT/'out/lip-profile/nato.csv').open()));runs=[];refinement=[];input_hashes={};clearances={}
    for name,start,pressure in [('closing',1284,0.),('opening',1374,700.)]:
        previous_pressure=None;changes=[]
        for subdivisions in ([16] if args.sanitizer else [16,32,64,128,256,512,1024,2048,4096,8192]):
            previous=None;previous_oracle=None;full=subdivisions<=64
            for interval in range(2):
                index=start+interval;p,q,paths=setup(index/1000);profiles=[vector(rows[i]) for i in [index,index+1]]
                clearances[f'{name}-{interval}']=clearance(p,q)
                native,receipt=run(f'{name}-{interval}',p,q,profiles,subdivisions,pressure,previous,args.binary,args.output,compact=not full)
                checked,previous_oracle=audit(native,p,q,profiles,pressure,previous_oracle,full)
                receipt.update(checked);runs.append(receipt);previous=native['steps'][-1]
                input_hashes.update({str(path.relative_to(ROOT)):sha(path) for path in paths})
            final=np.array(previous['pressure'])
            if previous_pressure is not None:changes.append(dict(subdivisions=subdivisions,maximum_endpoint_difference_Pa=float(max(abs(final-previous_pressure)))))
            previous_pressure=final
            print(json.dumps(dict(name=name,subdivisions=subdivisions,final_pressure_Pa=final.tolist(),change=None if not changes else changes[-1],maxima=runs[-1]['maxima'])),flush=True)
        refinement.append(dict(name=name,changes=changes))
    sources=['tools/planar_lip_profile_reference.py','tools/planar_inertive_motion_probe.cpp','tools/airspace_geometry_probe.cpp',
             'tools/lip_profile_reference.py','src/lip_profile.cpp','include/ducklips/lip_profile.hpp','src/moving_graph.cpp',
             'include/ducklips/moving_graph.hpp','src/planar_airspace.cpp','src/planar_face_motion.cpp','tools/planar_inertive_motion_reference.py',
             'tools/planar_regions_reference.py','tools/planar_region_motion_reference.py','tools/planar_pressure_motion_reference.py',
             'tools/planar_contact_partition_reference.py','tools/planar_boundary_reference.py','tools/planar_open_boundary_reference.py',
             'bodies/reference.body','out/nato-airspace-geometry.csv','out/lip-profile/nato.csv']
    report=dict(use='coupled planar mouth and integrated lip geometry, with interval continuation; no new speech',
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],runs=runs,refinement=refinement,
                input_poses=input_hashes,cap_clearance_m=clearances,
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),args.binary,ROOT/'build/ducklips_airspace_geometry_probe']},
                versions=dict(numpy=np.__version__,shapely=shapely.__version__,geos=shapely.geos_version_string))
    args.report.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
