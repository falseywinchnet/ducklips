#!/usr/bin/env python3
"""Native moving-pressure replay of the two complete Foxtrot contact intervals.

Five spatial pressure regions, sealed exterior and prescribed linear mobility.
Independent Boolean volumes, temporal apertures, state correspondence and a
separately evolved linear pressure solve. This is not speech or calibrated flow.
"""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import scipy
from scipy.integrate import quad
import shapely
from shapely.geometry import Polygon,LineString
from shapely.ops import unary_union
from planar_regions_reference import check_geometry
from planar_region_motion_reference import exact_regions,correspondence
from planar_contact_partition_reference import contours

ROOT=Path(__file__).resolve().parents[1]
BULK=1.14*350**2
MOBILITY=1e-7
PORTS=[(0,1,2,False),(1,4,4,True),(2,3,1,True),(3,4,3,False)]
INITIAL=np.array([600.,600.,0.,0.,0.])


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def sampled(p,q,fraction):
    c=contours(p,q,fraction)
    return dict(p,time=p['time']+fraction*(q['time']-p['time']),domain=c[0],solids=c[1:])


def write_input(p,q,subdivisions,path):
    lines=[f'2 2 {len(p["solids"])} {p["depth"]} {q["time"]-p["time"]:.17g} {subdivisions}',
           ' '.join(map(str,p['x_edges'])),' '.join(map(str,p['y_edges']))]
    for pose in [p,q]:
        for poly in [pose['domain'],*pose['solids']]:
            lines.append(str(len(poly)));lines.extend(' '.join(map(str,point)) for point in poly)
    lines+=['5 4',' '.join(map(str,INITIAL))]
    lines.extend(f'{i} {j} {face} {int(vertical)}' for i,j,face,vertical in PORTS)
    path.write_text('\n'.join(lines)+'\n')


def port_lines(p):
    lines=[];distances=[]
    for _,_,face,vertical in PORTS:
        if vertical:
            row,col=divmod(face,3);lines.append(LineString([(p['x_edges'][col],p['y_edges'][row]),(p['x_edges'][col],p['y_edges'][row+1])]))
            distances.append(.5*(p['x_edges'][col+1]-p['x_edges'][col-1]))
        else:
            row,col=divmod(face,2);lines.append(LineString([(p['x_edges'][col],p['y_edges'][row]),(p['x_edges'][col+1],p['y_edges'][row])]))
            distances.append(.5*(p['y_edges'][row+1]-p['y_edges'][row-1]))
    return lines,np.asarray(distances)


def audit(name,p,q,event,subdivisions,binary,out):
    source=out/f'{name}-{subdivisions}.txt';destination=source.with_suffix('.json');log=source.with_suffix('.log')
    write_input(p,q,subdivisions,source)
    with source.open() as f,destination.open('w') as g:
        result=subprocess.run([str(binary)],stdin=f,stdout=g,stderr=subprocess.PIPE,text=True)
    log.write_text(result.stderr);assert result.returncode==0,(name,subdivisions,result.stderr)
    native=json.loads(destination.read_text());states=native['states'];cuts=native['contact_cuts']
    assert cuts and max(abs(t-event['fraction']) for t in cuts)<1e-10,(name,cuts,event)
    assert states[0]['fraction']==0 and states[-1]['fraction']==1
    lines,distance=port_lines(p)
    @lru_cache(maxsize=None)
    def apertures(fraction):
        pose=sampled(p,q,fraction);polygons=[Polygon(c) for c in [pose['domain'],*pose['solids']]]
        assert all(g.is_valid for g in polygons)
        air=polygons[0].difference(unary_union(polygons[1:]))
        return np.array([line.intersection(air).length*pose['depth'] for line in lines])
    # No oracle pressure state is reset from a native intermediate result.
    oracle=INITIAL.copy();previous=None;records=[];mass0=None;sealed_mass=None
    maxima=dict(region_volume_error_m3=0.,port_geometry_error_m2=0.,interval_area_error_m2=0.,
                native_oracle_pressure_error_Pa=0.,total_mass_error_m3=0.,local_mass_error_m3=0.,
                sealed_group_mass_error_m3=0.,sealed_throat_flow_m3_per_s=0.)
    density=[];quadrature_evaluations=0
    for step,state in enumerate(states):
        t=state['fraction'];pose=sampled(p,q,t)
        pose.update(domain=state['domain'],solids=state['solids'])
        geometry=check_geometry(pose,state);regions,cells,error=exact_regions(pose,state)
        assert len(regions)==5 and [r[0] for r in state['regions']]==[0,2,0,1,3]
        exact_v=np.array([g.area*pose['depth'] for g in regions]);v=np.array([r[1] for r in state['regions']]);c=v/BULK
        pressure=np.array(state['pressure']);charge=v+c*pressure
        density.extend((1+pressure/BULK).tolist())
        maxima['region_volume_error_m3']=max(maxima['region_volume_error_m3'],error)
        maxima['port_geometry_error_m2']=max(maxima['port_geometry_error_m2'],geometry['maximum_port_area_error_m2'])
        if previous is None:
            mass0=sum(charge);assert np.array_equal(pressure,INITIAL)
        else:
            old_pose,old_regions,old_cells,old_exact_v,old_v,old_pressure,old_t=previous
            match=correspondence(old_pose,pose,old_regions,old_cells,regions,cells)
            assert sorted(match['one_to_one'])==[[i,i] for i in range(5)]
            assert not match['split_candidates'] and not match['merge_candidates'] and not match['old_without_witness'] and not match['new_without_witness']
            width=t-old_t;dt=width*(q['time']-p['time']);area=np.array(state['interval_area'])
            independent=[];quadrature_attempts=[]
            for edge in range(4):
                for divisions in [4,16,64,256,1024]:
                    integral,err,info,*warning=quad(lambda f:apertures(float(f))[edge],old_t,t,epsabs=1e-13*width,epsrel=1e-10,
                                                   points=np.linspace(old_t,t,divisions+1),limit=2*divisions+512,full_output=True)
                    quadrature_evaluations+=info['neval']
                    quadrature_attempts.append(dict(edge=edge,divisions=divisions,average_area_m2=integral/width,
                                                    error_estimate_m2=err/width,warning=warning))
                    if not warning:break
                assert not warning,(name,subdivisions,step,edge,warning)
                independent.append(integral/width)
            independent=np.asarray(independent)
            maxima['interval_area_error_m2']=max(maxima['interval_area_error_m2'],float(max(abs(area-independent))))
            g=MOBILITY*independent/distance;matrix=np.diag(exact_v/BULK/dt)
            for (i,j,_,_),conductance in zip(PORTS,g):
                matrix[i,i]+=conductance;matrix[j,j]+=conductance;matrix[i,j]-=conductance;matrix[j,i]-=conductance
            right=(old_exact_v/BULK*oracle-(exact_v-old_exact_v))/dt
            oracle=np.linalg.solve(matrix,right)
            maxima['native_oracle_pressure_error_Pa']=max(maxima['native_oracle_pressure_error_Pa'],float(max(abs(pressure-oracle))))
            balance=charge-(old_v+old_v/BULK*old_pressure)
            for (i,j,_,_),flow in zip(PORTS,state['flow']):balance[i]+=dt*flow;balance[j]-=dt*flow
            maxima['local_mass_error_m3']=max(maxima['local_mass_error_m3'],float(max(abs(balance))))
            closed_interval=(old_t>event['fraction']+1e-10) if name=='closure' else (t<event['fraction']-1e-10)
            if closed_interval:
                maxima['sealed_throat_flow_m3_per_s']=max(maxima['sealed_throat_flow_m3_per_s'],abs(state['flow'][1]))
                groups=np.array([sum(charge[:2]),sum(charge[2:])])
                if sealed_mass is None:sealed_mass=np.array([sum((old_v+old_v/BULK*old_pressure)[:2]),sum((old_v+old_v/BULK*old_pressure)[2:])])
                maxima['sealed_group_mass_error_m3']=max(maxima['sealed_group_mass_error_m3'],float(max(abs(groups-sealed_mass))))
            records.append(dict(fraction_first=old_t,fraction_last=t,independent_interval_area_m2=independent.tolist(),
                                oracle_pressure_Pa=oracle.tolist(),witness_pairs=match['one_to_one'],intermediate_witness_checks=match['intermediate_witness_checks'],
                                quadrature_attempts=quadrature_attempts))
        maxima['total_mass_error_m3']=max(maxima['total_mass_error_m3'],abs(sum(charge)-mass0))
        previous=pose,regions,cells,exact_v,v,pressure,t
    assert maxima['region_volume_error_m3']<1e-16 and maxima['port_geometry_error_m2']<1e-14
    assert maxima['interval_area_error_m2']<2e-12 and maxima['native_oracle_pressure_error_Pa']<1e-3,maxima
    assert maxima['total_mass_error_m3']<1e-17 and maxima['local_mass_error_m3']<1e-17 and maxima['sealed_group_mass_error_m3']<1e-17,maxima
    assert maxima['sealed_throat_flow_m3_per_s']==0,maxima
    audit_path=out/f'{name}-{subdivisions}-audit.json';audit_path.write_text(json.dumps(records,indent=2)+'\n')
    receipt=dict(name=name,subdivisions=subdivisions,steps=len(states)-1,contact_cuts=cuts,maxima=maxima,
                 relative_density_range=[min(density),max(density)],final_pressure_Pa=states[-1]['pressure'],
                 independent_quadrature_evaluations=quadrature_evaluations,command=[str(binary)],
                 input=str(source.relative_to(ROOT)),input_sha256=sha(source),output=str(destination.relative_to(ROOT)),output_sha256=sha(destination),
                 audit=str(audit_path.relative_to(ROOT)),audit_sha256=sha(audit_path),log=str(log.relative_to(ROOT)),log_sha256=sha(log))
    print(json.dumps({k:receipt[k] for k in ['name','subdivisions','steps','maxima','relative_density_range']}),flush=True)
    return receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_planar_pressure_motion_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/planar-pressure-motion')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_pressure_motion_reference.json')
    parser.add_argument('--sanitizer',action='store_true')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    contacts_path=ROOT/'reports/contact_partition_reference.json';contacts=json.loads(contacts_path.read_text())
    runs=[];refinement=[];inputs={}
    for name,start,contact in zip(['closure','release'],[5.031,5.062],contacts['events']):
        paths=[ROOT/f'out/planar-airspace/motion/pose-{t:.3f}.json' for t in [start,start+.001]]
        p,q=[json.loads(path.read_text()) for path in paths];event=contact['event']
        for pose in [p,q]:pose.update(x_edges=[pose['x_edges'][0],event['point_m'][0],pose['x_edges'][-1]],y_edges=[pose['y_edges'][0],.09,pose['y_edges'][-1]])
        previous=None;differences=[]
        for subdivisions in ([1,16] if args.sanitizer else [1,2,4,8,16,32,64,128]):
            result=audit(name,p,q,event,subdivisions,args.binary,args.output);runs.append(result)
            if previous is not None:differences.append(dict(subdivisions=subdivisions,maximum_endpoint_difference_Pa=float(max(abs(np.asarray(result['final_pressure_Pa'])-previous)))))
            previous=np.asarray(result['final_pressure_Pa'])
        refinement.append(dict(name=name,comparisons=differences))
        inputs.update({str(path.relative_to(ROOT)):sha(path) for path in paths})
    sources=['tools/planar_pressure_motion_probe.cpp','tools/planar_pressure_motion_reference.py',
             'include/ducklips/planar_airspace.hpp','src/planar_airspace.cpp','include/ducklips/planar_face_motion.hpp','src/planar_face_motion.cpp',
             'include/ducklips/cavity_network.hpp','src/cavity_network.cpp','tools/planar_regions_reference.py',
             'tools/planar_region_motion_reference.py','tools/planar_contact_partition_reference.py']
    report=dict(use='two separately initialized full-mouth moving pressure intervals on a coarse five-region graph; no speech or calibrated mobility',
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],runs=runs,refinement=refinement,input_poses=inputs,
                source_contacts=str(contacts_path.relative_to(ROOT)),source_contacts_sha256=sha(contacts_path),
                bulk_modulus_Pa=BULK,mobility_m2_per_Pa_s=MOBILITY,initial_pressure_Pa=INITIAL.tolist(),
                boundaries='sealed exterior; no glottal or lip port',
                hashes={str(path.relative_to(ROOT)):sha(path) for path in [*(ROOT/s for s in sources),args.binary]},
                versions=dict(numpy=np.__version__,scipy=scipy.__version__,shapely=shapely.__version__,geos=shapely.geos_version_string))
    args.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(refinement,indent=2),flush=True)


if __name__=='__main__':main()
