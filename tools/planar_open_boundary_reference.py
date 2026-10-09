#!/usr/bin/env python3
"""Pressure-release end caps and imposed inlet pressure on the planar replay.

Neither cap is a completed lip/glottal model; no radiation losses are present.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from planar_inertive_motion_reference import ROOT,run,audit_motion,PORTS,INITIAL
from planar_region_motion_reference import swept_boundaries


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load_cases(p,q):
    mouth=p['boundary_edges']['mouth_end'];inlet=p['boundary_edges']['supraglottal_entry']
    mouth_length=p['domain'][mouth][0]-.5*(p['x_edges'][1]+p['x_edges'][2])
    inlet_length=.5*(p['y_edges'][0]+p['y_edges'][1])-p['domain'][inlet][1]
    assert mouth_length>0 and inlet_length>0
    mouth_loads=[(mouth,3,mouth_length,0.),(mouth,4,mouth_length,0.)]
    return dict(mouth_open=mouth_loads,both_open=mouth_loads+[(inlet,0,inlet_length,600.)])


def clearance(p,q):
    # Convex hulls of all corresponding solid-edge sweeps conservatively
    # contain those moving surfaces. Include initial/final solid interiors.
    empty=dict(p,domain=p['solids'][0],solids=p['solids'][1:])
    other=dict(q,domain=q['solids'][0],solids=q['solids'][1:])
    occupied=unary_union([swept_boundaries(empty,other),*[Polygon(c) for pose in [p,q] for c in pose['solids']]])
    result={}
    for name,edge in p['boundary_edges'].items():
        assert q['boundary_edges'][name]==edge
        corners=[pose['domain'][k] for pose in [p,q] for k in [edge,(edge+1)%len(pose['domain'])]]
        cap=shapely.convex_hull(shapely.multipoints(corners));distance=cap.distance(occupied)
        assert distance>0,(name,distance)
        result[name]=distance
    return result


def controls(binary,out):
    results=[]
    # Closed inlet / ideal pressure-release outlet: quarter-wave duct.
    for n in [8,16,32]:
        length=.12;height=.02;dx=length/n;duration=.00037
        p=dict(time=0.,depth=.022,domain=[[0,0],[length,0],[length,height],[0,height]],solids=[],
               x_edges=np.linspace(0,length,n+1).tolist(),y_edges=[0,height])
        q=dict(p,time=duration);initial=100*np.cos(np.pi*(np.arange(n)+.5)/(2*n));ports=[(i,i+1,i+1,True) for i in range(n-1)]
        omega=2*350/dx*np.sin(np.pi/(4*n));previous=None
        for subdivisions in [16,32,64]:
            native,receipt=run(f'open-tube-{n}',p,q,subdivisions,initial,ports,binary,out,[(1,n-1,dx/2,0.)])
            error=float(max(abs(np.array(native['steps'][-1]['pressure'])-initial*np.cos(omega*duration))))
            e0=native['initial']['energy'];energy=max(abs(s['energy']/e0-1) for s in native['steps'])
            assert energy<1e-12
            if previous is not None:assert previous/error>14
            previous=error
            receipt.update(exact_discrete_pressure_error_Pa=error,maximum_relative_energy_error=energy,
                           discrete_frequency_Hz=omega/(2*np.pi),continuum_frequency_Hz=350/(4*length))
            results.append(receipt)
        assert error<1e-5
    # One compliant volume driven through an inertive cap by held external
    # pressure: closed-form forced oscillator fixes sign and source work.
    p=dict(time=0.,depth=.022,domain=[[0,0],[.04,0],[.04,.02],[0,.02]],solids=[],x_edges=[0,.04],y_edges=[0,.02])
    duration=.00005;q=dict(p,time=duration);pressure=600.;omega=350*np.sqrt(2)/.04
    for subdivisions in [8,16,32]:
        native,receipt=run('driven-cap',p,q,subdivisions,[0.],[],binary,out,[(1,0,.02,pressure)])
        final=native['steps'][-1];expected=pressure*(1-np.cos(omega*duration))
        error=abs(final['pressure'][0]-expected)
        work=-sum(pressure*s['average_flow'][0]*s['dt'] for s in native['steps'])
        energy_error=abs(final['energy']-work)
        assert final['pressure'][0]>0 and final['flow'][0]<0 and energy_error<1e-16
        receipt.update(exact_pressure_error_Pa=error,source_work_J=work,source_energy_error_J=energy_error);results.append(receipt)
    assert error<1e-6
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_planar_inertive_motion_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/planar-open-boundary')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_open_boundary_reference.json')
    parser.add_argument('--sanitizer',action='store_true')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    analytic=controls(args.binary,args.output)
    contact_path=ROOT/'reports/contact_partition_reference.json';contacts=json.loads(contact_path.read_text())
    runs=[];refinement=[];clearances={};inputs={}
    for name,start,contact in zip(['closure','release'],[5.031,5.062],contacts['events']):
        paths=[ROOT/f'out/planar-boundary/poses/pose-{t:.3f}.json' for t in [start,start+.001]]
        p,q=[json.loads(path.read_text()) for path in paths];event=contact['event']
        for pose in [p,q]:pose.update(x_edges=[pose['x_edges'][0],event['point_m'][0],pose['x_edges'][-1]],y_edges=[pose['y_edges'][0],.09,pose['y_edges'][-1]])
        clearances[name]=clearance(p,q)
        for case,loads in load_cases(p,q).items():
            differences=[];previous=None
            for subdivisions in ([16] if args.sanitizer else [16,32,64,128]):
                receipt=audit_motion(name,p,q,event,subdivisions,args.binary,args.output,loads,case);runs.append(receipt)
                if previous is not None:differences.append(dict(subdivisions=subdivisions,maximum_endpoint_difference_Pa=float(max(abs(np.array(receipt['final_pressure_Pa'])-previous)))))
                previous=np.array(receipt['final_pressure_Pa'])
            refinement.append(dict(name=name,case=case,comparisons=differences))
        # Default mode must still reproduce the previous sealed control.
        if not args.sanitizer:
            fresh,receipt=run('sealed-regression-'+name,p,q,128,INITIAL,PORTS,args.binary,args.output)
            old_path=ROOT/f'out/planar-inertive-motion/{name}-128.json';old=json.loads(old_path.read_text())
            worst=max(abs(a-b) for new,prior in zip(fresh['steps'],old['steps'],strict=True) for k in ['charge','momentum','pressure','flow','average_flow'] for a,b in zip(new[k],prior[k],strict=True))
            assert worst==0,worst
            receipt.update(previous_output=str(old_path.relative_to(ROOT)),previous_output_sha256=sha(old_path),maximum_state_difference=worst)
            analytic.append(receipt)
        inputs.update({str(path.relative_to(ROOT)):sha(path) for path in paths})
    sources=['tools/planar_open_boundary_reference.py','tools/planar_inertive_motion_reference.py','tools/planar_boundary_reference.py',
             'tools/planar_inertive_motion_probe.cpp','tools/planar_pressure_motion_reference.py','tools/planar_contact_partition_reference.py',
             'tools/planar_regions_reference.py','tools/planar_region_motion_reference.py','src/planar_airspace.cpp','include/ducklips/planar_airspace.hpp',
             'src/planar_face_motion.cpp','include/ducklips/planar_face_motion.hpp','src/moving_tree.cpp','include/ducklips/moving_tree.hpp']
    report=dict(use='ideal pressure boundaries at tongue/floor domain caps; no anatomical lip/glottal constriction or radiation losses',
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],runs=runs,controls=analytic,refinement=refinement,
                conservative_cap_to_solid_sweep_clearance_m=clearances,input_poses=inputs,contact_receipt=str(contact_path.relative_to(ROOT)),contact_sha256=sha(contact_path),
                hashes={str(path.relative_to(ROOT)):sha(path) for path in [*(ROOT/s for s in sources),args.binary]})
    args.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(refinement,indent=2),flush=True)


if __name__=='__main__':main()
