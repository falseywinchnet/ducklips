#!/usr/bin/env python3
"""Inviscid moving-airspace reference: geometry, mass, momentum and wave limits.

This replaces the earlier arbitrary Darcy mobility only in a diagnostic.
Five pressure regions and grid-centre inertances do not resolve a human mouth.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import shapely
from planar_pressure_motion_reference import sampled,port_lines,PORTS,BULK,INITIAL
from planar_regions_reference import check_geometry
from planar_region_motion_reference import exact_regions,correspondence
from planar_boundary_reference import check_boundaries

ROOT=Path(__file__).resolve().parents[1]
RHO=1.14
A=np.array([[.25,.25-np.sqrt(3)/6],[.25+np.sqrt(3)/6,.25]])


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def write_input(p,q,subdivisions,initial,ports,path):
    nx,ny=len(p['x_edges'])-1,len(p['y_edges'])-1
    lines=[f'{nx} {ny} {len(p["solids"])} {p["depth"]} {q["time"]-p["time"]:.17g} {subdivisions}',
           ' '.join(map(str,p['x_edges'])),' '.join(map(str,p['y_edges']))]
    for pose in [p,q]:
        for poly in [pose['domain'],*pose['solids']]:
            lines.append(str(len(poly)));lines.extend(' '.join(map(str,point)) for point in poly)
    lines += [f'{len(initial)} {len(ports)}',' '.join(map(str,initial))]
    lines.extend(f'{i} {j} {face} {int(vertical)}' for i,j,face,vertical in ports)
    path.write_text('\n'.join(lines)+'\n')


def run(name,p,q,subdivisions,initial,ports,binary,out,boundary_loads=()):
    source=out/f'{name}-{subdivisions}.txt';destination=source.with_suffix('.json');log=source.with_suffix('.log')
    write_input(p,q,subdivisions,initial,ports,source)
    command=[str(binary)]
    if boundary_loads:
        command.append('--boundary-loads')
        with source.open('a') as f:
            f.write(str(len(boundary_loads))+'\n')
            f.writelines(' '.join(map(str,row))+'\n' for row in boundary_loads)
    with source.open() as f,destination.open('w') as g:
        result=subprocess.run(command,stdin=f,stdout=g,stderr=subprocess.PIPE,text=True)
    log.write_text(result.stderr);assert result.returncode==0,(name,result.stderr)
    native=json.loads(destination.read_text())
    receipt=dict(name=name,subdivisions=subdivisions,steps=len(native['steps']),command=command,
                 **{k:str(v.relative_to(ROOT)) for k,v in [('input',source),('output',destination),('log',log)]},
                 **{k+'_sha256':sha(v) for k,v in [('input',source),('output',destination),('log',log)]})
    return native,receipt


def dense_step(z,c1,c2,g1,g2,drive,dt,incidence,series=None):
    n,m=incidence.shape;size=n+m
    def operator(c,g):
        return np.block([[np.zeros((n,n)),-incidence*g],[incidence.T/c,np.zeros((m,m))]])
    f1,f2=operator(c1,g1),operator(c2,g2);identity=np.eye(size)
    matrix=np.block([[identity-dt*A[0,0]*f1,-dt*A[0,1]*f2],[-dt*A[1,0]*f1,identity-dt*A[1,1]*f2]])
    source=np.r_[drive,np.zeros(m) if series is None else series]
    right=np.concatenate([z+dt*sum(row)*source for row in A])
    stage=np.linalg.solve(matrix,right).reshape(2,size)
    result=z+dt*(.5*(f1@stage[0]+f2@stage[1])+source)
    flow=.5*(g1*stage[0,n:]+g2*stage[1,n:])
    pressure=.5*(stage[0,:n]/c1+stage[1,:n]/c2)
    return result,flow,pressure


def audit_motion(name,p,q,event,subdivisions,binary,out,boundary_loads=(),case=None):
    native,receipt=run(name if case is None else case+"-"+name,p,q,subdivisions,INITIAL,PORTS,binary,out,boundary_loads)
    lines,distance=port_lines(p);incidence=np.zeros((5,4+len(boundary_loads)))
    distance=np.r_[distance,[row[2] for row in boundary_loads]]
    series=np.r_[np.zeros(4),[-row[3] for row in boundary_loads]]
    for k,(_,region,_,_) in enumerate(boundary_loads):incidence[region,4+k]=1
    selected=sorted({row[0] for row in boundary_loads})
    for edge,(i,j,_,_) in enumerate(PORTS):incidence[i,edge]=1;incidence[j,edge]=-1
    maxima=dict(volume_error_m3=0.,port_area_error_m2=0.,oracle_pressure_error_Pa=0.,oracle_momentum_error_Pa_s=0.,
                oracle_average_flow_error_m3_s=0.,oracle_average_pressure_error_Pa=0.,
                local_mass_error_m3=0.,global_mass_error_m3=0.,sealed_group_mass_error_m3=0.,closed_flow_m3_s=0.)
    snapshots=0
    def geometry(record):
        nonlocal snapshots
        pose=sampled(p,q,record['fraction'])
        if 'domain' in record:pose.update(domain=record['domain'],solids=record['solids'])
        checked=check_geometry(pose,record);regions,cells,error=exact_regions(pose,record)
        maxima['volume_error_m3']=max(maxima['volume_error_m3'],error)
        maxima['port_area_error_m2']=max(maxima['port_area_error_m2'],checked['maximum_port_area_error_m2'])
        v=np.array([r.area*p['depth'] for r in regions])
        area=np.array([regions[i].intersection(regions[j]).intersection(line).length*p['depth']
                       for (i,j,_,_),line in zip(PORTS,lines)])
        if boundary_loads:
            boundary,_=check_boundaries(pose,record,selected,regions)
            area=np.r_[area,[boundary.get((edge,region),(0.,None,None))[0] for edge,region,_,_ in boundary_loads]]
        snapshots+=1
        return pose,regions,cells,v,v/BULK,area/(RHO*distance)
    first=native['initial'];old_pose,old_regions,old_cells,old_v,c,g=geometry(first['geometry'])
    z=np.r_[c*INITIAL,np.zeros(incidence.shape[1])];old_native=first;mass0=sum(np.array(first['volume'])+first['charge'])
    density=[1+x/BULK for x in INITIAL];matches=0;sealed_mass=None;records=[];maximum_speed=0.;discharged=np.zeros(5)
    for state in native['steps']:
        _,_,_,_,c1,g1=geometry(state['stage1']);_,_,_,_,c2,g2=geometry(state['stage2'])
        pose,regions,cells,v,c,g=geometry(state['geometry']);dt=state['dt']
        match=correspondence(old_pose,pose,old_regions,old_cells,regions,cells)
        assert sorted(match['one_to_one'])==[[i,i] for i in range(5)]
        assert not any(match[k] for k in ['split_candidates','merge_candidates','old_without_witness','new_without_witness'])
        matches+=match['intermediate_witness_checks']
        z[5:][(g1==0)&(g2==0)]=0
        z,average_flow,average_pressure=dense_step(z,c1,c2,g1,g2,-(v-old_v)/dt,dt,incidence,series)
        z[5:][g==0]=0
        pressure=np.array(state['pressure']);momentum=np.array(state['momentum'])
        maxima['oracle_pressure_error_Pa']=max(maxima['oracle_pressure_error_Pa'],float(max(abs(pressure-z[:5]/c))))
        maxima['oracle_momentum_error_Pa_s']=max(maxima['oracle_momentum_error_Pa_s'],float(max(abs(momentum-z[5:]))))
        maxima['oracle_average_flow_error_m3_s']=max(maxima['oracle_average_flow_error_m3_s'],float(max(abs(np.array(state['average_flow'])-average_flow))))
        maxima['oracle_average_pressure_error_Pa']=max(maxima['oracle_average_pressure_error_Pa'],float(max(abs(np.array(state['average_pressure'])-average_pressure))))
        # Endpoint normal velocity under this coarse inertance model. This
        # reports model validity concerns; it is not a calibrated jet speed.
        opened=np.array(state['conductance'])>0
        maximum_speed=max(maximum_speed,float(max(abs(momentum[opened]/(RHO*distance[opened])),default=0.)))
        charge=np.array(state['volume'])+state['charge'];old_charge=np.array(old_native['volume'])+old_native['charge']
        old_discharged=discharged.copy()
        discharged+=dt*incidence[:,4:]@np.array(state['average_flow'][4:])
        balance=charge-old_charge+dt*incidence@np.array(state['average_flow'])
        maxima['local_mass_error_m3']=max(maxima['local_mass_error_m3'],float(max(abs(balance))))
        maxima['global_mass_error_m3']=max(maxima['global_mass_error_m3'],abs(sum(charge)+sum(discharged)-mass0))
        t=state['geometry']['fraction'];old_t=old_native['geometry']['fraction']
        closed=(old_t>event['fraction']+1e-10) if name=='closure' else (t<event['fraction']-1e-10)
        if closed:
            maxima['closed_flow_m3_s']=max(maxima['closed_flow_m3_s'],abs(state['average_flow'][1]))
            groups=np.array([sum((charge+discharged)[:2]),sum((charge+discharged)[2:])])
            if sealed_mass is None:sealed_mass=np.array([sum((old_charge+old_discharged)[:2]),sum((old_charge+old_discharged)[2:])])
            maxima['sealed_group_mass_error_m3']=max(maxima['sealed_group_mass_error_m3'],float(max(abs(groups-sealed_mass))))
        density.extend((1+pressure/BULK).tolist())
        records.append(dict(fraction=t,oracle_pressure_Pa=(z[:5]/c).tolist(),oracle_momentum_Pa_s=z[5:].tolist(),
                            oracle_average_flow_m3_s=average_flow.tolist(),oracle_average_pressure_Pa=average_pressure.tolist()))
        old_pose,old_regions,old_cells,old_v=pose,regions,cells,v;old_native=state
    assert max(abs(t-event['fraction']) for t in native['contact_cuts'])<1e-10
    assert maxima['oracle_pressure_error_Pa']<1e-3 and maxima['oracle_momentum_error_Pa_s']<1e-7,maxima
    assert maxima['oracle_average_flow_error_m3_s']<1e-9,maxima
    assert maxima['oracle_average_pressure_error_Pa']<1e-3,maxima
    assert maxima['local_mass_error_m3']<1e-17 and maxima['global_mass_error_m3']<1e-17 and maxima['sealed_group_mass_error_m3']<1e-17,maxima
    assert maxima['closed_flow_m3_s']==0,maxima
    final=native['steps'][-1];v0=np.array(first['volume']);vend=np.array(final['volume'])
    mass=np.array(final['volume'])+final['charge']
    audit_path=out/f'{receipt["name"]}-{subdivisions}-audit.json';audit_path.write_text(json.dumps(records,indent=2)+'\n')
    receipt.update(maxima=maxima,geometry_snapshots=snapshots,boundary_loads=boundary_loads,discharged_reference_volume_by_region_m3=discharged.tolist(),intermediate_witness_checks=matches,
                   group_balance_includes_external_ports=bool(boundary_loads),
                   relative_density_range=[min(density),max(density)],final_pressure_Pa=final['pressure'],final_flow_m3_s=final['flow'],
                   maximum_sampled_normal_speed_m_s=maximum_speed,maximum_sampled_normal_mach=maximum_speed/350,
                   volume_change_m3=(vend-v0).tolist(),transferred_mass_m3=(mass-v0*(1+INITIAL/BULK)).tolist(),
                   front_uniform_pressure_at_same_final_mass_Pa=BULK*(sum(mass[2:])/sum(vend[2:])-1),
                   audit=str(audit_path.relative_to(ROOT)),audit_sha256=sha(audit_path))
    print(json.dumps({k:receipt[k] for k in ['name','subdivisions','maxima','relative_density_range','final_pressure_Pa']}),flush=True)
    return receipt


def controls(binary,out):
    receipts=[]
    # Exact discrete standing wave, plus continuum-frequency spatial limit.
    for n in [8,16,32]:
        length=.12;height=.02;depth=.022;dt=.00037;dx=length/n
        p=dict(time=0.,depth=depth,domain=[[0,0],[length,0],[length,height],[0,height]],solids=[],
               x_edges=np.linspace(0,length,n+1).tolist(),y_edges=[0,height])
        q=dict(p,time=dt);initial=100*np.cos(np.pi*(np.arange(n)+.5)/n);ports=[(i,i+1,i+1,True) for i in range(n-1)]
        omega=2*350/dx*np.sin(np.pi/(2*n));previous=None
        for subdivisions in [8,16,32]:
            native,receipt=run(f'tube-{n}',p,q,subdivisions,initial,ports,binary,out)
            final=native['steps'][-1];expected=initial*np.cos(omega*dt)
            error=float(max(abs(np.array(final['pressure'])-expected)))
            initial_energy=native['initial']['energy'];energy_error=max(abs(s['energy']/initial_energy-1) for s in native['steps'])
            assert energy_error<1e-12 and all(max(abs(x) for x in s['wall_flow'])==0 for s in native['steps'])
            if previous is not None:assert previous/error>14,(n,subdivisions,previous,error)
            previous=error
            receipt.update(exact_discrete_pressure_error_Pa=error,maximum_relative_energy_error=energy_error,
                           discrete_frequency_Hz=omega/(2*np.pi),continuum_frequency_Hz=350/(2*length))
            receipts.append(receipt)
        assert previous<1e-4
    # Sealed moving piston: exact mass law; no branch momentum or transport.
    p=dict(time=0.,depth=.022,domain=[[0,0],[.04,0],[.04,.02],[0,.02]],solids=[],x_edges=[0,.05],y_edges=[0,.03])
    q=dict(p,time=.001,domain=[[0,0],[.04,0],[.04,.021],[0,.021]])
    native,receipt=run('sealed-expansion',p,q,16,[600.],[],binary,out)
    mass0=native['initial']['volume'][0]*(1+600/BULK)
    error=max(abs(s['pressure'][0]-BULK*(mass0/s['volume'][0]-1)) for s in native['steps'])
    assert error<1e-8
    receipt.update(exact_sealed_pressure_error_Pa=error);receipts.append(receipt)
    return receipts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_planar_inertive_motion_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/planar-inertive-motion')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_inertive_motion_reference.json')
    parser.add_argument('--sanitizer',action='store_true')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    control_receipts=controls(args.binary,args.output)
    contact_path=ROOT/'reports/contact_partition_reference.json';contacts=json.loads(contact_path.read_text())
    runs=[];refinement=[];inputs={}
    for name,start,contact in zip(['closure','release'],[5.031,5.062],contacts['events']):
        paths=[ROOT/f'out/planar-airspace/motion/pose-{t:.3f}.json' for t in [start,start+.001]]
        p,q=[json.loads(path.read_text()) for path in paths];event=contact['event']
        for pose in [p,q]:pose.update(x_edges=[pose['x_edges'][0],event['point_m'][0],pose['x_edges'][-1]],y_edges=[pose['y_edges'][0],.09,pose['y_edges'][-1]])
        previous=None;differences=[]
        for subdivisions in ([16] if args.sanitizer else [8,16,32,64,128]):
            result=audit_motion(name,p,q,event,subdivisions,args.binary,args.output);runs.append(result)
            if previous is not None:differences.append(dict(subdivisions=subdivisions,maximum_endpoint_difference_Pa=float(max(abs(np.array(result['final_pressure_Pa'])-previous)))))
            previous=np.array(result['final_pressure_Pa'])
        refinement.append(dict(name=name,comparisons=differences));inputs.update({str(path.relative_to(ROOT)):sha(path) for path in paths})
    sources=['tools/planar_inertive_motion_probe.cpp','tools/planar_inertive_motion_reference.py','include/ducklips/moving_tree.hpp','src/moving_tree.cpp',
             'tools/planar_boundary_reference.py','tools/planar_pressure_motion_reference.py','tools/planar_contact_partition_reference.py','tools/planar_region_motion_reference.py','tools/planar_regions_reference.py',
             'include/ducklips/planar_airspace.hpp','src/planar_airspace.cpp','include/ducklips/planar_face_motion.hpp','src/planar_face_motion.cpp']
    report=dict(use='inviscid canonical pressure/momentum control on five moving regions; two separately initialized intervals; no speech',
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],runs=runs,controls=control_receipts,refinement=refinement,input_poses=inputs,
                contact_receipt=str(contact_path.relative_to(ROOT)),contact_sha256=sha(contact_path),
                boundaries='sealed exterior',bulk_modulus_Pa=BULK,density_kg_m3=RHO,inverse_inertance='area / (density * grid-centre distance)',
                hashes={str(path.relative_to(ROOT)):sha(path) for path in [*(ROOT/s for s in sources),args.binary]},
                versions=dict(numpy=np.__version__,shapely=shapely.__version__,geos=shapely.geos_version_string))
    args.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(refinement,indent=2),flush=True)


if __name__=='__main__':main()
