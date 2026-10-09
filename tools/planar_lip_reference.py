#!/usr/bin/env python3
"""Connect trimmed native mouth regions to the assembled articulator lip tube.

An inviscid, open-lip integration reference, not a new speech recording.
"""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import shapely
from shapely.geometry import Polygon, LineString, box
from planar_inertive_motion_reference import ROOT, RHO, BULK, PORTS, INITIAL, write_input, dense_step, sha, run as run_without_lips
from planar_pressure_motion_reference import sampled, port_lines
from planar_regions_reference import check_geometry
from planar_region_motion_reference import exact_regions, correspondence, swept_boundaries
from planar_boundary_reference import check_boundaries
from planar_open_boundary_reference import clearance


def run(name,p,q,subdivisions,initial,ports,sections,attachments,binary,out):
    source=out/f'{name}-{subdivisions}.txt';destination=source.with_suffix('.json');log=source.with_suffix('.log')
    write_input(p,q,subdivisions,initial,ports,source)
    with source.open('a') as f:
        f.write(f'{len(sections)} {len(attachments)}\n')
        f.writelines(' '.join(map(str,row))+'\n' for row in sections)
        f.writelines(' '.join(map(str,row))+'\n' for row in attachments)
    command=[str(binary),'--lip-tube']
    with source.open() as f,destination.open('w') as g:
        result=subprocess.run(command,stdin=f,stdout=g,stderr=subprocess.PIPE,text=True)
    log.write_text(result.stderr);assert result.returncode==0,(name,result.stderr)
    native=json.loads(destination.read_text())
    receipt=dict(name=name,subdivisions=subdivisions,steps=len(native['steps']),command=command,sections=sections,attachments=attachments,
                 **{k:str(v.relative_to(ROOT)) for k,v in [('input',source),('output',destination),('log',log)]},
                 **{k+'_sha256':sha(v) for k,v in [('input',source),('output',destination),('log',log)]})
    return native,receipt


def check_partition(poses):
    maximum=0.;hashes={};volumes=[]
    for path in sorted(poses.glob('*.json')):
        p=json.loads(path.read_text());old_path=ROOT/'out/planar-boundary/poses'/path.name
        old=json.loads(old_path.read_text());partition=p['partition']
        # Reference-body fixture: both seams lie on straight centreline parts.
        assert partition==dict(oral_sections=64,epilarynx_sections=8,lip_sections=4,length_scale=1,reference_core_arc_cm=[2.1875,16.40625])
        mouth=p['boundary_edges']['mouth_end'];entry=p['boundary_edges']['supraglottal_entry']
        start=p['domain'][entry][1];end=p['domain'][mouth][0]
        cut=Polygon(old['domain']).intersection(box(-1,start,end,1))
        difference=Polygon(p['domain']).symmetric_difference(cut).area
        maximum=max(maximum,difference);assert difference<1e-17
        assert p['solids']==old['solids']
        core_length=(partition['reference_core_arc_cm'][1]-partition['reference_core_arc_cm'][0])*.01
        assert abs(core_length+(8+4)*.175/64-.175)<1e-16
        hashes.update({str(f.relative_to(ROOT)):sha(f) for f in [path,old_path]})
        volumes.append(dict(time=p['time'],core_volume_m3=sum(p['volume']),removed_volume_m3=sum(old['volume'])-sum(p['volume'])))
    assert len(volumes)==66
    return dict(poses=66,maximum_domain_symmetric_difference_m2=maximum,volumes=volumes,hashes=hashes)


def audit(name,p,q,sections,event,subdivisions,binary,out):
    n=5;count=len(sections);mouth=p['boundary_edges']['mouth_end']
    ports=PORTS+([(1,2,2,False)] if name=='lip-narrowing' else [])
    distance=p['domain'][mouth][0]-.5*(p['x_edges'][1]+p['x_edges'][2]);assert distance>0
    attachments=[(mouth,3,distance),(mouth,4,distance)]
    native,receipt=run(name,p,q,subdivisions,INITIAL,ports,sections,attachments,binary,out)
    edges=[(i,j) for i,j,_,_ in ports]+[(3,n),(4,n)]+[(n+i-1,n+i) for i in range(1,count)]+[(n+count-1,-1)]
    incidence=np.zeros((n+count,len(edges)))
    for k,(i,j) in enumerate(edges):
        incidence[i,k]=1
        if j>=0:incidence[j,k]=-1
    lines,core_distances=port_lines(p);sections=np.array(sections)
    if name=='lip-narrowing':lines.append(lines[0]);core_distances=np.r_[core_distances,core_distances[0]]
    maxima=dict(volume_m3=0.,internal_port_area_m2=0.,cap_area_m2=0.,oracle_pressure_Pa=0.,oracle_momentum_Pa_s=0.,
                average_flow_m3_s=0.,local_mass_m3=0.,global_mass_m3=0.,closed_neck_flow_m3_s=0.,metric_conductance=0.)
    snapshots=0;witnesses=0
    def metrics(record):
        nonlocal snapshots
        f=record['fraction'];pose=sampled(p,q,f)
        # std::lerp preserves equal endpoints exactly. Preserve the fixed seam
        # coordinates too: rounding a stationary cap outside its grid edge
        # by one ulp would erase its boundary intersection in the oracle.
        contours=[]
        for a,b,c in zip([p['domain'],*p['solids']],[q['domain'],*q['solids']],[pose['domain'],*pose['solids']]):
            a=np.array(a);b=np.array(b);contours.append(np.where(a==b,a,c).tolist())
        pose.update(domain=contours[0],solids=contours[1:])
        if 'domain' in record:pose.update(domain=record['domain'],solids=record['solids'])
        checked=check_geometry(pose,record);regions,cells,error=exact_regions(pose,record)
        boundaries,errors=check_boundaries(pose,record,[mouth],regions)
        maxima['volume_m3']=max(maxima['volume_m3'],error)
        maxima['internal_port_area_m2']=max(maxima['internal_port_area_m2'],checked['maximum_port_area_error_m2'])
        maxima['cap_area_m2']=max(maxima['cap_area_m2'],errors['area_m2'])
        a=(1-f)*sections[:,0]+f*sections[:,2];length=(1-f)*sections[:,1]+f*sections[:,3]
        v=np.r_[[r.area*p['depth'] for r in regions],a*length]
        core_area=np.array([regions[i].intersection(regions[j]).intersection(line).length*p['depth'] for (i,j,_,_),line in zip(ports,lines)])
        cap_areas=np.array([boundaries[(mouth,i)][0] for i in [3,4]])
        # Independent serial air-column inertances; cap shares allocate the
        # first half of ONE lip tube across the incoming region connections.
        branch_mass=RHO*(distance/cap_areas+.5*length[0]/(a[0]*cap_areas/sum(cap_areas)))
        lip_mass=RHO*.5*(length[:-1]/a[:-1]+length[1:]/a[1:])
        g=np.r_[core_area/(RHO*core_distances),1/branch_mass,1/lip_mass,2*a[-1]/(RHO*length[-1])]
        snapshots+=1
        return pose,regions,cells,v,v/BULK,g
    first=native['initial'];old_pose,old_regions,old_cells,old_v,c,g=metrics(first['geometry'])
    continuous_open=[]
    if name=='lip-narrowing':
        # Whole-face events are insufficient where two region pairs share a
        # grid face. A segment untouched by ANY swept polygon boundary is a
        # conservative, positive aperture witness for the complete interval.
        final_regions,_,_=exact_regions(q,native['steps'][-1]['geometry']);sweep=swept_boundaries(p,q)
        for (i,j,_,_),line in zip(ports,lines):
            opening=old_regions[i].intersection(old_regions[j]).intersection(final_regions[i]).intersection(final_regions[j]).intersection(line)
            area=opening.difference(sweep).length*p['depth'];assert area>0,(i,j,area)
            continuous_open.append(dict(pair=[i,j],minimum_certified_open_area_m2=area))
    z=np.r_[c*np.r_[INITIAL,np.zeros(count)],np.zeros(len(edges))]
    old_native=first;mass0=sum(old_v+z[:n+count]);outflow=0.;density=[1.];records=[]
    for s in native['steps']:
        _,_,_,_,c1,g1=metrics(s['stage1']);_,_,_,_,c2,g2=metrics(s['stage2'])
        pose,regions,cells,v,c,g=metrics(s['geometry']);dt=s['dt']
        match=correspondence(old_pose,pose,old_regions,old_cells,regions,cells)
        assert sorted(match['one_to_one'])==[[i,i] for i in range(n)]
        assert not any(match[k] for k in ['split_candidates','merge_candidates','old_without_witness','new_without_witness'])
        witnesses+=match['intermediate_witness_checks']
        z[n+count:][(g1==0)&(g2==0)]=0
        z,average_flow,_=dense_step(z,c1,c2,g1,g2,-(v-old_v)/dt,dt,incidence)
        z[n+count:][g==0]=0
        maxima['oracle_pressure_Pa']=max(maxima['oracle_pressure_Pa'],float(max(abs(np.array(s['pressure'])-z[:n+count]/c))))
        maxima['oracle_momentum_Pa_s']=max(maxima['oracle_momentum_Pa_s'],float(max(abs(np.array(s['momentum'])-z[n+count:]))))
        maxima['average_flow_m3_s']=max(maxima['average_flow_m3_s'],float(max(abs(np.array(s['average_flow'])-average_flow))))
        maxima['volume_m3']=max(maxima['volume_m3'],float(max(abs(np.array(s['volume'])-v))))
        maxima['metric_conductance']=max(maxima['metric_conductance'],float(max(abs(np.array(s['conductance'])-g))))
        mass=np.array(s['volume'])+s['charge'];old_mass=np.array(old_native['volume'])+old_native['charge']
        balance=mass-old_mass+dt*incidence@np.array(s['average_flow'])
        maxima['local_mass_m3']=max(maxima['local_mass_m3'],float(max(abs(balance))))
        outflow+=dt*s['average_flow'][-1]
        maxima['global_mass_m3']=max(maxima['global_mass_m3'],abs(sum(mass)+outflow-mass0))
        f=s['geometry']['fraction'];old_f=old_native['geometry']['fraction']
        if event is not None:
            closed=old_f>event['fraction']+1e-10 if name=='closure' else f<event['fraction']-1e-10
            if closed:maxima['closed_neck_flow_m3_s']=max(maxima['closed_neck_flow_m3_s'],abs(s['average_flow'][1]))
        density.extend((1+np.array(s['pressure'])/BULK).tolist())
        records.append(dict(fraction=f,oracle_pressure_Pa=(z[:n+count]/c).tolist()))
        old_pose,old_regions,old_cells,old_v=pose,regions,cells,v;old_native=s
    if event is not None:assert native['contact_cuts'] and max(abs(t-event['fraction']) for t in native['contact_cuts'])<1e-10
    assert maxima['oracle_pressure_Pa']<1e-3 and maxima['oracle_momentum_Pa_s']<1e-7,maxima
    assert maxima['average_flow_m3_s']<1e-9 and maxima['metric_conductance']<1e-12,maxima
    assert maxima['local_mass_m3']<1e-17 and maxima['global_mass_m3']<1e-17,maxima
    assert maxima['closed_neck_flow_m3_s']==0,maxima
    destination=out/f'{name}-{subdivisions}-audit.json';destination.write_text(json.dumps(records,indent=2)+'\n')
    receipt.update(maxima=maxima,geometry_snapshots=snapshots,intermediate_witness_checks=witnesses,continuously_open_port_witnesses=continuous_open,
                   relative_density_range=[min(density),max(density)],final_pressure_Pa=native['steps'][-1]['pressure'],
                   net_ambient_outflow_m3=outflow,lip_volume_change_m3=sum(native['steps'][-1]['volume'][n:])-sum(first['volume'][n:]),
                   audit=str(destination.relative_to(ROOT)),audit_sha256=sha(destination))
    print(json.dumps({k:receipt[k] for k in ['name','subdivisions','relative_density_range','maxima']}),flush=True)
    return receipt


def controls(binary,out):
    # Splitting a uniform cap and core across width must not duplicate the
    # first half-tube inertia. Symmetry has an exact unsplit counterpart.
    results=[]
    for moving in [False,True]:
        native=[]
        for split in [False,True]:
            p=dict(time=0.,depth=.022,domain=[[0,0],[.04,0],[.04,.02],[0,.02]],solids=[],x_edges=[0,.04],y_edges=[0,.01,.02] if split else [0,.02])
            q=dict(p,time=.00005);sections=[[.00044,.004,.0003 if moving else .00044,.006 if moving else .004]]*2
            initial=[100.,100.] if split else [100.];ports=[(0,1,1,False)] if split else []
            attachments=[(1,i,.02) for i in range(len(initial))]
            data,receipt=run(f'symmetry-{int(moving)}-{int(split)}',p,q,64,initial,ports,sections,attachments,binary,out)
            if not moving:
                receipt['maximum_relative_energy_error']=max(abs(s['energy']/data['initial']['energy']-1) for s in data['steps'])
                assert receipt['maximum_relative_energy_error']<1e-12
            native.append(data);results.append(receipt)
        worst=0.
        for a,b in zip(native[0]['steps'],native[1]['steps'],strict=True):
            worst=max(worst,max(abs(np.array(a['pressure'])-np.array(b['pressure'])[1:])),abs(b['pressure'][0]-b['pressure'][1]))
            assert abs(a['average_flow'][-1]-b['average_flow'][-1])<1e-15
        assert worst<1e-9
        results[-1]['split_unsplit_maximum_pressure_difference_Pa']=float(worst)
    # No nonzero floor is allowed to hide a fully closed lip section.
    source=out/'closed-lip.txt';write_input(p,q,4,[100.,100.],ports,source)
    with source.open('a') as f:f.write('1 2\n0 .004 0 .004\n1 0 .02\n1 1 .02\n')
    with source.open() as f:result=subprocess.run([str(binary),'--lip-tube'],stdin=f,capture_output=True,text=True)
    assert result.returncode==2 and 'closed-lip storage transport' in result.stderr and not result.stdout
    log=source.with_suffix('.log');log.write_text(result.stderr)
    results.append(dict(name='closed-lip-rejection',input=str(source.relative_to(ROOT)),input_sha256=sha(source),log=str(log.relative_to(ROOT)),log_sha256=sha(log),exit_code=result.returncode))
    return results


def geometry_controls(binary,out):
    geometry_binary=binary.parent/'ducklips_airspace_geometry_probe'
    trajectory=ROOT/'out/nato-airspace-geometry.csv';body=ROOT/'bodies/reference.body'
    directory=out/'geometry-checks';directory.mkdir(exist_ok=True)
    commands=[]
    for core in [False,True]:
        target=directory/('core' if core else 'legacy')
        command=[str(geometry_binary),str(trajectory),'5.031','5.032','512',str(target.with_suffix('.csv')),'--frames',str(target)]
        if core:command.extend(['--core',str(body)])
        result=subprocess.run(command,capture_output=True,text=True);assert result.returncode==0,result.stderr
        target.with_suffix('.log').write_text(result.stderr);commands.append(command)
        for path in target.glob('*.json'):
            previous=ROOT/('out/articulated-lip/poses' if core else 'out/planar-boundary/poses')/path.name
            assert json.loads(path.read_text())==json.loads(previous.read_text())
    scaled_body=directory/'scaled.body';scaled_body.write_text(body.read_text().replace('tract_length=0.175','tract_length=0.1925'))
    destination=directory/'scaled.json'
    command=[str(geometry_binary),str(trajectory),'5.031','5.031','512',str(directory/'scaled.csv'),str(destination),'--core',str(scaled_body)]
    result=subprocess.run(command,capture_output=True,text=True);assert result.returncode==0,result.stderr
    (directory/'scaled.log').write_text(result.stderr);commands.append(command)
    original=json.loads((directory/'core/pose-5.031.json').read_text());scaled=json.loads(destination.read_text());scale=scaled['partition']['length_scale'];assert scale==1.1
    errors={}
    for key,power in [('domain',1),('solids',1),('x_edges',1),('y_edges',1),('volume',3),('vertical_face',2),('horizontal_face',2)]:
        if key=='solids':error=max(float(np.max(abs(np.array(b)-np.array(a)*scale))) for a,b in zip(original[key],scaled[key]))
        else:error=float(np.max(abs(np.array(scaled[key])-np.array(original[key])*scale**power)))
        errors[key]=error;assert error<(1e-16 if power==1 else 1e-17),(key,error)
    return dict(commands=commands,length_scale=scale,maximum_scaling_errors=errors,
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [geometry_binary,trajectory,body,*[p for p in directory.rglob('*') if p.is_file()]]})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_planar_inertive_motion_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/articulated-lip')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_lip_reference.json')
    parser.add_argument('--sanitizer',action='store_true');args=parser.parse_args()
    args.binary=args.binary.resolve();args.output=args.output.resolve();args.report=args.report.resolve();args.output.mkdir(parents=True,exist_ok=True)
    poses=ROOT/'out/articulated-lip/poses';partition=check_partition(poses)
    airway_path=ROOT/'out/nato-spatial-airway.csv'
    airway={round(float(r['time'])*1000):r for r in csv.DictReader(airway_path.open())}
    contact_path=ROOT/'reports/contact_partition_reference.json';contacts=json.loads(contact_path.read_text())['events']
    runs=[];refinements=[];clearances={};used_inputs={};fine_runs=[];regressions=[]
    for name,start,contact in [('closure',5.031,contacts[0]['event']),('release',5.062,contacts[1]['event']),('lip-narrowing',5.084,None)]:
        paths=[poses/f'pose-{t:.3f}.json' for t in [start,start+.001]]
        p,q=[json.loads(path.read_text()) for path in paths]
        # End the coarse grid at the true seam, keeping its centre behind it.
        cut=(contact or contacts[1]['event'])['point_m'][0]
        for pose in [p,q]:
            mouth=pose['boundary_edges']['mouth_end'];entry=pose['boundary_edges']['supraglottal_entry']
            pose.update(x_edges=[-.01,cut,pose['domain'][mouth][0]],y_edges=[pose['domain'][entry][1],.09,.105])
        clearances[name]=clearance(p,q)
        rows=[airway[round(t*1000)] for t in [start,start+.001]];layout=p['partition'];first=layout['oral_sections']-layout['lip_sections']
        sections=[[float(r[k]) for r in rows for k in [f'area_{i}',f'length_{i}']] for i in range(first,layout['oral_sections'])]
        previous=None;differences=[]
        for subdivisions in ([16] if args.sanitizer else [16,32,64,128]):
            receipt=audit(name,p,q,sections,contact,subdivisions,args.binary,args.output);runs.append(receipt)
            pressure=np.array(receipt['final_pressure_Pa'])
            if previous is not None:differences.append(dict(subdivisions=subdivisions,maximum_endpoint_difference_Pa=float(max(abs(pressure-previous)))))
            previous=pressure
        if not args.sanitizer:
            # Short lip sections add fast undamped modes. Check temporal
            # convergence beyond the full geometry/oracle audit above; these
            # additional runs are explicitly native-only refinement evidence.
            for subdivisions in ([256,512,1024,2048] if name=='lip-narrowing' else [256,512,1024]):
                data,fine=run(name,p,q,subdivisions,INITIAL,PORTS+([(1,2,2,False)] if name=='lip-narrowing' else []),
                              sections,receipt['attachments'],args.binary,args.output)
                pressure=np.array(data['steps'][-1]['pressure']);difference=float(max(abs(pressure-previous)))
                fine.update(final_pressure_Pa=pressure.tolist(),maximum_endpoint_difference_Pa=difference,
                            validation='native-only temporal refinement; full independent geometry/state audit at 16..128 subdivisions')
                fine_runs.append(fine);differences.append(dict(subdivisions=subdivisions,maximum_endpoint_difference_Pa=difference));previous=pressure
                print(json.dumps(dict(name=name,subdivisions=subdivisions,native_refinement_difference_Pa=difference)),flush=True)
            assert differences[-1]['maximum_endpoint_difference_Pa']<.02,differences
        refinements.append(dict(name=name,comparisons=differences))
        used_inputs.update({str(path.relative_to(ROOT)):sha(path) for path in paths})
        if name in ['closure','release']:
            old_p,old_q=[json.loads((ROOT/'out/planar-boundary/poses'/path.name).read_text()) for path in paths]
            for pose in [old_p,old_q]:pose.update(x_edges=[-.01,contact['point_m'][0],.1],y_edges=[0,.09,.105])
            data,regression=run_without_lips('original-'+name,old_p,old_q,128,INITIAL,PORTS,args.binary,args.output)
            old_path=ROOT/f'out/planar-inertive-motion/{name}-128.json';old=json.loads(old_path.read_text())
            worst=max(abs(a-b) for new,prior in zip(data['steps'],old['steps'],strict=True)
                      for key in ['charge','momentum','pressure','flow','average_flow','volume','compliance','conductance']
                      for a,b in zip(new[key],prior[key],strict=True))
            assert worst==0
            regression.update(previous_output=str(old_path.relative_to(ROOT)),previous_output_sha256=sha(old_path),maximum_state_difference=worst)
            regressions.append(regression)
    analytic=controls(args.binary,args.output)
    geometric=geometry_controls(args.binary,args.output)
    sources=['tools/planar_lip_reference.py','tools/planar_inertive_motion_probe.cpp','tools/airspace_geometry_probe.cpp',
             'tools/planar_inertive_motion_reference.py','tools/planar_regions_reference.py','tools/planar_region_motion_reference.py',
             'tools/planar_boundary_reference.py','tools/planar_pressure_motion_reference.py','tools/planar_contact_partition_reference.py',
             'tools/planar_open_boundary_reference.py','src/planar_airspace.cpp','src/planar_face_motion.cpp','src/moving_graph.cpp',
             'include/ducklips/moving_graph.hpp','src/articulation.cpp','include/ducklips/articulation.hpp','bodies/reference.body']
    report=dict(use='trimmed planar core with native assembled open-lip sections; inviscid ideal outlet, no new speech',
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],partition=partition,runs=runs,controls=analytic,refinement=refinements,
                fine_native_runs=fine_runs,unchanged_tree_controls=regressions,geometry_controls=geometric,
                cap_clearance_m=clearances,input_poses=used_inputs,
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),airway_path,contact_path,args.binary]},
                versions=dict(numpy=np.__version__,shapely=shapely.__version__,geos=shapely.geos_version_string))
    args.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(refinements,indent=2))


if __name__=='__main__':main()
