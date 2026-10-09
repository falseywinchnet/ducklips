#!/usr/bin/env python3
"""Actual Bravo mouth, integrated lips and epilarynx with an implicit inlet.

Held aperture controls: no moving folds, losses, radiation or new speech.
"""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import brentq
import shapely
from shapely.geometry import LineString
from planar_lip_profile_reference import setup,pose_at,PORTS
from lip_profile_reference import oracle as lip_metrics,vector
from planar_inertive_motion_reference import ROOT,RHO,BULK,write_input,dense_step,sha
from planar_regions_reference import check_geometry
from planar_region_motion_reference import exact_regions,correspondence,swept_boundaries
from planar_boundary_reference import check_boundaries
from planar_open_boundary_reference import clearance

EDGES=[(0,1),(0,2),(1,2),(2,3),(3,4),(4,5),(5,-1),(0,6),*[(i,i+1) for i in range(6,13)]]
B=np.zeros((14,len(EDGES)))
for k,(i,j) in enumerate(EDGES):
    B[i,k]=1
    if j>=0:B[j,k]=-1


def run(name,p,q,profiles,epi,subdivisions,pressure,area,reservoir,previous,binary,out,compact):
    source=out/f'{name}-{subdivisions}.txt';destination=source.with_suffix('.json');log=source.with_suffix('.log')
    write_input(p,q,subdivisions,[pressure]*2,PORTS,source)
    mouth=p['boundary_edges']['mouth_end'];entry=p['boundary_edges']['supraglottal_entry']
    mouth_distance=.5*(p['x_edges'][-1]-p['x_edges'][0]);entry_distance=.5*(p['y_edges'][1]-p['y_edges'][0])
    with source.open('a') as f:
        f.write('4 2\n');f.writelines(' '.join(map(str,v))+'\n' for v in profiles)
        f.write(' '.join(map(str,[pressure]*4))+'\n')
        f.writelines(f'{mouth} {i} {mouth_distance}\n' for i in range(2))
        f.write('8 1\n');f.writelines(' '.join(map(str,v))+'\n' for v in np.column_stack([epi[0,0],epi[0,1],epi[1,0],epi[1,1]]))
        f.write(' '.join(map(str,[pressure]*8))+'\n')
        f.write(f'{entry} 0 {entry_distance}\n{area} {reservoir}\n')
        if previous is not None:
            f.writelines(' '.join(map(str,previous[k]))+'\n' for k in ['charge','momentum'])
            f.write(str(previous['source_momentum'])+'\n')
    command=[str(binary),'--lip-profile','--epilarynx','--glottal-orifice']
    if previous is not None:command.append('--initial-state')
    if compact:command.append('--compact-geometry')
    with source.open() as f,destination.open('w') as g:r=subprocess.run(command,stdin=f,stdout=g,stderr=subprocess.PIPE,text=True)
    log.write_text(r.stderr);assert r.returncode==0,(name,r.stderr)
    native=json.loads(destination.read_text())
    if previous is not None:
        for key in ['volume','compliance','conductance','charge','momentum','source_momentum']:
            assert native['initial'][key]==previous[key],('continuation',key)
    receipt=dict(name=name,subdivisions=subdivisions,command=command,steps=len(native['steps']),aperture_m2=area,reservoir_Pa=reservoir,
                 continued=previous is not None,
                 **{k:str(v.relative_to(ROOT)) for k,v in [('input',source),('output',destination),('log',log)]},
                 **{k+'_sha256':sha(v) for k,v in [('input',source),('output',destination),('log',log)]})
    return native,receipt


def audit(native,p,q,profiles,epi,pressure,area,reservoir,previous,full):
    n=14;mouth=p['boundary_edges']['mouth_end'];entry=p['boundary_edges']['supraglottal_entry']
    md=.5*(p['x_edges'][-1]-p['x_edges'][0]);ed=.5*(p['y_edges'][1]-p['y_edges'][0]);cd=.5*(p['y_edges'][-1]-p['y_edges'][0])
    line=LineString([(p['x_edges'][0],.09),(p['x_edges'][-1],.09)])
    maxima={k:0. for k in ['volume_m3','cap_area_m2','internal_area_m2','conductance','pressure_Pa','momentum_Pa_s','average_flow_m3_s',
                            'source_flow_m3_s','source_pressure_Pa','source_momentum_Pa_s','orifice_law_Pa','local_mass_m3','global_mass_m3']}
    checks=0;witnesses=0
    def geometry(record):
        nonlocal checks
        fraction=record['fraction'];pose=pose_at(p,q,record)
        checked=check_geometry(pose,record);regions,cells,error=exact_regions(pose,record)
        boundary,errors=check_boundaries(pose,record,[mouth,entry],regions)
        maxima['volume_m3']=max(maxima['volume_m3'],error)
        maxima['cap_area_m2']=max(maxima['cap_area_m2'],errors['area_m2'])
        maxima['internal_area_m2']=max(maxima['internal_area_m2'],checked['maximum_port_area_error_m2'])
        lip=lip_metrics((1-fraction)*profiles[0]+fraction*profiles[1]);ea,el=(1-fraction)*epi[0]+fraction*epi[1]
        v=np.r_[[r.area*p['depth'] for r in regions],lip['volume'],ea*el]
        cap=np.array([boundary[(mouth,i)][0] for i in range(2)]);share=cap/sum(cap)
        mouth_g=np.zeros(2) if lip['conductance'][0]==0 else 1/(RHO*md/cap+1/(lip['conductance'][0]*share))
        entry_area=boundary[(entry,0)][0]
        assert set(boundary)=={(mouth,0),(mouth,1),(entry,0)}
        entry_g=1/(RHO*(ed/entry_area+.5*el[0]/ea[0]))
        internal=regions[0].intersection(regions[1]).intersection(line).length*p['depth']
        g=np.r_[internal/(RHO*cd),mouth_g,lip['conductance'][1:],entry_g,1/(RHO*(.5*el[:-1]/ea[:-1]+.5*el[1:]/ea[1:]))]
        checks+=1
        return pose,regions,cells,v,v/BULK,g
    old=native['initial'];mass0=sum(np.array(old['volume'])+old['charge']);inlet=0.;outlet=0.;z=None;lm=old['source_momentum']
    continuous_open=None
    if full:
        old_pose,old_regions,old_cells,old_v,c,g=geometry(old['geometry'])
        z=np.r_[c*pressure,np.zeros(len(EDGES))] if previous is None else previous[0].copy()
        lm=0. if previous is None else previous[1]
        last_regions,_,_=exact_regions(q,native['steps'][-1]['geometry'])
        opening=old_regions[0].intersection(old_regions[1]).intersection(last_regions[0]).intersection(last_regions[1]).intersection(line)
        continuous_open=opening.difference(swept_boundaries(p,q)).length*p['depth'];assert continuous_open>0
    p0,p1=profiles;cuts=[]
    if p0[-1]==p1[-1]==0 and (p0[2]<1)!=(p1[2]<1):
        t=(1-p0[2])/(p1[2]-p0[2])
        if 0<t<1:cuts.append(t)
    assert native['lip_contact_cuts']==cuts and not native['contact_cuts']
    for state in native['steps']:
        dt=state['dt'];f=.5*(old['geometry']['fraction']+state['geometry']['fraction'])
        ea,el=(1-f)*epi[0]+f*epi[1];inertance=RHO*.5*el[-1]/ea[-1]
        assert abs(state['source_inertance']-inertance)<1e-12
        if full:
            _,_,_,_,c1,g1=geometry(state['stage1']);_,_,_,_,c2,g2=geometry(state['stage2'])
            pose,regions,cells,v,c,g=geometry(state['geometry'])
            match=correspondence(old_pose,pose,old_regions,old_cells,regions,cells)
            assert match['one_to_one']==[[0,0],[1,1]] and not any(match[k] for k in ['split_candidates','merge_candidates','old_without_witness','new_without_witness'])
            witnesses+=match['intermediate_witness_checks'];z[n:][(g1==0)&(g2==0)]=0
            free,free_flow,free_pressure=dense_step(z,c1,c2,g1,g2,-(v-old_v)/dt,dt,B)
            unit=np.zeros(n);unit[-1]=1
            response,response_flow,response_pressure=dense_step(np.zeros_like(z),c1,c2,g1,g2,unit,dt,B)
            effective=reservoir-free_pressure[-1]+2*lm/dt
            impedance=response_pressure[-1]+2*inertance/dt
            if area==0:u=0.
            else:
                # Independent bracketed root of the source pressure law.
                k=RHO/(2*area**2);bound=np.sqrt(abs(effective)/k)+1e-12
                u=brentq(lambda u:effective-impedance*u-k*u*abs(u),-bound,bound,xtol=1e-22,rtol=1e-14)
            lm=2*inertance*u-lm if area>0 else 0.
            z=free+u*response;z[n:][g==0]=0
            maxima['pressure_Pa']=max(maxima['pressure_Pa'],float(max(abs(np.array(state['pressure'])-z[:n]/c))))
            maxima['momentum_Pa_s']=max(maxima['momentum_Pa_s'],float(max(abs(np.array(state['momentum'])-z[n:]))))
            maxima['average_flow_m3_s']=max(maxima['average_flow_m3_s'],float(max(abs(np.array(state['average_flow'])-free_flow-u*response_flow))))
            maxima['source_flow_m3_s']=max(maxima['source_flow_m3_s'],abs(state['source_flow']-u))
            maxima['source_pressure_Pa']=max(maxima['source_pressure_Pa'],abs(state['source_pressure']-free_pressure[-1]-u*response_pressure[-1]))
            maxima['source_momentum_Pa_s']=max(maxima['source_momentum_Pa_s'],abs(state['source_momentum']-lm))
            maxima['volume_m3']=max(maxima['volume_m3'],float(max(abs(np.array(state['volume'])-v))))
            maxima['conductance']=max(maxima['conductance'],float(max(abs(np.array(state['conductance'])-g))))
            old_pose,old_regions,old_cells,old_v=pose,regions,cells,v
        flow=state['source_flow'];load=state['source_pressure']
        if area>0:
            residual=reservoir-load-RHO*flow*abs(flow)/(2*area**2)-(state['source_momentum']-old['source_momentum'])/dt
            maxima['orifice_law_Pa']=max(maxima['orifice_law_Pa'],abs(residual))
        else:assert flow==state['source_momentum']==0
        current=np.array(state['volume'])+state['charge'];prior=np.array(old['volume'])+old['charge']
        balance=current-prior+dt*B@np.array(state['average_flow']);balance[-1]-=dt*flow
        maxima['local_mass_m3']=max(maxima['local_mass_m3'],float(max(abs(balance))))
        inlet+=dt*flow;outlet+=dt*state['average_flow'][6]
        maxima['global_mass_m3']=max(maxima['global_mass_m3'],abs(sum(current)+outlet-inlet-mass0))
        if state['conductance'][6]==old['conductance'][6]==0:assert state['average_flow'][6]==0
        old=state
    assert maxima['local_mass_m3']<1e-17 and maxima['global_mass_m3']<1e-17 and maxima['orifice_law_Pa']<1e-6,maxima
    if full:assert maxima['pressure_Pa']<1e-3 and maxima['source_flow_m3_s']<1e-10 and maxima['conductance']<1e-12,maxima
    return dict(independent_audit=full,maxima=maxima,geometry_checks=checks,intermediate_witness_checks=witnesses,
                continuously_open_internal_witness_m2=continuous_open,net_inlet_volume_m3=inlet,net_outlet_volume_m3=outlet,
                final_pressure_Pa=old['pressure'],final_source_flow_m3_s=old['source_flow']), (z,lm)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_planar_inertive_motion_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/planar-source')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_source_reference.json')
    parser.add_argument('--subdivisions',type=int,nargs='+',default=[16,32,64,128,256,512])
    args=parser.parse_args();args.binary=args.binary.resolve();args.output=args.output.resolve();args.report=args.report.resolve();args.output.mkdir(parents=True,exist_ok=True)
    lip_rows=list(csv.DictReader((ROOT/'out/lip-profile/nato.csv').open()))
    airway=list(csv.DictReader((ROOT/'out/nato-spatial-airway.csv').open()));runs=[];inputs={};refinements=[]
    for name,start,pressure,aperture,reservoir in [('closing',1284,0.,3e-6,700.),('opening',1374,700.,3e-6,700.)]:
        last=None
        for subdivisions in args.subdivisions:
            previous=None;previous_oracle=None;mass0=None;inlet=0;outlet=0
            for interval in range(2):
                index=start+interval;p,q,paths=setup(index/1000)
                profiles=[vector(lip_rows[i]) for i in [index,index+1]]
                epi=np.array([[[float(airway[t][f'{field}_{i}']) for i in reversed(range(8))] for field in ['area','length']] for t in [index,index+1]])
                assert sum(epi[0,1])==p['domain'][p['boundary_edges']['supraglottal_entry']][1]
                native,receipt=run(f'{name}-{interval}',p,q,profiles,epi,subdivisions,pressure,aperture,reservoir,previous,args.binary,args.output,subdivisions>64)
                checked,previous_oracle=audit(native,p,q,profiles,epi,pressure,aperture,reservoir,previous_oracle,subdivisions<=64)
                receipt.update(checked);receipt['cap_clearance_m']=clearance(p,q);runs.append(receipt)
                if mass0 is None:mass0=sum(np.array(native['initial']['volume'])+native['initial']['charge'])
                inlet+=checked['net_inlet_volume_m3'];outlet+=checked['net_outlet_volume_m3']
                previous=native['steps'][-1];inputs.update({str(path.relative_to(ROOT)):sha(path) for path in paths})
            mass_error=abs(sum(np.array(previous['volume'])+previous['charge'])+outlet-inlet-mass0)
            assert mass_error<1e-17
            final=np.array(previous['pressure']);change=None if last is None else float(max(abs(final-last)))
            refinements.append(dict(name=name,subdivisions=subdivisions,maximum_endpoint_change_Pa=change,whole_control_mass_error_m3=mass_error))
            last=final
            print(json.dumps(dict(name=name,subdivisions=subdivisions,change_Pa=change,source_flow=previous['source_flow'],maxima=runs[-1]['maxima'])),flush=True)
    sources=['src/graph_orifice.cpp','include/ducklips/graph_orifice.hpp','tools/graph_orifice_probe.cpp','src/moving_graph.cpp','include/ducklips/moving_graph.hpp',
             'tools/planar_inertive_motion_probe.cpp','tools/planar_source_reference.py','tools/planar_lip_profile_reference.py','tools/lip_profile_reference.py',
             'tools/planar_inertive_motion_reference.py','tools/planar_regions_reference.py','tools/planar_region_motion_reference.py','tools/planar_pressure_motion_reference.py',
             'tools/planar_boundary_reference.py','tools/planar_open_boundary_reference.py','src/planar_airspace.cpp','src/planar_face_motion.cpp','src/lip_profile.cpp',
             'src/sources.cpp','bodies/reference.body','out/lip-profile/nato.csv','out/nato-spatial-airway.csv']
    args.report.write_text(json.dumps(dict(use='held glottal aperture plus epilarynx connected to actual moving mouth/lips; no new speech',
        command=[sys.executable,*sys.argv],runs=runs,refinements=refinements,input_poses=inputs,
        hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),args.binary]},
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,shapely=shapely.__version__,geos=shapely.geos_version_string)),indent=2)+'\n')


if __name__=='__main__':main()
