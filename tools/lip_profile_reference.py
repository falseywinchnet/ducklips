#!/usr/bin/env python3
"""Independent axial lip geometry and closed-port storage audit. No audio."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import scipy
from scipy.integrate import quad
from planar_inertive_motion_reference import dense_step

ROOT=Path(__file__).resolve().parents[1]
FIELDS=['length','base_area','activation','place','core','width','target_area']
RHO=1.14;BULK=RHO*350**2
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def vector(row):return np.array([float(row[k]) for k in FIELDS])


def oracle(p,n=4):
    length,base,activation,place,core,width,target=p
    activation=np.clip(activation,0,1);target=min(base,target)
    minimum=target+(base-target)*(1-activation)**2*(1+2*activation)
    tail=np.sqrt(-np.log(1e-4))
    def area(x):
        u=max(abs(x-place)-core,0)/width
        if u>tail:return base
        return minimum+(base-minimum)*(-np.expm1(-u*u))
    def integrate(a,b,kind):
        cuts=sorted(set([a,b,*[x for x in [place-core,place+core,place-core-tail*width,place+core+tail*width] if a<x<b]]))
        if kind==2 and minimum==0 and max(a,place-core)<=min(b,place+core):return np.inf
        def f(x):
            v=area(x)
            return 1/v if kind==2 else v*(x if kind==1 else 1)
        return sum(quad(f,x,y,epsabs=0,epsrel=2e-11,limit=150)[0] for x,y in zip(cuts[:-1],cuts[1:]))
    v=[];centre=[];mid=[]
    for i in range(n):
        a,b=length*(n-i-1)/n,length*(n-i)/n
        volume=integrate(a,b,0);v.append(volume);centre.append(integrate(a,b,1)/volume if volume>0 else (a+b)/2);mid.append(area((a+b)/2))
    paths=[(centre[0],length),*[(centre[i],centre[i-1]) for i in range(1,n)],(0,centre[-1])]
    g=np.array([1/(RHO*integrate(a,b,2)) for a,b in paths]);v=np.array(v)
    return dict(volume=v,centroid=np.array(centre),compliance=v/BULK,conductance=g,midpoint_area=np.array(mid))


def invoke(binary,mode,source,destination):
    with source.open() as f,destination.open('w') as g:result=subprocess.run([str(binary),mode],stdin=f,stdout=g,stderr=subprocess.PIPE,text=True)
    log=source.with_suffix('.log');log.write_text(result.stderr);assert result.returncode==0,(source,result.stderr)
    return dict(command=[str(binary),mode],**{k:str(v.relative_to(ROOT)) for k,v in [('input',source),('output',destination),('log',log)]},
                **{k+'_sha256':sha(v) for k,v in [('input',source),('output',destination),('log',log)]})


def geometry(rows,binary,out,sanitizer):
    airway=ROOT/'out/lip-profile/airway.csv';previous=ROOT/'out/nato-spatial-airway.csv'
    assert sha(airway)==sha(previous),'exporting the profile changed assembled tract geometry'
    assert sha(ROOT/'out/lip-profile/geometry.csv')==sha(ROOT/'out/nato-airspace-geometry.csv'),'profile replay changed tongue state'
    old=list(csv.DictReader(airway.open()));assert len(rows)==len(old)==23239
    selected=list(range(len(rows))) if not sanitizer else [0,1283,1284,1285,1286,1374,1375,1376,13585,13672,13749,13802]
    source=out/'nato-metrics.txt';destination=source.with_suffix('.jsonl')
    source.write_text('\n'.join('4 '+' '.join(map(str,vector(rows[i]))) for i in selected)+'\n')
    receipt=invoke(binary,'--metrics',source,destination);records=[json.loads(s) for s in destination.read_text().splitlines()]
    midpoint_error=0.;min_closed_volume=np.inf;closed=[];worst={k:0. for k in ['volume','centroid','conductance']}
    for i,actual in zip(selected,records,strict=True):
        midpoint_error=max(midpoint_error,max(abs(np.array(actual['midpoint_area'])-[float(old[i][f'area_{j}']) for j in range(60,64)])))
        assert min(actual['volume'])>0
        if float(old[i]['area_63'])==0:
            assert actual['conductance'][-1]==0 and actual['conductance'][-2]>0
            min_closed_volume=min(min_closed_volume,actual['volume'][-1]);closed.append(i)
        if i%401==0 or i in [1283,1284,1285,1286,1374,1375,1376,13585,13672,13749,13802]:
            reference=oracle(vector(rows[i]))
            for k in worst:worst[k]=max(worst[k],float(max(abs(np.array(actual[k])-reference[k]))))
    assert midpoint_error<1e-18 and closed and worst['volume']<1e-17 and worst['centroid']<1e-12 and worst['conductance']<1e-10,(midpoint_error,worst)
    # Geometric storage is independent of how the axial interval is partitioned,
    # including a finer discretization with wholly empty terminal cells.
    p=vector(rows[1285]);source2=out/'partition.txt';destination2=source2.with_suffix('.jsonl')
    source2.write_text('\n'.join(str(n)+' '+' '.join(map(str,p)) for n in [1,2,4,8,16,32])+'\n')
    partition=invoke(binary,'--metrics',source2,destination2);parts=[json.loads(s) for s in destination2.read_text().splitlines()]
    volumes=[sum(s['volume']) for s in parts];assert max(volumes)-min(volumes)<1e-19
    assert parts[2]['volume'][-1]>0 and parts[3]['volume'][-1]==0
    receipt.update(poses=len(selected),maximum_midpoint_area_difference_m2=midpoint_error,closed_frames=len(closed),
                   minimum_last_section_air_volume_at_closed_midpoint_m3=min_closed_volume,
                   independent_maxima=worst,partition_control=partition,partition_volume_range_m3=[min(volumes),max(volumes)])
    return receipt


def sequence(name,profiles,durations,subdivisions,initial,inlet,binary,out,audit=True):
    n=4;source=out/f'{name}-{subdivisions}.txt';destination=source.with_suffix('.jsonl')
    lines=[f'{n} {len(durations)} {subdivisions} {initial} {inlet}',' '.join(map(str,profiles[0]))]
    for duration,p in zip(durations,profiles[1:],strict=True):lines.extend([str(duration),' '.join(map(str,p))])
    source.write_text('\n'.join(lines)+'\n');receipt=invoke(binary,'--sequence',source,destination)
    states=[json.loads(s) for s in destination.read_text().splitlines()]
    contacts=[]
    for k,(p,q) in enumerate(zip(profiles[:-1],profiles[1:])):
        if p[-1]==q[-1]==0 and (p[2]<1)!=(q[2]<1):
            fraction=(1-p[2])/(q[2]-p[2])
            if 0<fraction<1:contacts.append((k,fraction))
    assert len(states)>=1+len(durations)*subdivisions
    for k,f in contacts:assert any(s.get('interval')==k and s.get('fraction')==f for s in states)
    first=states[0];old=first;z=np.r_[first['charge'],first['momentum']]
    if audit:
        old_oracle_volume=oracle(profiles[0])['volume'];z=np.r_[old_oracle_volume/BULK*initial,np.zeros(n+1)]
    b=np.zeros((n,n+1));b[0,0]=1;b[-1,-1]=1
    for k in range(1,n):b[k-1,k]=1;b[k,k]=-1
    series=np.r_[-inlet,np.zeros(n)];mass0=sum(np.array(first['geometry']['volume'])+first['charge']);external=0
    maxima=dict(pressure_Pa=0.,momentum_Pa_s=0.,average_flow_m3_s=0.,local_mass_m3=0.,global_mass_m3=0.,sealed_flow_m3_s=0.)
    metric_errors=dict(volume=0.,centroid=0.,conductance=0.);density=[];profile_checks=0
    for state in states[1:]:
        k=state['interval'];f=state['fraction'];a=state['start_fraction'];dt=state['dt']
        if audit:
            geometries=[]
            for key,t in [('stage1',a+.21132486540518711775*(f-a)),('stage2',a+.78867513459481288225*(f-a)),('geometry',f)]:
                reference=oracle((1-t)*profiles[k]+t*profiles[k+1]);geometries.append(reference);profile_checks+=1
                for metric in metric_errors:metric_errors[metric]=max(metric_errors[metric],float(max(abs(np.array(state[key][metric])-reference[metric]))))
            g1,g2,end=geometries;drive=-(end['volume']-old_oracle_volume)/dt
            z[n:][(g1['conductance']==0)&(g2['conductance']==0)]=0
            z,flow,_=dense_step(z,g1['compliance'],g2['compliance'],g1['conductance'],g2['conductance'],drive,dt,b,series)
            z[n:][end['conductance']==0]=0
            maxima['pressure_Pa']=max(maxima['pressure_Pa'],float(max(abs(np.array(state['pressure'])-z[:n]/end['compliance']))))
            maxima['momentum_Pa_s']=max(maxima['momentum_Pa_s'],float(max(abs(np.array(state['momentum'])-z[n:]))))
            maxima['average_flow_m3_s']=max(maxima['average_flow_m3_s'],float(max(abs(np.array(state['average_flow'])-flow))))
            old_oracle_volume=end['volume']
        current=np.array(state['geometry']['volume'])+state['charge'];prior=np.array(old['geometry']['volume'])+old['charge']
        balance=current-prior+dt*b@np.array(state['average_flow']);external+=dt*(state['average_flow'][0]+state['average_flow'][-1])
        maxima['local_mass_m3']=max(maxima['local_mass_m3'],float(max(abs(balance))))
        maxima['global_mass_m3']=max(maxima['global_mass_m3'],abs(sum(current)+external-mass0))
        if state['stage1']['conductance'][-1]==state['stage2']['conductance'][-1]==0:
            maxima['sealed_flow_m3_s']=max(maxima['sealed_flow_m3_s'],abs(state['average_flow'][-1]))
        density.extend((1+np.array(state['pressure'])/BULK).tolist());old=state
    assert maxima['local_mass_m3']<1e-18 and maxima['global_mass_m3']<1e-18 and maxima['sealed_flow_m3_s']==0,maxima
    if audit:assert maxima['pressure_Pa']<1e-3 and maxima['average_flow_m3_s']<1e-9,maxima
    receipt.update(independent_audit=audit,independent_geometry_checks=profile_checks,contact_cuts=contacts,maxima=maxima,metric_errors=metric_errors,
                   relative_density_range=[min(density),max(density)],final_pressure_Pa=states[-1]['pressure'],
                   minimum_endpoint_storage_m3=min(min(s['geometry']['volume']) for s in states))
    print(json.dumps({k:receipt[k] for k in ['input','independent_audit','maxima','relative_density_range']}),flush=True)
    return receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_lip_profile_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/lip-profile/audit')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/lip_profile_reference.json')
    parser.add_argument('--sanitizer',action='store_true');args=parser.parse_args()
    args.binary=args.binary.resolve();args.output=args.output.resolve();args.report=args.report.resolve();args.output.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((ROOT/'out/lip-profile/nato.csv').open()));geometric=geometry(rows,args.binary,args.output,args.sanitizer)
    trajectories=[];refinements=[]
    for name,start,initial,inlet in [('closing',1284,0.,0.),('opening',1374,700.,700.)]:
        profiles=[vector(rows[i]) for i in range(start,start+3)];previous=None;changes=[]
        for subdivisions in ([16] if args.sanitizer else [16,32,64,128,256,512,1024,2048,4096,8192]):
            r=sequence(name,profiles,[.001,.001],subdivisions,initial,inlet,args.binary,args.output,audit=subdivisions<=64);trajectories.append(r)
            pressure=np.array(r['final_pressure_Pa'])
            if previous is not None:changes.append(dict(subdivisions=subdivisions,maximum_endpoint_difference_Pa=float(max(abs(pressure-previous)))))
            previous=pressure
        refinements.append(dict(name=name,changes=changes))
    sealed=sequence('sealed-storage',[vector(rows[1285])]*2,[.001],64,700.,700.,args.binary,args.output)
    assert max(abs(np.array(sealed['final_pressure_Pa'])-700))<1e-8
    sources=['src/lip_profile.cpp','include/ducklips/lip_profile.hpp','src/articulation.cpp','include/ducklips/articulation.hpp',
             'tools/lip_profile_probe.cpp','tools/lip_profile_reference.py','tools/tongue_reach_audit.cpp','src/moving_graph.cpp','include/ducklips/moving_graph.hpp',
             'tools/planar_inertive_motion_reference.py','bodies/reference.body','samples/nato.dict','out/lip-profile/nato.csv','out/lip-profile/airway.csv','out/lip-profile/geometry.csv']
    report=dict(use='existing continuous axial lip law; volume and inertance reference, no new speech',
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],geometry=geometric,trajectories=trajectories,
                sealed_storage=sealed,refinement=refinements,
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/f for f in sources),args.binary,ROOT/'build/ducklips_tongue_reach_audit']},
                versions=dict(numpy=np.__version__,scipy=scipy.__version__))
    args.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(refinements,indent=2))


if __name__=='__main__':main()
