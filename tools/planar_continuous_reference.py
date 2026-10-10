#!/usr/bin/env python3
"""Connect the shared continuous source to actual moving Bravo geometry."""
import argparse
import csv
import json
import subprocess
from pathlib import Path
import numpy as np
from shapely.geometry import LineString
from scipy.linalg import expm
from slit_thermal_reference import spectrum
from slit_viscous_reference import spectrum as viscous_spectrum
from lip_viscous_reference import lip_segments
from radiation_boundary_reference import inverse as radiation_inverse, COEFF as radiation_coefficient
from planar_source_reference import (ROOT,RHO,BULK,B,EDGES,setup,pose_at,PORTS,
    lip_metrics,vector,write_input,dense_step,sha,check_geometry,exact_regions,check_boundaries)


def source_vector(s):
    return np.array(sum([s[k] for k in ['effort','flow','lower_position','lower_velocity','upper_position','upper_velocity']],[]))

def decay_exponential(matrix):
    """Full-matrix reference: scale a Taylor expm1, then double its increment.

    Repeated squaring of I+D loses the small slow-mode increment when the
    fastest viscous modes demand many scalings. D <- 2D+D@D retains it.
    With norm <= 1/8 the omitted Taylor tail is below 1e-30 before doubling.
    This is independent of the native symmetric eigensystem propagator.
    """
    norm=np.linalg.norm(matrix,1)
    if norm==0:return np.eye(len(matrix))
    scalings=max(0,int(np.ceil(np.log2(norm/.125))))
    scaled=matrix/(2.**scalings);increment=scaled.copy();term=scaled.copy()
    for k in range(2,18):
        term=term@scaled/k;increment+=term
    for _ in range(scalings):increment=2*increment+increment@increment
    return np.eye(len(matrix))+increment

def run(name,p,q,profiles,epi,n,pressure,previous,binary,out,radiation=False,thermal_width=0.,viscous_width=0.,lip_viscous_width=0.,lip_viscous_order=8):
    path=out/f'{name}-{n}.txt'
    write_input(p,q,n,[pressure]*2,PORTS,path)
    mouth=p['boundary_edges']['mouth_end'];entry=p['boundary_edges']['supraglottal_entry']
    md=.5*(p['x_edges'][-1]-p['x_edges'][0]);ed=.5*(p['y_edges'][1]-p['y_edges'][0])
    with path.open('a') as f:
        f.write('4 2\n');f.writelines(' '.join(map(str,v))+'\n' for v in profiles)
        f.write(' '.join(map(str,[pressure]*4))+'\n')
        f.writelines(f'{mouth} {i} {md}\n' for i in range(2))
        f.write('8 1\n');f.writelines(' '.join(map(str,v))+'\n' for v in np.column_stack([epi[0,0],epi[0,1],epi[1,0],epi[1,1]]))
        f.write(' '.join(map(str,[pressure]*8))+'\n')
        f.write(f'{entry} 0 {ed}\n3e-6 700\n1\n')
        if previous is not None:
            for k in ['charge','momentum']:f.write(' '.join(map(str,previous[k]))+'\n')
            f.write(' '.join(map(str,source_vector(previous['continuous_source'])))+'\n')
            if radiation:f.write(str(previous['radiation_state']['momentum'])+'\n')
            if thermal_width:
                for site in previous['thermal_sites']:f.write(' '.join(map(str,site['charge']))+'\n')
            if viscous_width or lip_viscous_width:
                for site in previous['viscous_sites']:f.write(' '.join(map(str,site['momentum']))+'\n')
    command=[str(binary),'--lip-profile','--epilarynx','--continuous-folds',str(ROOT/'bodies/reference.body'),'--compact-geometry']
    if radiation:command+=['--radiation']
    if thermal_width:command+=['--thermal-slit-width',str(thermal_width)]
    if viscous_width:command+=['--viscous-slit-width',str(viscous_width)]
    if lip_viscous_width:command+=['--lip-viscous-width',str(lip_viscous_width),'--lip-viscous-order',str(lip_viscous_order)]
    if previous is not None:command+=['--initial-state']
    with path.open() as f,path.with_suffix('.json').open('w') as g:
        result=subprocess.run(command,stdin=f,stdout=g,stderr=subprocess.PIPE,text=True)
    path.with_suffix('.log').write_text(result.stderr)
    assert result.returncode==0,(name,n,result.stderr)
    data=json.loads(path.with_suffix('.json').read_text())
    continuation=None
    if previous is not None:
        for k in ['volume','compliance','conductance','charge','momentum']:
            assert data['initial'][k]==previous[k],k
        if radiation:assert data['initial']['radiation_state']==previous['radiation_state']
        if thermal_width:assert data['initial']['thermal_sites']==previous['thermal_sites']
        if viscous_width or lip_viscous_width:assert data['initial']['viscous_sites']==previous['viscous_sites']
        a=data['initial']['continuous_source'];b=previous['continuous_source']
        mass_error=float(np.max(abs(np.array(a['mass'])-b['mass'])))
        momentum_error=float(np.max(abs(np.array(a['canonical_momentum'])-b['canonical_momentum'])))
        assert mass_error<1e-20 and momentum_error<1e-14,(mass_error,momentum_error)
        continuation=dict(source_mass_roundtrip_m3=mass_error,source_momentum_roundtrip=momentum_error,
            source_energy_roundtrip_J=abs(a['energy']-b['energy']))
    return data,dict(name=name,subdivisions=n,command=command,continuation=continuation,
        files={str(path.with_suffix(ext).relative_to(ROOT)):sha(path.with_suffix(ext)) for ext in ['.txt','.json','.log']})

def audit(data,p,q,profiles,epi,geometry_check,radiation=False,thermal_width=0.,viscous_width=0.,lip_viscous_width=0.,lip_viscous_order=8):
    mouth=p['boundary_edges']['mouth_end'];entry=p['boundary_edges']['supraglottal_entry']
    md=.5*(p['x_edges'][-1]-p['x_edges'][0]);ed=.5*(p['y_edges'][1]-p['y_edges'][0]);cd=.5*(p['y_edges'][-1]-p['y_edges'][0])
    line=LineString([(p['x_edges'][0],.09),(p['x_edges'][-1],.09)])
    checks=0
    def metrics(record):
        nonlocal checks
        fraction=record['fraction'];pose=pose_at(p,q,record)
        check_geometry(pose,record)
        regions,cells,error=exact_regions(pose,record)
        boundary,errors=check_boundaries(pose,record,[mouth,entry],regions)
        lip=lip_metrics((1-fraction)*profiles[0]+fraction*profiles[1]);ea,el=(1-fraction)*epi[0]+fraction*epi[1]
        v=np.r_[[r.area*p['depth'] for r in regions],lip['volume'],ea*el]
        cap=np.array([boundary[(mouth,i)][0] for i in range(2)]);share=cap/sum(cap)
        mouth_g=np.zeros(2) if lip['conductance'][0]==0 else 1/(RHO*md/cap+1/(lip['conductance'][0]*share))
        entry_g=1/(RHO*(ed/boundary[(entry,0)][0]+.5*el[0]/ea[0]))
        internal=regions[0].intersection(regions[1]).intersection(line).length*p['depth']
        g=np.r_[internal/(RHO*cd),mouth_g,lip['conductance'][1:],entry_g,1/(RHO*(.5*el[:-1]/ea[:-1]+.5*el[1:]/ea[1:]))]
        checks+=1
        return v,v/BULK,g
    def exit_area(fraction):
        length,base,activation,place,core,width,target=(1-fraction)*profiles[0]+fraction*profiles[1]
        activation=np.clip(activation,0,1);target=min(base,target)
        minimum=target+(base-target)*(1-activation)**2*(1+2*activation)
        u=max(abs(place)-core,0)/width
        return base if u>np.sqrt(-np.log(1e-4)) else minimum+(base-minimum)*(-np.expm1(-u*u))
    def exchange(duct,memory,g,area,dt):
        if area==0:
            assert g==0 and memory==0
            return duct,memory,0.
        h=radiation_inverse(area);rate=(g+h)/(radiation_coefficient*h*h)
        # Independent equilibrium/deviation propagation, rather than the
        # native transferred-momentum update. Very stiff modes decay to zero.
        common=(duct+memory)/(g+h);decay=np.exp(-dt*rate)
        a=h*common+(duct-h*common)*decay;b=g*common+(memory-g*common)*decay
        heat=.5*g*(duct*duct-a*a)+.5*h*(memory*memory-b*b)
        return a,b,heat
    def stored(record):
        result=np.array(record['volume'])+record['charge']
        if thermal_width:
            for site in record['thermal_sites']:
                assert abs(sum(site['charge'])-site['stored_air'])<1e-20
                result[site['node']]+=site['stored_air']
        return result
    thermal_k,thermal_p=spectrum();thermal_k*=.4
    thermal_metric=np.r_[1.,np.sqrt(thermal_k)]
    thermal_matrix=np.zeros((25,25));thermal_matrix[0,0]=-np.sum(thermal_k*thermal_p)
    thermal_matrix[0,1:]=thermal_matrix[1:,0]=thermal_p*np.sqrt(thermal_k);thermal_matrix[1:,1:]=-np.diag(thermal_p)
    exponentials={}
    def thermal_exchange(primary,memory,compliance,gap,dt):
        time=dt*2.3e-5/gap**2
        if time not in exponentials:exponentials[time]=expm(thermal_matrix*time)
        energy_state=np.r_[primary,memory]/thermal_metric
        next_state=exponentials[time]@energy_state
        heat=.5*(energy_state@energy_state-next_state@next_state)/compliance
        physical=next_state*thermal_metric
        return physical[0],physical[1:],heat
    visc_w,visc_p=viscous_spectrum();visc_exponentials={}
    def viscous_exchange(primary,memory,site,dt):
        l=site['segment_inertance'];g=site['branch_conductance'];a=site['half_gap'];time=dt*(1.86e-5/RHO)/a**2;ratio=l*g
        key=(ratio,time)
        if key not in visc_exponentials:
            matrix=np.zeros((25,25));matrix[0,0]=-ratio*(3+sum(visc_w));matrix[0,1:]=matrix[1:,0]=np.sqrt(ratio*visc_w*visc_p);matrix[1:,1:]=-np.diag(visc_p)
            visc_exponentials[key]=decay_exponential(matrix*time)
        metric=np.sqrt(np.r_[g,visc_p/(l*visc_w)]);before=np.r_[primary,memory]*metric;after=visc_exponentials[key]@before
        heat=.5*(before@before-after@after);physical=after/metric
        return physical[0],physical[1:],heat
    old=data['initial'];initial_mass=sum(stored(old))+sum(old['continuous_source']['mass'])
    thermal_memory=np.array([site['charge'] for site in old.get('thermal_sites',[])],dtype=float)
    max_thermal=max_thermal_heat=max_thermal_work=thermal_heat=thermal_work=0.
    viscous_memory=np.array([site['momentum'] for site in old.get('viscous_sites',[])],dtype=float)
    max_viscous=max_viscous_heat=max_viscous_work=viscous_heat=viscous_work=0.
    transfer=0.;max_mass=0.;max_local=0.;max_graph=0.;max_work=0.;max_port=0.;max_column=0.
    max_radiation=0.;max_radiation_heat=0.;max_radiation_work=0.;radiated=0.;radiation_work=0.
    rad_momentum=old['radiation_state']['momentum'] if radiation else 0.
    z=np.r_[old['charge'],old['momentum']]
    if geometry_check:
        ov,c,g=metrics(old['geometry'])
        assert np.max(abs(ov-old['volume']))<1e-17
    for s in data['steps']:
        dt=s['dt'];flow=np.array(s['average_flow'])
        if radiation:
            assert s['radiated_energy']>=0
            radiated+=s['radiated_energy'];radiation_work+=s['radiation_geometry_work']
            assert abs(s['radiation_work_residual'])<1e-17
            assert abs(s['radiation_state']['area']-exit_area(s['geometry']['fraction']))<1e-18
        network=stored(s);prior=stored(old)
        if viscous_width or lip_viscous_width:
            assert s['viscous_dissipation']>=0 and abs(s['viscous_work_residual'])<1e-16
            viscous_heat+=s['viscous_dissipation'];viscous_work+=s['viscous_geometry_work']
            fraction=s['geometry']['fraction'];areas,lengths=(1-fraction)*epi[0]+fraction*epi[1]
            throat_count=14 if viscous_width else 0
            lip_reference=lip_segments((1-fraction)*profiles[0]+fraction*profiles[1],lip_viscous_order,lip_viscous_width) if lip_viscous_width else []
            for i,site in enumerate(s['viscous_sites']):
                if i<throat_count:
                    j=i//2+i%2;branch=8+i//2
                    assert site['branch']==branch and abs(site['half_gap']-areas[j]/(2*viscous_width))<1e-15
                    assert abs(site['segment_inertance']-RHO*.5*lengths[j]/areas[j])<1e-10
                else:
                    reference=lip_reference[i-throat_count];branch=reference['branch'];assert site['branch']==branch
                    for key in ['position','half_gap','segment_inertance','branch_conductance']:
                        assert abs(site[key]/reference[key]-1)<1e-10,(key,site[key],reference[key])
                assert site['branch_conductance']==s['conductance'][branch]
        if thermal_width:
            assert s['thermal_dissipation']>=0 and abs(s['thermal_work_residual'])<1e-16 and abs(s['thermal_air_residual'])<1e-20
            thermal_heat+=s['thermal_dissipation'];thermal_work+=s['thermal_geometry_work']
            for i,site in enumerate(s['thermal_sites']):
                fraction=s['geometry']['fraction'];a,l=(1-fraction)*epi[0]+fraction*epi[1]
                assert site['node']==6+i and abs(site['half_gap']-a[i]/(2*thermal_width))<1e-15
                assert abs(site['compliance']-a[i]*l[i]/BULK)<1e-22
        local=network-prior+dt*B@flow;local[-1]-=dt*s['source_flow']
        max_local=max(max_local,float(np.max(abs(local))))
        source_mass=sum(s['continuous_source']['mass']);old_source_mass=sum(old['continuous_source']['mass'])
        max_local=max(max_local,abs(source_mass-old_source_mass-dt*(s['source_reservoir_flow']-s['source_flow'])))
        transfer+=dt*(s['source_reservoir_flow']-flow[6])
        max_mass=max(max_mass,abs(sum(network)+source_mass-initial_mass-transfer))
        max_work=max(max_work,abs(s['source_work_residual']));max_port=max(max_port,abs(s['source_port_residual']))
        t=s['geometry']['fraction'];ea,el=(1-t)*epi[0]+t*epi[1]
        max_column=max(max_column,abs(s['continuous_source']['outlet_inertance']-RHO*.5*el[-1]/ea[-1]))
        if geometry_check:
            first_g=g[6];old_c=c.copy()
            _,c1,g1=metrics(s['stage1']);_,c2,g2=metrics(s['stage2']);v,c,g=metrics(s['geometry'])
            assert np.max(abs(v-s['volume']))<1e-17
            assert np.max(abs(c-s['compliance']))<1e-20
            assert np.max(abs(g-s['conductance']))<1e-10
            z[14:][(g1==0)&(g2==0)]=0
            step_thermal_heat=step_thermal_work=0.;step_viscous_heat=step_viscous_work=0.
            if thermal_width:
                for i,site in enumerate(old['thermal_sites']):
                    z[6+i],thermal_memory[i],h=thermal_exchange(z[6+i],thermal_memory[i],old_c[6+i],site['half_gap'],dt/2)
                    step_thermal_heat+=h
            if viscous_width or lip_viscous_width:
                for i,site in enumerate(old['viscous_sites']):
                    j=14+site['branch'];z[j],viscous_memory[i],h=viscous_exchange(z[j],viscous_memory[i],site,dt/2);step_viscous_heat+=h
            if radiation:
                a0=exit_area(old['geometry']['fraction']);am=exit_area(.5*(old['geometry']['fraction']+s['geometry']['fraction']));a1=exit_area(s['geometry']['fraction'])
                work=.5*(radiation_inverse(am)-radiation_inverse(a0))*rad_momentum**2
                if am==0:rad_momentum=0.
                z[20],rad_momentum,heat=exchange(z[20],rad_momentum,first_g,am,dt/2)
            drive=-(v-ov)/dt;drive[-1]+=s['source_flow']
            z,af,ap=dense_step(z,c1,c2,g1,g2,drive,dt,B)
            if radiation:
                z[20],rad_momentum,second_heat=exchange(z[20],rad_momentum,g[6],am,dt/2)
                work+=.5*(radiation_inverse(a1)-radiation_inverse(am))*rad_momentum**2
                if a1==0:rad_momentum=0.
                max_radiation=max(max_radiation,abs(rad_momentum-s['radiation_state']['momentum']))
                max_radiation_heat=max(max_radiation_heat,abs(heat+second_heat-s['radiated_energy']))
                max_radiation_work=max(max_radiation_work,abs(work-s['radiation_geometry_work']))
            if viscous_width or lip_viscous_width:
                for i in reversed(range(len(s['viscous_sites']))):
                    site=s['viscous_sites'][i];previous=old['viscous_sites'][i]
                    step_viscous_work+=.5*np.sum(visc_p/visc_w*viscous_memory[i]**2)*(1/site['segment_inertance']-1/previous['segment_inertance'])
                    j=14+site['branch'];z[j],viscous_memory[i],h=viscous_exchange(z[j],viscous_memory[i],site,dt/2);step_viscous_heat+=h
                    max_viscous=max(max_viscous,float(np.max(abs(viscous_memory[i]-site['momentum']))))
                max_viscous_heat=max(max_viscous_heat,abs(step_viscous_heat-s['viscous_dissipation']));max_viscous_work=max(max_viscous_work,abs(step_viscous_work-s['viscous_geometry_work']))
            if thermal_width:
                for i,site in enumerate(s['thermal_sites']):
                    step_thermal_work+=.5*np.sum(thermal_memory[i]**2/thermal_k)*(1/c[6+i]-1/old_c[6+i])
                    z[6+i],thermal_memory[i],h=thermal_exchange(z[6+i],thermal_memory[i],c[6+i],site['half_gap'],dt/2)
                    step_thermal_heat+=h
                    max_thermal=max(max_thermal,float(np.max(abs(thermal_memory[i]-site['charge']))))
                max_thermal_heat=max(max_thermal_heat,abs(step_thermal_heat-s['thermal_dissipation']))
                max_thermal_work=max(max_thermal_work,abs(step_thermal_work-s['thermal_geometry_work']))
            z[14:][g==0]=0
            max_graph=max(max_graph,float(np.max(abs(z[:14]/c-s['pressure']))),float(np.max(abs(ap-s['average_pressure']))))
            assert np.max(abs(z[14:]-s['momentum']))<1e-9,(s['geometry']['fraction'],float(np.max(abs(z[14:]-s['momentum']))),max_graph,max_thermal)
            ov=v
        # A step ending at closure can transport air before its endpoint.
        # Exact zero average flow applies to an entirely sealed interval.
        if old['conductance'][6]==0 and s['conductance'][6]==0:assert flow[6]==0
        old=s
    assert max_radiation<1e-12 and max_radiation_heat<1e-17 and max_radiation_work<1e-17,(max_radiation,max_radiation_heat,max_radiation_work)
    assert max_viscous<1e-10 and max_viscous_heat<1e-15 and max_viscous_work<1e-15,(max_viscous,max_viscous_heat,max_viscous_work)
    assert max_thermal<1e-18 and max_thermal_heat<1e-16 and max_thermal_work<1e-16,(max_thermal,max_thermal_heat,max_thermal_work)
    assert max_local<1e-17 and max_mass<1e-17,(max_local,max_mass)
    assert max_work<1e-12 and max_port<1e-18,(max_work,max_port)
    assert max_column<1e-10 and max_graph<1e-5,(max_column,max_graph)
    return dict(lip_viscous_width=lip_viscous_width,lip_viscous_order=lip_viscous_order,viscous_slit_width=viscous_width,viscous_dissipation_J=viscous_heat,viscous_geometry_work_J=viscous_work,maximum_independent_viscous_momentum_error_Pa_s=max_viscous if geometry_check else None,maximum_independent_viscous_energy_error_J=max_viscous_heat if geometry_check else None,maximum_independent_viscous_geometry_work_error_J=max_viscous_work if geometry_check else None,
        thermal_slit_width=thermal_width,thermal_dissipation_J=thermal_heat,thermal_geometry_work_J=thermal_work,
        maximum_independent_thermal_charge_error_m3=max_thermal if geometry_check else None,maximum_independent_thermal_energy_error_J=max_thermal_heat if geometry_check else None,maximum_independent_thermal_geometry_work_error_J=max_thermal_work if geometry_check else None,
        steps=len(data['steps']),geometry_snapshots_checked=checks,maximum_local_mass_residual_m3=max_local,
        maximum_combined_mass_residual_m3=max_mass,maximum_source_work_residual_J=max_work,
        maximum_port_work_residual_J=max_port,maximum_column_metric_error=max_column,
        maximum_independent_graph_pressure_error_Pa=max_graph if geometry_check else None,
        radiation_enabled=radiation,radiated_energy_J=radiated,radiation_geometry_work_J=radiation_work,
        maximum_radiation_momentum_error_Pa_s=max_radiation if geometry_check else None,maximum_radiated_energy_error_J=max_radiation_heat if geometry_check else None,maximum_radiation_geometry_work_error_J=max_radiation_work if geometry_check else None,
        net_external_transfer_m3=transfer,initial_total_mass_m3=initial_mass,
        final_total_mass_m3=sum(stored(old))+sum(old['continuous_source']['mass']))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--subdivisions',nargs='+',type=int,default=[256,512,1024])
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_planar_inertive_motion_probe')
    parser.add_argument('--radiation',action='store_true')
    parser.add_argument('--thermal-slit-width',type=float,default=0.)
    parser.add_argument('--viscous-slit-width',type=float,default=0.)
    parser.add_argument('--lip-viscous-width',type=float,default=0.)
    parser.add_argument('--lip-viscous-order',type=int,choices=[4,8,16,32],default=8)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    assert np.isfinite(args.lip_viscous_width) and args.lip_viscous_width>=0
    assert np.isfinite(args.viscous_slit_width) and args.viscous_slit_width>=0
    assert np.isfinite(args.thermal_slit_width) and args.thermal_slit_width>=0
    if args.output is None:args.output=ROOT/('out/planar-lip-viscous' if args.lip_viscous_width else 'out/planar-viscous' if args.viscous_slit_width else 'out/planar-thermal' if args.thermal_slit_width else 'out/planar-radiation' if args.radiation else 'out/planar-continuous')
    if args.report is None:args.report=ROOT/('reports/planar_lip_viscous_reference.json' if args.lip_viscous_width else 'reports/planar_viscous_reference.json' if args.viscous_slit_width else 'reports/planar_thermal_reference.json' if args.thermal_slit_width else 'reports/planar_radiation_reference.json' if args.radiation else 'reports/planar_continuous_reference.json')
    args.binary=args.binary.resolve();args.output=args.output.resolve();args.report=args.report.resolve();args.output.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((ROOT/'out/lip-profile/nato.csv').open()))
    airway=list(csv.DictReader((ROOT/'out/nato-spatial-airway.csv').open()))
    runs=[];refinements=[];inputs={}
    for name,start,pressure in [('closing',1284,0.),('opening',1374,700.)]:
        last=None
        for n in args.subdivisions:
            previous=None;first_mass=None;net=0.
            for interval in range(2):
                index=start+interval;p,q,paths=setup(index/1000)
                profiles=np.array([vector(rows[i]) for i in [index,index+1]])
                epi=np.array([[[float(airway[t][f'{field}_{i}']) for i in reversed(range(8))] for field in ['area','length']] for t in [index,index+1]])
                data,receipt=run(f'{name}-{interval}',p,q,profiles,epi,n,pressure,previous,args.binary,args.output,args.radiation,args.thermal_slit_width,args.viscous_slit_width,args.lip_viscous_width,args.lip_viscous_order)
                receipt.update(audit(data,p,q,profiles,epi,n==args.subdivisions[0],args.radiation,args.thermal_slit_width,args.viscous_slit_width,args.lip_viscous_width,args.lip_viscous_order));runs.append(receipt)
                if first_mass is None:first_mass=receipt['initial_total_mass_m3']
                net+=receipt['net_external_transfer_m3'];previous=data['steps'][-1]
                inputs.update({str(path.relative_to(ROOT)):sha(path) for path in paths})
            balance=abs(receipt['final_total_mass_m3']-first_mass-net);assert balance<1e-17,balance
            final=np.array(previous['pressure'])
            refinement=dict(case=name,subdivisions=n,endpoint_pressure_difference_Pa=None if last is None else float(np.max(abs(final-last))),whole_control_mass_residual_m3=balance,
                final_source=previous['continuous_source'])
            refinements.append(refinement);last=final
            print(json.dumps(dict(case=name,subdivisions=n,endpoint_pressure_difference_Pa=refinement['endpoint_pressure_difference_Pa'],mass_residual_m3=balance)),flush=True)
    paths=[Path(__file__),args.binary,ROOT/'tools/planar_inertive_motion_probe.cpp',ROOT/'src/continuous_glottis.cpp',ROOT/'src/graph_continuous_glottis.cpp',
        ROOT/'bodies/reference.body',ROOT/'out/lip-profile/nato.csv',ROOT/'out/nato-spatial-airway.csv',ROOT/'tools/planar_source_reference.py',ROOT/'tools/planar_lip_profile_reference.py',ROOT/'tools/lip_profile_reference.py',ROOT/'tools/planar_inertive_motion_reference.py']
    paths += [ROOT/'src/radiation_boundary.cpp',ROOT/'include/ducklips/radiation_boundary.hpp',ROOT/'tools/radiation_boundary_reference.py']
    paths += [ROOT/'src/slit_thermal.cpp',ROOT/'include/ducklips/slit_thermal.hpp',ROOT/'include/ducklips/thermal_spectrum.hpp',ROOT/'tools/slit_thermal_reference.py']
    paths += [ROOT/'src/slit_viscous.cpp',ROOT/'include/ducklips/slit_viscous.hpp',ROOT/'include/ducklips/viscous_spectrum.hpp',ROOT/'tools/slit_viscous_reference.py']
    paths += [ROOT/'src/lip_viscous.cpp',ROOT/'include/ducklips/lip_viscous.hpp',ROOT/'tools/lip_viscous_reference.py']
    report=dict(lip_viscous_width=args.lip_viscous_width,lip_viscous_order=args.lip_viscous_order,viscous_slit_width=args.viscous_slit_width,thermal_slit_width=args.thermal_slit_width,radiation_enabled=args.radiation,scope='Continuous shared-inertia gas/tissue source connected to moving Bravo mouth, lip and epilarynx controls. Independent graph solve uses the recorded coupled source flux; it is not an independent whole-source motion oracle. Yielding-wall losses are absent; open internal lip-path viscous quadrature is enabled by lip_viscous_width; explicit throat-internal viscous segments are enabled by viscous_slit_width; explicit throat slit thermal storage is enabled by thermal_slit_width, and aperture radiation by radiation_enabled. No new speech.',runs=runs,refinements=refinements,input_poses=inputs,
        hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})
    args.report.write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
