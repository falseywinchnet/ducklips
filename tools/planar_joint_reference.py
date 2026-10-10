#!/usr/bin/env python3
"""Independent unreduced sparse Radau stages for the actual moving airway."""
import numpy as np
from scipy.sparse import lil_matrix,csc_matrix,diags
from scipy.sparse.linalg import splu
from shapely.geometry import LineString
from planar_source_reference import RHO,BULK,B,pose_at,lip_metrics,check_geometry,exact_regions,check_boundaries
from lip_viscous_reference import lip_segments
from lip_outlet_viscous_reference import outlet_segments
from slit_thermal_reference import spectrum as thermal_spectrum
from slit_viscous_reference import spectrum as viscous_spectrum
from radiation_boundary_reference import inverse as radiation_inverse,COEFF
A=np.array([[5/12,-1/12],[3/4,1/4]]);AI=np.linalg.inv(A);WEIGHT=np.array([.75,.25])

def audit_joint(data,p,q,profiles,epi,geometry_check,radiation=False,thermal_width=0.,viscous_width=0.,lip_viscous_width=0.,lip_viscous_order=8,outlet_viscous=False,stages=2):
    assert stages in [2,3]
    if stages==3:
        r=np.sqrt(6);a=np.array([[(88-7*r)/360,(296-169*r)/1800,(-2+3*r)/225],[(296+169*r)/1800,(88+7*r)/360,(-2-3*r)/225],[(16-r)/36,(16+r)/36,1/9]])
    else:a=A
    ai=np.linalg.inv(a);weight=a[-1];stability=weight[:,None]*a+weight[None,:]*a.T-np.outer(weight,weight)
    mouth=p['boundary_edges']['mouth_end'];entry=p['boundary_edges']['supraglottal_entry'];md=.5*(p['x_edges'][-1]-p['x_edges'][0]);ed=.5*(p['y_edges'][1]-p['y_edges'][0]);cd=.5*(p['y_edges'][-1]-p['y_edges'][0])
    line=LineString([(p['x_edges'][0],.09),(p['x_edges'][-1],.09)]);checks=0
    tw,tp=thermal_spectrum();tw*=.4;vw,vp=viscous_spectrum()
    def geometry(record):
        nonlocal checks
        t=record['fraction'];pose=pose_at(p,q,record);check_geometry(pose,record);regions,cells,error=exact_regions(pose,record)
        boundary,errors=check_boundaries(pose,record,[mouth,entry],regions);lip=lip_metrics((1-t)*profiles[0]+t*profiles[1]);area,length=(1-t)*epi[0]+t*epi[1]
        volume=np.r_[[r.area*p['depth'] for r in regions],lip['volume'],area*length]
        cap=np.array([boundary[(mouth,i)][0] for i in range(2)]);share=cap/sum(cap)
        mouth_g=np.zeros(2) if lip['conductance'][0]==0 else 1/(RHO*md/cap+1/(lip['conductance'][0]*share))
        entry_g=1/(RHO*(ed/boundary[(entry,0)][0]+.5*length[0]/area[0]))
        internal=regions[0].intersection(regions[1]).intersection(line).length*p['depth']
        g=np.r_[internal/(RHO*cd),mouth_g,lip['conductance'][1:],entry_g,1/(RHO*(.5*length[:-1]/area[:-1]+.5*length[1:]/area[1:]))]
        checks+=1;return volume,volume/BULK,g
    def config(t,c,g):
        area,length=(1-t)*epi[0]+t*epi[1];no=[];bo=[];nk=[];nl=[];bk=[];bl=[];resistance=np.zeros(15)
        if thermal_width:
            for i in range(8):no.extend([6+i]*24);nk.extend(c[6+i]*tw);nl.extend(2.3e-5/(area[i]/(2*thermal_width))**2*tp)
        sites=[]
        if viscous_width:
            for i in range(7):
                for j in [i,i+1]:sites.append(dict(branch=8+i,segment_inertance=RHO*.5*length[j]/area[j],half_gap=area[j]/(2*viscous_width)))
        if lip_viscous_width:sites+=lip_segments((1-t)*profiles[0]+t*profiles[1],lip_viscous_order,lip_viscous_width)
        if outlet_viscous:sites+=outlet_segments((1-t)*profiles[0]+t*profiles[1],lip_viscous_order,lip_viscous_width)
        for site in sites:
            l=site['segment_inertance'];gap=site['half_gap'];owner=site['branch'];rate=(1.86e-5/RHO)/gap**2 if l else 0.
            bo.extend([owner]*24);bk.extend(l*vw/vp);bl.extend(rate*vp);resistance[owner]+=3*l*rate
        if radiation:
            profile=(1-t)*profiles[0]+t*profiles[1];a=np.clip(profile[2],0,1);minimum=min(profile[1],profile[6])+(profile[1]-min(profile[1],profile[6]))*(1-a)**2*(1+2*a)
            h=radiation_inverse(minimum);bo.append(6);bk.append(1/h if h else 0.);bl.append(1/(COEFF*h) if h else 0.)
        return np.array(no,dtype=int),np.array(bo,dtype=int),np.array(nk),np.array(nl),np.array(bk),np.array(bl),resistance
    def memories(record):
        nm=np.array([z for site in record.get('thermal_sites',[]) for z in site['charge']],dtype=float)
        bm=np.array([z for site in record.get('viscous_sites',[]) for z in site['momentum']]+([record['radiation_state']['momentum']] if radiation else []),dtype=float)
        return nm,bm
    def metric(c,g,nk,bk):
        return np.r_[1/c,g,np.divide(1,nk,out=np.zeros_like(nk),where=nk>0),np.divide(1,bk,out=np.zeros_like(bk),where=bk>0)]
    def state(record):
        nm,bm=memories(record);return np.r_[record['charge'],record['momentum'],nm,bm]
    def stored(record):
        result=np.array(record['volume'])+record['charge']
        for site in record.get('thermal_sites',[]):result[site['node']]+=sum(site['charge'])
        return result
    def step_reference(x,c,g,configuration,dt,drive):
        no,bo,nk,nl,bk,bl,resistance=configuration;n=14;m=15;count=n+m+len(no)+len(bo)
        matrix=lil_matrix((stages*count,stages*count));right=np.zeros((count,stages));sqc=np.sqrt(c);sqg=np.sqrt(g)
        def add(i,j,block):matrix[stages*i:stages*(i+1),stages*j:stages*(j+1)]+=block
        for i in range(n):add(i,i,ai);right[i]=ai@np.ones(stages)*x[i]/sqc[i]+dt*drive[i]/sqc[i]
        for e in range(m):
            if g[e]==0:add(n+e,n+e,np.eye(stages));continue
            add(n+e,n+e,ai+dt*resistance[e]*g[e]*np.eye(stages));right[n+e]=ai@np.ones(stages)*x[n+e]*sqg[e]
            for i in np.nonzero(B[:,e])[0]:
                block=dt*B[i,e]*sqg[e]/sqc[i]*np.eye(stages);add(i,n+e,block);add(n+e,i,-block)
        for node,owners,storage,rates,offset in [(True,no,nk,nl,n+m),(False,bo,bk,bl,n+m+len(no))]:
            for j,(owner,k,rate) in enumerate(zip(owners,storage,rates)):
                index=offset+j
                if k==0:
                    assert x[index]==0;add(index,index,np.eye(stages));continue
                add(index,index,a+np.eye(stages)/(dt*rate));right[index]=-x[index]/np.sqrt(k)
                if node:
                    scale=np.sqrt(k)/sqc[owner];add(index,owner,-scale*np.eye(stages));add(owner,index,scale*np.eye(stages))
                elif g[owner]>0:
                    scale=np.sqrt(k)*sqg[owner];add(index,n+owner,-scale*np.eye(stages));add(n+owner,index,scale*np.eye(stages))
        matrix=csc_matrix(matrix);scale=1/np.asarray(abs(matrix).max(axis=1).toarray()).ravel()
        solution=splu((diags(scale)@matrix).tocsc()).solve(right.ravel()*scale).reshape(count,stages)
        pressure=solution[:n]/sqc[:,None];flow=solution[n:n+m]*sqg[:,None]
        result=x.copy();result[:n]=c*pressure[:,-1];result[n:n+m]=np.divide(flow[:,-1],g,out=np.zeros(m),where=g>0)
        heat_node=heat_branch=0.;df=np.zeros((count,stages));df[:n]=drive[:,None]-B@flow
        drop=B.T@pressure;df[n:n+m]=drop-resistance[:,None]*flow
        for node,owners,storage,rates,offset in [(True,no,nk,nl,n+m),(False,bo,bk,bl,n+m+len(no))]:
            for j,(owner,k,rate) in enumerate(zip(owners,storage,rates)):
                index=offset+j
                if k==0:result[index]=0;continue
                f=solution[index]*np.sqrt(k)/dt;df[index]=f
                value=pressure[owner] if node else flow[owner]
                result[index]=k*value[-1]-f[-1]/rate if dt*rate>1 else x[index]+dt*(a[-1]@f)
                heat=dt*(weight@(f*f))/k/rate
                if node:df[owner]-=f;heat_node+=heat
                else:df[n+owner]-=f;heat_branch+=heat
        h=metric(c,g,nk,bk);heat_dc=dt*np.sum(resistance*(flow*flow@weight));numeric=.5*dt*dt*np.sum(h*np.einsum("ni,ij,nj->n",df,stability,df))
        work=dt*np.sum(drive*(pressure@weight));return result,pressure@weight,flow@weight,heat_node+heat_branch+heat_dc,numeric,work
    old=data['initial'];reference_state=state(old);initial_mass=sum(stored(old))+sum(old['continuous_source']['mass']);transfer=0.;max_mass=max_local=max_pressure=max_memory=max_heat=max_numeric=max_work=0.;energy=physical=numerical=geometry_work=0.
    if geometry_check:geometry(old['geometry'])
    for s in data['steps']:
        dt=s['dt'];flow=np.array(s['average_flow']);before=stored(old);after=stored(s);source_drive=np.zeros(14);source_drive[-1]=s['source_flow']
        max_local=max(max_local,float(max(abs(after-before-dt*(-B@flow+source_drive)))))
        transfer+=dt*(s['source_reservoir_flow']-flow[6]);mass=sum(after)+sum(s['continuous_source']['mass']);max_mass=max(max_mass,abs(mass-initial_mass-transfer))
        receipt=s['joint_passive'];assert abs(receipt['energy_residual'])<1e-15 and abs(receipt['air_residual'])<1e-18
        assert receipt['numerical_dissipation']>=0 and min(receipt[k] for k in ['viscous_dissipation','thermal_dissipation','radiated_energy'])>=0
        assert abs(s['source_work_residual'])<1e-12 and abs(s['source_port_residual'])<1e-17
        physical+=sum(receipt[k] for k in ['viscous_dissipation','thermal_dissipation','radiated_energy']);numerical+=receipt['numerical_dissipation'];geometry_work+=sum(receipt[k] for k in ['graph_geometry_work','viscous_geometry_work','thermal_geometry_work','radiation_geometry_work'])
        if geometry_check:
            volume,c,g=geometry(s['stage1']);endv,endc,endg=geometry(s['geometry']);assert s['stage1']['fraction']==s['stage2']['fraction']
            assert max(abs(endv-s['volume']))<1e-17 and max(abs(endg-s['conductance']))<1e-10
            middle=config(s['stage1']['fraction'],c,g);old_conf=config(old['geometry']['fraction'],np.array(old['compliance']),np.array(old['conductance']));end_conf=config(s['geometry']['fraction'],endc,endg)
            x=reference_state.copy();midmetric=metric(c,g,middle[2],middle[4]);oldmetric=metric(np.array(old['compliance']),np.array(old['conductance']),old_conf[2],old_conf[4]);endmetric=metric(endc,endg,end_conf[2],end_conf[4])
            first_work=.5*np.sum(x*x*(midmetric-oldmetric));x[midmetric==0]=0.
            drive=-(endv-np.array(old['volume']))/dt+source_drive
            y,pavg,qavg,heat,num,work=step_reference(x,c,g,middle,dt,drive)
            second_work=.5*np.sum(y*y*(endmetric-midmetric));y[endmetric==0]=0.
            reference_state=y.copy();actual=state(s);max_pressure=max(max_pressure,float(max(abs(y[:14]/endc-s['pressure']))),float(max(abs(pavg-s['average_pressure']))))
            max_memory=max(max_memory,float(max(abs((y-actual)*np.sqrt(np.maximum(endmetric,0))))))
            max_heat=max(max_heat,abs(heat-sum(receipt[k] for k in ['viscous_dissipation','thermal_dissipation','radiated_energy'])))
            max_numeric=max(max_numeric,abs(num-receipt['numerical_dissipation']))
            max_work=max(max_work,abs(first_work+second_work-sum(receipt[k] for k in ['graph_geometry_work','viscous_geometry_work','thermal_geometry_work','radiation_geometry_work'])))
            expected_change=.5*np.sum(actual*actual*endmetric-state(old)**2*oldmetric)
            assert abs(expected_change-receipt['energy_change'])<1e-15
        old=s
    assert max_mass<1e-17 and max_local<1e-17,(max_mass,max_local)
    assert max_pressure<1e-5 and max_memory<1e-10 and max_heat<1e-15 and max_numeric<1e-15 and max_work<1e-15,(max_pressure,max_memory,max_heat,max_numeric,max_work)
    return dict(steps=len(data['steps']),geometry_snapshots_checked=checks,maximum_local_mass_residual_m3=max_local,maximum_combined_mass_residual_m3=max_mass,
        maximum_independent_pressure_error_Pa=max_pressure if geometry_check else None,maximum_independent_energy_coordinate_error=max_memory if geometry_check else None,
        maximum_independent_physical_dissipation_error_J=max_heat if geometry_check else None,maximum_independent_numerical_dissipation_error_J=max_numeric if geometry_check else None,maximum_independent_geometry_work_error_J=max_work if geometry_check else None,
        physical_dissipation_J=physical,numerical_dissipation_J=numerical,geometry_work_J=geometry_work,net_external_transfer_m3=transfer,initial_total_mass_m3=initial_mass,final_total_mass_m3=sum(stored(old))+sum(old['continuous_source']['mass']))
