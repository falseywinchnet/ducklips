#!/usr/bin/env python3
"""Check flat-core/shoulder outlet quadrature through exact closure."""
import csv,hashlib,json,subprocess
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad,quad_vec
from lip_viscous_reference import area_function,RHO,NU,ROOT
from lip_profile_reference import vector

def outlet_geometry(profile):
    length,base,activation,place,core,width,target=profile
    assert place==0 and core>0 and base>0
    activation=np.clip(activation,0,1);minimum=min(base,target)+(base-min(base,target))*(1-activation)**2*(1+2*activation)
    area,_=area_function(profile)
    end=length/4
    volume=quad(lambda x:float(area(x)),0,end,points=[core],epsabs=0,epsrel=2e-12)[0]
    centroid=quad(lambda x:x*float(area(x)),0,end,points=[core],epsabs=0,epsrel=2e-12)[0]/volume
    assert centroid>core
    # Resolve the shrinking shoulder scale without sampling over it. This
    # coordinate is used only by the integration oracle, not state identities.
    scale=width*np.sqrt(minimum/base) if minimum>0 else 0.
    angle=np.arctan((centroid-core)/scale) if minimum>0 else 0.
    def transform(theta):
        u=(scale/width)*np.tan(theta)
        value=minimum+(base-minimum)*(-np.expm1(-u*u))
        return value,scale/np.cos(theta)**2
    if minimum>0:
        shoulder=quad(lambda t:transform(t)[1]/transform(t)[0],0,angle,epsabs=0,epsrel=2e-11)[0]
        conductance=1/(RHO*(core/minimum+shoulder))
    else:conductance=0.
    return area,minimum,centroid,conductance,transform,angle

def outlet_segments(profile,order,width):
    area,minimum,last,g,_,_=outlet_geometry(profile);core=profile[4]
    x,w=leggauss(order);positive=x[order//2:];weights=w[order//2:]
    x=np.column_stack([-positive,positive]).ravel();w=np.repeat(weights,2)
    positions=np.r_[core/2,(core+last)/2+(last-core)/2*x]
    areas=np.r_[minimum,area(positions[1:])]
    inertias=np.r_[RHO*core/minimum if minimum else 0.,RHO*(last-core)/2*w/areas[1:]]
    if g>0:inertias/=sum(inertias)*g
    return [dict(branch=6,position=float(x),segment_inertance=float(l),branch_conductance=float(g),half_gap=float(a/(2*width))) for x,l,a in zip(positions,inertias,areas)]

def correction(time):
    # Stable exact slit correction after removing the plug term. Its Taylor
    # remainder below |time|=1e-4 is below binary64 precision relative to 3.
    time=np.asarray(time,dtype=complex);small=abs(time)<1e-4;u=np.sqrt(np.where(small,1.,time));t=np.tanh(u)
    return np.where(small,3+time/5-time*time/175,time*t/(u-t))

def main():
    out=ROOT/'out/lip-outlet-viscous-reference';out.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((ROOT/'out/lip-profile/nato.csv').open()));profiles=[]
    for start in [1284,1374]:
        for k in range(9):
            t=start+k/4;i=int(t);f=t-i;profiles.append((str(t),(1-f)*vector(rows[i])+f*vector(rows[i+1])))
    for k in range(1,8):
        profile=vector(rows[1284]);profile[2]=1-10.**(-k);profiles.append((f'near-seal-{k}',profile))
    lines=[];cases=[];s=2j*np.pi*np.geomspace(20,200000,129);width=.03
    for label,profile in profiles:
        area,minimum,last,g,transform,angle=outlet_geometry(profile)
        if minimum>0:
            def local(value):
                rate=NU/(value/(2*width))**2
                return RHO/value*rate*correction(s/rate)
            exact=profile[4]*local(minimum)+quad_vec(lambda t:transform(t)[1]*local(transform(t)[0]),0,angle,epsabs=0,epsrel=2e-10)[0]
        else:exact=None
        for order in [4,8,16,32]:
            lines.append(str(order+1)+' '+' '.join(map(str,np.r_[profile,0,last,width])))
            cases.append((label,profile,order,g,exact))
    inp=out/'geometry.input';inp.write_text('\n'.join(lines)+'\n')
    result=subprocess.run([str(ROOT/'build/ducklips_lip_viscous_probe'),'--outlet-fixture'],input=inp.read_text(),text=True,capture_output=True)
    assert result.returncode==0,result.stderr;dest=out/'geometry.txt';dest.write_text(result.stdout)
    errors={n:dict(geometry=0.,correction=0.,inertia_normalization=0.) for n in [4,8,16,32]}
    for line,(label,profile,n,g,exact) in zip(result.stdout.splitlines(),cases,strict=True):
        v=np.fromstring(line,sep=' ');actual=v[3:].reshape(n+1,3);expected=outlet_segments(profile,n,width)
        target=np.array([[r[k] for k in ['position','segment_inertance','half_gap']] for r in expected]);positive=target!=0
        assert np.all(actual[~positive]==0)
        error=float(max(abs(actual[positive]/target[positive]-1)));e=errors[n];e['geometry']=max(e['geometry'],error)
        assert error<2e-9,(label,n,error)
        assert (g==0 and v[0]==0) or abs(v[0]/g-1)<2e-9
        if g>0:
            assert abs(sum(actual[:,1])*v[0]-1)<1e-13
            rate=NU/actual[:,2]**2
            approximate=np.sum(actual[:,1]*rate*correction(s[:,None]/rate),axis=1)
            e['correction']=max(e['correction'],float(max(abs(approximate/exact-1))))
            e['inertia_normalization']=max(e['inertia_normalization'],abs(v[2]))
    assert errors[8]['correction']<.002 and errors[16]['correction']<.0004 and errors[32]['correction']<.0002,errors
    paths=[Path(__file__),ROOT/'tools/lip_viscous_reference.py',ROOT/'tools/lip_profile_reference.py',ROOT/'src/lip_viscous.cpp',ROOT/'include/ducklips/lip_viscous.hpp',ROOT/'src/lip_profile.cpp',ROOT/'tools/lip_viscous_probe.cpp',ROOT/'build/ducklips_lip_viscous_probe',ROOT/'out/lip-profile/nato.csv',inp,dest]
    report=dict(scope='Actual outlet flat core plus fixed-reference shoulder quadrature; exact closure and near-seal spatial controls. No moving-wall PDE or speech claim.',profiles=len(profiles),fixtures=len(cases),frequency_Hz=[20,200000],maximum_relative_errors_by_order=errors,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/lip_outlet_viscous_reference.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(errors,indent=2))
if __name__=='__main__':main()
