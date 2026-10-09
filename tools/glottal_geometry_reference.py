#!/usr/bin/env python3
"""Independent kinetic-energy and continuous articulated-passage checks."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import scipy
import mpmath as mp
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()



def kinetic_audit(binary,out):
    rng=np.random.default_rng(9031);fixtures=[]
    for i in range(240):
        m0=10**rng.uniform(-12,-7,2)
        m1=m0*(1+rng.uniform(-1e-12,1e-12,2)) if i<80 else m0*10**rng.uniform(-1,1,2)
        p0=rng.uniform(-1,1,3)*10**rng.uniform(-5,-1)
        p1=rng.uniform(-1,1,3)*10**rng.uniform(-5,-1)
        fixtures.append(np.r_[m0,m1,p0,p1].tolist())
    fixtures.append([7.5e-9,1.5e-9,7.5e-9,1.5e-9,0,0,0,0,0,0])
    invalid=[[0,1e-9,1e-9,1e-9,0,0,0,0,0,0],[-1e-9,1e-9,1e-9,1e-9,0,0,0,0,0,0]]
    inp=out/'kinetic-input.txt';inp.write_text(''.join(' '.join(map(str,v))+'\n' for v in fixtures+invalid))
    r=subprocess.run([str(binary),'--pairs'],input=inp.read_text(),capture_output=True,text=True);assert r.returncode==0,r.stderr
    output=out/'kinetic-output.jsonl';output.write_text(r.stdout);native=[json.loads(s) for s in r.stdout.splitlines()]
    assert len(native)==len(fixtures)+len(invalid) and all(not s['accepted'] for s in native[len(fixtures):])
    mp.mp.dps=100;rho=mp.mpf(1.14);h=list(map(mp.mpf,[.0025,.0005]))
    def g(m):
        half=[rho*h[i]**2/(2*m[i]) for i in range(2)]
        return [1/half[0],1/sum(half),1/half[1]]
    def energy(m,p):return sum(p[j]**2*g(m)[j]/2 for j in range(3))
    maxima=dict(energy_relative=0.,instant_effort_relative=0.,transfer_effort_relative=0.,flow_scaled=0.,chain_scaled=0.)
    for row,s in zip(fixtures,native):
        assert s['accepted'];row=list(map(mp.mpf,row));m0=row[:2];m1=row[2:4];p0=row[4:7];p1=row[7:]
        a=energy(m0,p0);b=energy(m1,p1);scale=max(a,b,mp.mpf('1e-50'))
        maxima['energy_relative']=max(maxima['energy_relative'],float(max(abs(mp.mpf(s['energy0'])-a),abs(mp.mpf(s['energy1'])-b))/scale))
        for i in range(2):
            effort=mp.diff(lambda z:energy([z,m0[1]] if i==0 else [m0[0],z],p0),m0[i])
            maxima['instant_effort_relative']=max(maxima['instant_effort_relative'],float(abs(mp.mpf(s['instant_effort0'][i])-effort)/max(abs(effort),mp.mpf('1e-50'))))
        # Independent symmetrized coordinate increments of the whole energy
        # at fixed endpoint momenta; average both momentum endpoints.
        middle0=[m1[0],m0[1]];middle1=[m0[0],m1[1]]
        for i in range(2):
            if m0[i]==m1[i]:
                effort=sum(mp.diff(lambda z:energy([z,m0[1]] if i==0 else [m0[0],z],p),m0[i]) for p in [p0,p1])/2
            else:
                if i==0:effort=sum(energy(middle0,p)-energy(m0,p)+energy(m1,p)-energy(middle1,p) for p in [p0,p1])/(4*(m1[0]-m0[0]))
                else:effort=sum(energy(middle1,p)-energy(m0,p)+energy(m1,p)-energy(middle0,p) for p in [p0,p1])/(4*(m1[1]-m0[1]))
            maxima['transfer_effort_relative']=max(maxima['transfer_effort_relative'],float(abs(mp.mpf(s['mass_effort'][i])-effort)/max(abs(effort),mp.mpf('1e-50'))))
        residual=b-a-sum(mp.mpf(s['mass_effort'][i])*(m1[i]-m0[i]) for i in range(2))-sum(mp.mpf(s['flow'][j])*(p1[j]-p0[j]) for j in range(3))
        maxima['chain_scaled']=max(maxima['chain_scaled'],float(abs(residual)/scale))
        for j in range(3):
            expected=(g(m0)[j]+g(m1)[j])*(p0[j]+p1[j])/4
            flow_scale=max((g(m0)[j]+g(m1)[j])*(abs(p0[j])+abs(p1[j]))/4,mp.mpf('1e-50'))
            maxima['flow_scaled']=max(maxima['flow_scaled'],float(abs(mp.mpf(s['flow'][j])-expected)/flow_scale))
    assert max(maxima.values())<1e-12,maxima
    checks=subprocess.run([str(binary)],capture_output=True,text=True);assert checks.returncode==0,checks.stderr
    log=out/'kinetic-native.log';log.write_text(checks.stdout)
    return dict(valid_cases=len(fixtures),invalid_cases=len(invalid),digits=100,maxima=maxima,input=str(inp.relative_to(ROOT)),input_sha256=sha(inp),output=str(output.relative_to(ROOT)),output_sha256=sha(output),binary=str(binary.relative_to(ROOT)),binary_sha256=sha(binary),native_checks=str(log.relative_to(ROOT)),native_checks_sha256=sha(log))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_graph_glottal_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/glottal-geometry')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/glottal_geometry_reference.json')
    parser.add_argument('--duration',type=float,default=.005)
    parser.add_argument('--steps',type=int,nargs='+',default=[120,240,480,960,1920,3840,7680,15360])
    parser.add_argument('--modes',nargs='+',choices=['forward','reverse','relaxation','rest_motion'],default=['forward','reverse','relaxation','rest_motion'])
    args=parser.parse_args();args.binary=args.binary.resolve();args.output=args.output.resolve();args.report=args.report.resolve();args.output.mkdir(parents=True,exist_ok=True)
    kinetics=kinetic_audit(ROOT/'build/ducklips_glottal_kinetics_probe',args.output)
    bulk=1.14*350**2;gamma=1.4;duration=args.duration;length=.014;thickness=np.array([.0025,.0005]);wall=length*thickness
    tissue_mass=np.array([1.25e-4]*2+[2.5e-5]*2);stiff=np.array([80.]*2+[8.]*2);damping=2*np.array([.1]*2+[.6]*2)*np.sqrt(tissue_mass*stiff)
    rho=1.14;viscosity=1.86e-5
    def metric(m):
        half=.5*rho*thickness**2/m
        return 1/np.array([half[0],sum(half),half[1]])
    c=np.array([2e-10,3e-10,1e-10]);g=np.array([.002,.003,.004,.002]);b=np.array([[1,0,1,1],[-1,1,0,0],[0,-1,-1,0]])
    runs=[];oracles=[];checks_by_mode=[]
    # Coordinates: gas mass, tissue position/velocity, passage momentum,
    # graph pressure/flow. No native finite-change efforts or implicit solver.
    for mode in args.modes:
        reservoir=0. if mode=='relaxation' else -200. if mode=='reverse' else 700.
        reservoir_mu=bulk/(gamma-1)*np.expm1((gamma-1)/gamma*np.log1p(gamma*reservoir/bulk))
        pressure=np.array([700.,300.]) if mode=='relaxation' else np.zeros(2)
        m0=thickness*3e-6*np.power(1+gamma*pressure/bulk,1/gamma)
        y0=np.r_[m0,np.zeros(18)]
        if mode=='relaxation':y0[10:13]=np.array([.001,-.0005,.0004]);y0[13:16]=np.array([2e-8,-1e-8,3e-8])/c
        def derivative(t,y):
            m=y[:2];x=y[2:6];v=y[6:10];momentum=y[10:13];p=y[13:16];flow=y[16:];q=metric(m)*momentum
            rest=3e-6-5e-7*t/duration if mode=='rest_motion' else 3e-6
            volume=thickness*(rest+length*np.array([x[0]+x[1],x[2]+x[3]]));log_r=np.log(m/volume)
            pressure=bulk/gamma*np.expm1(gamma*log_r);mu=bulk/(gamma-1)*np.expm1((gamma-1)*log_r)
            force=np.repeat(pressure*wall,2)-stiff*x-damping*v;coupled=25*(x[:2]-x[2:]);force[:2]-=coupled;force[2:]+=coupled
            # Continuous derivative of kinetic energy, computed from each
            # half-cell's velocity rather than the native rational secants.
            kinetic_head=.25*rho*(thickness/m)**2*np.array([q[0]**2+q[1]**2,q[1]**2+q[2]**2])
            total_mu=mu+kinetic_head
            area=volume/thickness;relative_density=m/volume
            half_r=6*viscosity*thickness*length**2/(area**3*relative_density**2)
            resistance=np.array([half_r[0],sum(half_r),half_r[1]])
            quadratic=np.array([0.,0.,.5*rho*(thickness[1]/m[1])**2])
            pdot=np.array([reservoir_mu-total_mu[0],total_mu[0]-total_mu[1],total_mu[1]-p[2]])-resistance*q-quadratic*q*abs(q)
            drive=-b@flow;drive[2]+=q[2]
            return np.r_[q[0]-q[1],q[1]-q[2],v,force/tissue_mass,pdot,drive/c,g*(b.T@p)]
        options=dict(method='DOP853',rtol=2e-12,atol=np.r_[[1e-23]*2,[1e-16]*4,[1e-10]*4,[1e-12]*3,[1e-8]*3,[1e-15]*4])
        ode=solve_ivp(derivative,(0,duration),y0,max_step=1e-6,**options);assert ode.success
        refined=solve_ivp(derivative,(0,duration),y0,max_step=5e-7,**options);assert refined.success
        difference=abs(ode.y[:,-1]-refined.y[:,-1]);scale=np.r_[[1e-8]*2,[1e-4]*4,[1.]*4,[.05]*3,[700.]*3,[1e-4]*4]
        assert max(difference/scale)<1e-8,difference
        oracles.append(dict(mode=mode,final_state=refined.y[:,-1].tolist(),normalized_refinement_difference=float(max(difference/scale))))
        previous=None
        for steps in args.steps:
            command=[str(args.binary),'--trajectory',str(steps),str(duration),mode,'--articulated']
            r=subprocess.run(command,capture_output=True,text=True);assert r.returncode==0,r.stderr
            path=args.output/f'{mode}-{steps}.json';path.write_text(r.stdout);native=json.loads(r.stdout)
            y=np.r_[native['mass'],native['position'],native['velocity'],native['passage_momentum'],np.array(native['charge'])/c,g*np.array(native['network_momentum'])]
            error=abs(y-refined.y[:,-1]);norm=float(max(error/scale))
            receipt=dict(mode=mode,steps=steps,normalized_error=norm,refinement_ratio=None if previous is None else previous/norm,
                maximum_mass_error_m3=float(max(error[:2])),maximum_position_error_m=float(max(error[2:6])),maximum_velocity_error_m_s=float(max(error[6:10])),
                maximum_passage_momentum_error_Pa_s=float(max(error[10:13])),maximum_network_pressure_error_Pa=float(max(error[13:16])),maximum_network_flow_error_m3_s=float(max(error[16:])),
                maximum_mass_balance_m3=native['maximum_mass_balance'],maximum_work_balance_J=native['maximum_work_balance'],maximum_law_residual_Pa=native['maximum_law_residual'],
                command=command,output=str(path.relative_to(ROOT)),output_sha256=sha(path))
            runs.append(receipt);previous=norm
        print(json.dumps(runs[-1]),flush=True)
        checks_by_mode.append(dict(mode=mode,second_order=all(3.5<r['refinement_ratio']<4.5 for r in runs[-2:]),normalized_endpoint_accuracy=norm<1e-4))
    checks=subprocess.run([str(args.binary),'--articulated'],capture_output=True,text=True);assert checks.returncode==0,checks.stderr
    log=args.output/'native-checks.log';log.write_text(checks.stdout+checks.stderr)
    sources=['include/ducklips/graph_glottal_storage.hpp','src/graph_glottal_storage.cpp','include/ducklips/glottal_storage.hpp','src/glottal_storage.cpp','include/ducklips/gas_storage.hpp','src/gas_storage.cpp','include/ducklips/moving_graph.hpp','src/moving_graph.cpp','include/ducklips/sources.hpp','include/ducklips/physics.hpp','tools/graph_glottal_probe.cpp','tools/glottal_geometry_reference.py','include/ducklips/glottal_kinetics.hpp','src/glottal_kinetics.cpp','tools/glottal_kinetics_probe.cpp']
    args.report.write_text(json.dumps(dict(scope='mass-derived half-passage kinetics, geometric viscous losses and gas/tissue/graph coupling; positive open chambers, no separation/contact or speech',duration_s=duration,kinetics=kinetics,runs=runs,oracles=oracles,acceptance=checks_by_mode,
        native_checks=str(log.relative_to(ROOT)),native_checks_sha256=sha(log),hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),args.binary]},versions=dict(numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__)),indent=2)+'\n')
    assert all(c['second_order'] and c['normalized_endpoint_accuracy'] for c in checks_by_mode),checks_by_mode

if __name__=='__main__':main()
