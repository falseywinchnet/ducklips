#!/usr/bin/env python3
"""Independent full-state Radau equations and continuous ODE for coupled losses."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from slit_viscous_reference import spectrum,ROOT
OUT=ROOT/'out/relaxation-graph-reference'
A=np.array([[5/12,-1/12],[3/4,1/4]]);BWEIGHT=np.array([.75,.25])
def fixture(c,g,r,edges,no,bo,nm,bm,state,nodal,series,dt,steps=1):
    n=len(c);lines=[f'{n} {len(g)} {len(no)} {len(bo)} {dt} {steps}']
    items=[np.array(edges).ravel(),no,bo,c,g,r,np.array(nm).ravel(),np.array(bm).ravel(),state,nodal,series]
    lines.extend(' '.join(map(str,np.asarray(x).ravel())) for x in items)
    return '\n'.join(lines)
def system(c,g,r,edges,no,bo,nm,bm):
    n=len(c);m=len(g);a=len(no);b=len(bo);dim=n+m+a+b;matrix=np.zeros((dim,dim));incidence=np.zeros((n,m))
    for i,(u,v) in enumerate(edges):
        incidence[u,i]=1
        if v<n:incidence[v,i]=-1
    matrix[:n,n:n+m]=-incidence*g;matrix[n:n+m,:n]=incidence.T/c
    matrix[n:n+m,n:n+m]=-np.diag(r*g)
    for j,(owner,(k,rate)) in enumerate(zip(no,nm)):
        z=n+m+j;matrix[owner,owner]-=rate*k/c[owner];matrix[owner,z]+=rate
        matrix[z,owner]=rate*k/c[owner];matrix[z,z]=-rate
    for j,(owner,(k,rate)) in enumerate(zip(bo,bm)):
        z=n+m+a+j;p=n+owner;matrix[p,p]-=rate*k*g[owner];matrix[p,z]+=rate
        matrix[z,p]=rate*k*g[owner];matrix[z,z]=-rate
    metric=np.sqrt(np.r_[1/c,g,1/np.array(nm)[:,0],1/np.array(bm)[:,0]])
    return matrix*metric[:,None]/metric[None,:],metric,incidence

def main():
    OUT.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(14102026);lines=[];references=[]
    edges=[(0,1),(1,2),(2,0),(2,3)];no=[0,2];bo=[1,1,3,3]
    for i in range(80):
        c=10**rng.uniform(-12,-10,3);g=10**rng.uniform(-4,-2,4);r=10**rng.uniform(2,6,4)
        nm=np.array([(c[o]*rng.uniform(.01,.4),10**rng.uniform(1,8)) for o in no]);bm=np.array([(rng.uniform(.01,.2)/g[o],10**rng.uniform(1,8)) for o in bo])
        matrix,metric,incidence=system(c,g,r,edges,no,bo,nm,bm);state=rng.normal(size=len(metric))*1e-5/metric
        nodal=rng.normal(size=3)*1e-5;series=rng.normal(size=4)*100;dt=10**rng.uniform(-8,-4)
        force=np.r_[nodal,series,np.zeros(6)]*metric;y=state*metric
        stages=np.linalg.solve(np.eye(2*len(y))-dt*np.kron(A,matrix),np.tile(y,2)+dt*np.kron(A.sum(axis=1),force)).reshape(2,-1)
        final=stages[1]/metric;physical=stages/metric
        pressure=BWEIGHT@physical[:,:3]/c;flow=(BWEIGHT@physical[:,3:7])*g
        lines.append(fixture(c,g,r,edges,no,bo,nm,bm,state,nodal,series,dt));references.append((final,pressure,flow,metric,dt))
    source=OUT/'random.input';source.write_text('\n'.join(lines)+'\n');binary=ROOT/'build/ducklips_relaxation_graph_probe'
    result=subprocess.run([str(binary),'--fixture'],input=source.read_text(),text=True,capture_output=True);assert result.returncode==0,result.stderr
    dest=OUT/'random.txt';dest.write_text(result.stdout);data=np.loadtxt(dest);worst=pressure_error=flow_error=ledger=0.
    for row,(final,pressure,flow,metric,dt) in zip(data,references,strict=True):
        dim=len(metric);worst=max(worst,float(max(abs((row[:dim]-final)*metric))/1e-5))
        pressure_error=max(pressure_error,float(max(abs(row[dim:dim+3]-pressure))));flow_error=max(flow_error,float(max(abs(row[dim+3:dim+7]-flow))))
        ledger=max(ledger,abs(row[-2]));assert row[-5]>=0 and row[-4]>=0
    assert worst<1e-8 and pressure_error<1e-5 and flow_error<1e-12 and ledger<1e-17,(worst,pressure_error,flow_error,ledger)
    # Continuous trajectory independent of the native stage elimination.
    c=np.array([1e-11,1.4e-11,8e-12]);g=np.array([.002,.003,.001,.0025]);r=np.array([0.,1e5,0,1e4])
    nm=np.array([[4e-12,1e5],[3e-12,2e7]]);bm=np.array([[12.,1e4],[18.,3e7],[100.,1e5],[50.,2e6]])
    matrix,metric,_=system(c,g,r,edges,no,bo,nm,bm);state=np.r_[[1e-9,-2e-9,5e-10],[.01,-.02,.015,.005],[2e-10,-1e-10],[.001,-.002,.003,.002]]
    nodal=np.array([0,1e-5,0]);series=np.zeros(4);force=np.r_[nodal,series,np.zeros(6)]*metric;duration=.0003
    augmented=np.zeros((len(metric)+1,len(metric)+1));augmented[:-1,:-1]=matrix;augmented[:-1,-1]=force
    exact=(expm(augmented*duration)@np.r_[state*metric,1.])[:-1]
    ode=solve_ivp(lambda t,y:matrix@y+force,(0,duration),state*metric,method='Radau',jac=matrix,rtol=2e-11,atol=1e-15)
    assert ode.success;oracle_error=float(max(abs(exact-ode.y[:,-1])));assert oracle_error<1e-12
    refinements=[]
    for steps in [16,32,64,128,256]:
        inp=OUT/f'trajectory-{steps}.input';inp.write_text(fixture(c,g,r,edges,no,bo,nm,bm,state,nodal,series,duration/steps,steps)+'\n')
        run=subprocess.run([str(binary),'--fixture'],input=inp.read_text(),text=True,capture_output=True);assert run.returncode==0,run.stderr
        target=inp.with_suffix('.txt');target.write_text(run.stdout);row=np.fromstring(run.stdout,sep=' ')
        refinements.append(dict(steps=steps,maximum_energy_coordinate_error=float(max(abs(row[:len(metric)]*metric-exact)))))
    # Same constant-pressure slit that exposed the separated update's defect.
    weights,poles=spectrum();inertia=100.;gap=1e-5;pressure=700.;nu=1.86e-5/1.14;steady=pressure/(3*inertia*nu/gap**2);dc=3*inertia*nu/gap**2
    forced=[]
    for initialization in ['rest','steady']:
        for stiffness in [.01,.1,1.,10.,100.,10000.]:
            dt=stiffness*gap**2/(3*nu);bm=np.column_stack([inertia*weights/poles,nu/gap**2*poles]);state=np.zeros(26)
            if initialization=='steady':state[1]=steady*inertia;state[2:]=bm[:,0]*steady
            inp=OUT/f'forced-{initialization}-{stiffness}.input';inp.write_text(fixture([1e20],[1/inertia],[dc],[(0,1)],[],[0]*24,[],bm,state,[0],[pressure],dt,8000)+'\n')
            run=subprocess.run([str(binary),'--fixture'],input=inp.read_text(),text=True,capture_output=True);assert run.returncode==0,run.stderr
            inp.with_suffix('.txt').write_text(run.stdout);row=np.fromstring(run.stdout,sep=' ');ratio=row[27]/steady
            assert abs(ratio-1)<1e-11 and abs(row[26])<1e-20,(initialization,stiffness,ratio,row[26])
            forced.append(dict(initialization=initialization,dc_decay_times_step=stiffness,average_flow_ratio=ratio,reservoir_pressure_deviation_Pa=float(row[26])))
    paths=[Path(__file__),ROOT/'src/relaxation_graph.cpp',ROOT/'include/ducklips/relaxation_graph.hpp',ROOT/'tools/relaxation_graph_probe.cpp',ROOT/'tools/slit_viscous_reference.py',binary,*sorted(OUT.iterdir())]
    report=dict(scope='Coupled fixed-metric pressure/transport/relaxation update. Independent full-state Radau equations, continuous matrix exponential and stiff ODE, and large-compliance constant-pressure reservoir controls. No speech claim.',random_cases=80,maximum_scaled_state_error=worst,maximum_pressure_error_Pa=pressure_error,maximum_flow_error_m3_s=flow_error,maximum_energy_residual_J=ledger,continuous_oracle_agreement=oracle_error,refinements=refinements,forced_flow=forced,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/relaxation_graph_reference.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='hashes'},indent=2))
if __name__=='__main__':main()
