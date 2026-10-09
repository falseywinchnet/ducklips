#!/usr/bin/env python3
"""Audit exponential-relaxation contour kinetics and contact cycles."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import mpmath as mp
from contour_sealed_reference import metrics
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/contour-relaxation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run_json(command,path):
    result=subprocess.run(list(map(str,command)),capture_output=True,text=True);assert result.returncode==0,result.stderr
    path.write_text(result.stdout);return json.loads(result.stdout)
def state(s,graph=False):
    if graph:return np.r_[s['mass'],s['position'],s['velocity'],s['passage_momentum'],np.array(s['charge'])/np.array([2e-10,3e-10,1e-10]),np.array([.002,.003,.004,.002])*np.array(s['network_momentum'])]
    return np.r_[s['mass'],s['position'],s['velocity'],s['momentum']]

def kernel():
    rng=np.random.default_rng(9517);rows=[]
    for n in range(144):
        a=10**rng.uniform(-7,-4,2);b=a*rng.uniform(.8,1.2,2)
        if 80<=n<128:
            cell=n%2
            if cell==0:a[0]=-.3*a[1];b[0]=-.3*b[1]
            else:a[1]=-.05*a[0];b[1]=-.05*b[0]
        if n>=128:b[n%2]=-.05*b[1-n%2]
        m0=np.array(metrics(np.r_[1e-9,1e-9,a])[0])*rng.uniform(.8,1.2,2);m1=np.array(metrics(np.r_[1e-9,1e-9,b])[0])*rng.uniform(.8,1.2,2)
        if n<40:m1=m0*(1+rng.uniform(-1e-10,1e-10,2));b=a*(1+rng.uniform(-1e-10,1e-10,2))
        p0=rng.uniform(-.003,.003,3);p1=rng.uniform(-.003,.003,3)
        for i in range(2):
            if a[i]<=0:p0[0 if i==0 else 2]=0
            if b[i]<=0:p1[0 if i==0 else 2]=0
        row=[*m0,*a,*p0,*m1,*b,*p1];rows.append(row)
        if n>=128:rows.append(row[7:]+row[:7])
    for n,row in enumerate(rows):row.extend([0.,0.,0.] if n==0 else [1.,1.,1.] if n==1 else rng.uniform(0,1,3).tolist())
    inp=OUT/'dissipation-input.txt';inp.write_text(''.join(' '.join(map(str,r))+'\n' for r in rows));binary=ROOT/'build/ducklips_contour_dissipation_probe';r=subprocess.run([str(binary),'--weighted'],input=inp.read_text(),capture_output=True,text=True);assert r.returncode==0,r.stderr
    output=OUT/'dissipation-output.jsonl';output.write_text(r.stdout);native=[json.loads(s) for s in r.stdout.splitlines()];assert len(native)==len(rows);mp.mp.dps=80
    maxima=dict(effort_scaled=0.,flow_scaled=0.,numerical_loss_scaled=0.,chain_scaled=0.)
    def energy(x,p):return sum(pi*pi*g/2 for pi,g in zip(p,metrics(x,mp.log)[1]))
    for row,s in zip(rows,native):
        row=list(map(mp.mpf,row));x=row[:4];p0=row[4:7];y=row[7:11];p1=row[11:14];weight=row[14:];g0=metrics(x,mp.log)[1];g1=metrics(y,mp.log)[1];e0=energy(x,p0);e1=energy(y,p1);scale=max(e0+e1,mp.mpf('1e-60'));expected=[mp.mpf(0)]*4
        def shape_energy(z):return sum(((1-w)*(a*a+b*b)/4+w*a*a/2)*g for a,b,w,g in zip(p0,p1,weight,metrics(z,mp.log)[1]))
        for order in [range(4),range(3,-1,-1)]:
            current=x.copy()
            for k in order:
                nxt=current.copy();nxt[k]=y[k];delta=nxt[k]-current[k]
                d=mp.diff(lambda z:shape_energy([z if j==k else current[j] for j in range(4)]),current[k])/2 if delta==0 else (shape_energy(nxt)-shape_energy(current))/(2*delta)
                expected[k]+=d;current=nxt
        for k in range(4):maxima['effort_scaled']=max(maxima['effort_scaled'],float(abs(s['effort'][k]-expected[k])*max(abs(x[k]),abs(y[k]))/scale))
        for j in range(3):
            target=(1-weight[j])*(g0[j]+g1[j])*(p0[j]+p1[j])/4+weight[j]*g1[j]*p1[j];maxima['flow_scaled']=max(maxima['flow_scaled'],float(abs(s['flow'][j]-target)/max(abs(target),mp.mpf('1e-60'))))
        loss=sum(w*g*(b-a)**2/2 for w,g,a,b in zip(weight,g1,p0,p1));maxima['numerical_loss_scaled']=max(maxima['numerical_loss_scaled'],float(abs(s['numerical_loss']-loss)/scale));maxima['chain_scaled']=max(maxima['chain_scaled'],float(abs(s['chain_residual'])/scale));assert s['numerical_loss']>=0
    assert max(maxima.values())<1e-9,maxima
    rates=np.r_[0.,np.geomspace(1e-12,1e12,97)];weight_input=OUT/'weights-input.txt';weight_input.write_text(''.join(str(x)+'\n' for x in rates))
    result=subprocess.run([str(binary),'--weights'],input=weight_input.read_text(),capture_output=True,text=True);assert result.returncode==0;weight_output=OUT/'weights-output.txt';weight_output.write_text(result.stdout)
    weights=list(map(float,result.stdout.splitlines()));assert len(weights)==len(rates);max_weight=0.;max_decay=0.
    for x,beta in zip(rates,weights):
        target=mp.mpf(0) if x==0 else mp.coth(mp.mpf(x)/2)-2/mp.mpf(x);max_weight=max(max_weight,float(abs(beta-target)))
        decay=(1-(1-beta)*x/2)/(1+(1+beta)*x/2);max_decay=max(max_decay,float(abs(decay-mp.exp(-mp.mpf(x)))));assert 0<=beta<=1
    assert max_weight<1e-13 and max_decay<1e-14
    maxima.update(relaxation_weight_absolute=max_weight,fixed_passage_decay_absolute=max_decay)
    return dict(cases=len(rows),digits=80,maxima=maxima,input=str(inp.relative_to(ROOT)),input_sha256=sha(inp),output=str(output.relative_to(ROOT)),output_sha256=sha(output))

def open_motion():
    commit='5bc0d1ace9273b158346921fd1a652edbb8c508f';ref=ROOT/'reports/contour_glottal_reference.json';original=subprocess.check_output(['git','show',commit+':reports/contour_glottal_reference.json'],cwd=ROOT);assert original==ref.read_bytes()
    for path in ['tools/contour_glottal_reference.py','src/gas_storage.cpp','include/ducklips/physics.hpp','include/ducklips/sources.hpp']:
        assert subprocess.check_output(['git','show',commit+':'+path],cwd=ROOT)==(ROOT/path).read_bytes(),path
    report=json.loads(original);scale=np.r_[[1e-8]*2,[1e-4]*4,[1.]*4,[.05]*3,[700.]*3,[1e-4]*4];runs=[];acceptance=[]
    for mode in ['forward','reverse','relaxation','rest_motion']:
        oracle=np.array(next(r['final_state'] for r in report['oracles'] if r['mode']==mode));previous=None
        for steps in [1920,3840,7680,15360]:
            command=[ROOT/'build/ducklips_contour_glottal_probe','--trajectory',str(steps),'.005',mode,'--relaxation'];path=OUT/f'open-{mode}-{steps}.json';s=run_json(command,path);error=abs(state(s,True)-oracle);norm=float(max(error/scale));runs.append(dict(mode=mode,steps=steps,normalized_error=norm,refinement_ratio=None if previous is None else previous/norm,numerical_loss_J=s['numerical_loss_J'],physical_loss_J=s['physical_loss_J'],reservoir_work_J=s['reservoir_work_J'],maximum_work_residual_J=s['maximum_work_balance'],output=str(path.relative_to(ROOT)),output_sha256=sha(path),command=list(map(str,command))));previous=norm
        acceptance.append(dict(mode=mode,second_order=all(3.5<r['refinement_ratio']<4.5 for r in runs[-2:]),endpoint_accuracy=previous<1e-4));print('open',json.dumps(runs[-1]),flush=True)
    return dict(scope='fresh fitted-relaxation runs against reused independent continuous open endpoints',oracle_commit=commit,oracle_report=str(ref.relative_to(ROOT)),oracle_report_sha256=sha(ref),runs=runs,acceptance=acceptance)

def cycles():
    scale=np.r_[[1e-8]*2,[1e-4]*4,[1.]*4,[.05]*3];runs=[];acceptance=[]
    for mode in ['lower','upper']:
        previous_state=None;previous_difference=None;previous_loss=None
        for steps in [500,1000,2000,4000,8000,16000]:
            command=[ROOT/'build/ducklips_contour_cycle_probe',str(steps),mode,'--relaxation'];path=OUT/f'cycle-{mode}-{steps}.json';s=run_json(command,path);y=state(s);difference=None if previous_state is None else float(max(abs(y-previous_state)/scale));ratio=None if previous_difference is None else previous_difference/difference
            receipt={key:s[key] for key in ['closures','reopenings','first_closed_endpoint_s','first_reopened_endpoint_s','initial_energy_J','numerical_loss_J','physical_loss_J','maximum_energy_before_sealing_J','maximum_work_residual_J','maximum_mass_residual_m3','minimum_area_m2','maximum_density_ratio']}
            receipt.update(mode=mode,steps=steps,normalized_successive_difference=difference,successive_difference_ratio=ratio,numerical_loss_ratio=None if previous_loss is None else previous_loss/s['numerical_loss_J'],output=str(path.relative_to(ROOT)),output_sha256=sha(path),command=list(map(str,command)));runs.append(receipt)
            assert s['closures']==1 and s['reopenings']==1;previous_state=y;previous_difference=difference;previous_loss=s['numerical_loss_J']
        acceptance.append(dict(mode=mode,successive_differences_decrease=all(r['successive_difference_ratio']>1 for r in runs[-2:]),final_successive_difference=difference<1e-4,numerical_loss_decreases=all(r['numerical_loss_ratio']>1 for r in runs[-2:]),numerical_loss_fraction=runs[-1]['numerical_loss_J']/runs[-1]['initial_energy_J']<1e-4));print('cycle',json.dumps(runs[-1]),flush=True)
    return dict(scope='native temporal self-refinement across closure and reopening; no independent full-cycle continuous oracle',runs=runs,acceptance=acceptance)

def main():
    OUT.mkdir(parents=True,exist_ok=True);constitutive=kernel();print('kernel',constitutive,flush=True);motion=open_motion();cycle=cycles()
    sources=['include/ducklips/contour_kinetics.hpp','src/contour_kinetics.cpp','include/ducklips/contour_glottal_storage.hpp','src/contour_glottal_storage.cpp','include/ducklips/fold_contact.hpp','src/fold_contact.cpp','src/glottal_contour.cpp','src/gas_storage.cpp','src/moving_graph.cpp','tools/contour_glottal_probe.cpp','tools/contour_cycle_probe.cpp','tools/contour_dissipation_probe.cpp','tools/contour_relaxation_reference.py','tools/contour_sealed_reference.py']
    bins=['build/ducklips_contour_glottal_probe','build/ducklips_contour_cycle_probe','build/ducklips_contour_dissipation_probe'];report=dict(scope='exponentially fitted dissipative kinetic update with cap closure/reopening; positive chamber storage only',kernel=constitutive,open_motion=motion,cycles=cycle,hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources+bins)]},versions=dict(numpy=np.__version__,mpmath=mp.__version__))
    (ROOT/'reports/contour_relaxation_reference.json').write_text(json.dumps(report,indent=2)+'\n')
    assert all(a['second_order'] and a['endpoint_accuracy'] for a in motion['acceptance']),motion['acceptance'];assert all(all(v for k,v in a.items() if k!='mode') for a in cycle['acceptance']),cycle['acceptance']
if __name__=='__main__':main()
