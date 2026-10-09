#!/usr/bin/env python3
"""Independent quadrature of signed fold contours and their open air regions."""
import hashlib
import json
from pathlib import Path
import subprocess
import warnings
import numpy as np
import scipy
from scipy.integrate import quad, IntegrationWarning
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    out=ROOT/'out/glottal-contact';out.mkdir(parents=True,exist_ok=True);binary=ROOT/'build/ducklips_glottal_contour_probe'
    rng=np.random.default_rng(930);fixtures=[]
    for i in range(120):
        h=10**rng.uniform(-4,-2,2);a=rng.uniform(-1,1,2)*10**rng.uniform(-7,-4);fixtures.append([*h,*a])
    for b in [3e-6,0.,3e-12,-3e-12,-.6e-6+3e-12,-.6e-6-3e-12,-1e-6]:fixtures.append([.0025,.0005,3e-6,b])
    fixtures.extend([[.0025,.0005,0.,3e-6],[.0025,.0005,0.,0.],[.0025,.0005,-1e-6,-2e-6]])
    motion=ROOT/'out/glottal-contact/motion-1000000.csv';poses=np.genfromtxt(motion,delimiter=',',names=True)
    fixtures.extend([[.0025,.0005,float(a),float(b)] for a,b in zip(poses['lower_area_m2'],poses['upper_area_m2'])])
    invalid=[[0.,.0005,3e-6,1e-6],[.0025,-.0005,3e-6,1e-6]]
    inp=out/'contour-input.txt';inp.write_text(''.join(' '.join(map(str,row))+'\n' for row in fixtures+invalid))
    run=subprocess.run([str(binary),'--batch'],input=inp.read_text(),capture_output=True,text=True);assert run.returncode==0,run.stderr
    output=out/'contour-output.jsonl';output.write_text(run.stdout);native=[json.loads(s) for s in run.stdout.splitlines()]
    assert len(native)==len(fixtures)+len(invalid) and all(not s['accepted'] for s in native[len(fixtures):])
    maxima=dict(volume_scaled=0.,centroid_scaled=0.,gradient_scaled=0.,inverse_area_relative=0.,inverse_cube_relative=0.,aperture_scaled=0.)
    # Work in normalized longitudinal/area coordinates; quadrature operates
    # on the interpolated profile, not the native clipped-trapezoid formula.
    with warnings.catch_warnings():
        warnings.simplefilter('error',IntegrationWarning)
        for row,actual in zip(fixtures,native):
            assert actual['accepted'],row
            h0,h1,a,b=row;H=h0+h1;S=max(abs(a),abs(b),1e-30);seam=h0/H
            nodes=np.array([0,h0/2,h0+h1/2,H])/H;values=np.array([a,a,b,b])/S
            def raw(t):return float(np.interp(t,nodes,values))
            knots=list(nodes[1:-1])+[seam]
            if a*b<0:knots.append(brentq(raw,float(nodes[1]),float(nodes[2]),xtol=5e-16))
            def integrate(fn,l,u,extra=()):
                if u<=l:return 0.
                points=sorted(set([l,u,*[x for x in [*knots,*extra] if l<x<u]]))
                return sum(quad(fn,x,y,epsabs=2e-13,epsrel=2e-12)[0] for x,y in zip(points[:-1],points[1:]))
            chambers=[]
            for index,(l,u) in enumerate([(0.,seam),(seam,1.)]):
                volume=integrate(lambda t:max(0.,raw(t)),l,u);moment=integrate(lambda t:t*max(0.,raw(t)),l,u)
                centroid=moment/volume if volume>0 else .5*(l+u)
                def gradient(t,k):return float(np.interp(t,nodes,[1.,1.,0.,0.] if k==0 else [0.,0.,1.,1.])) if raw(t)>0 else 0.
                grad=[integrate(lambda t:gradient(t,k),l,u) for k in range(2)]
                c=actual['chambers'][index];maxima['volume_scaled']=max(maxima['volume_scaled'],abs(c['volume']/(H*S)-volume))
                if volume>1e-18:maxima['centroid_scaled']=max(maxima['centroid_scaled'],abs(c['centroid']/H-centroid))
                maxima['gradient_scaled']=max(maxima['gradient_scaled'],max(abs(np.array(c['derivative'])/H-grad)))
                chambers.append((volume,centroid))
            maxima['aperture_scaled']=max(maxima['aperture_scaled'],max(abs(np.array(actual['aperture'])/S-np.maximum(0,[raw(0),raw(seam),raw(1)]))))
            paths=[(0.,chambers[0][1]),(chambers[0][1],chambers[1][1]),(chambers[1][1],1.)]
            for index,(l,u) in enumerate(paths):
                has_air=chambers[0][0]>0 if index==0 else chambers[1][0]>0 if index==2 else all(c[0]>0 for c in chambers)
                opened=has_air and min(raw(l),raw(u))>0
                p=actual['passages'][index];assert p['open']==opened,(row,index,p,opened)
                if not opened:assert p['inverse_area']==[0.,0.] and p['inverse_area_cubed']==[0.,0.];continue
                for cell,(cl,cu) in enumerate([(0.,seam),(seam,1.)]):
                    x=max(l,cl);y=min(u,cu)
                    if y<=x:
                        assert p['inverse_area'][cell]==0 and p['inverse_area_cubed'][cell]==0,(row,index,cell,l,u,x,y,p);continue
                    minimum=min(raw(x),raw(y))
                    for power,key,scale_key in [(1,'inverse_area','inverse_area_relative'),(3,'inverse_area_cubed','inverse_cube_relative')]:
                        # Resolve thin endpoint layers by geometric area
                        # subdivisions, then use ordinary adaptive quadrature.
                        extra=[];segments=sorted(set([x,y,*[t for t in knots if x<t<y]]))
                        for sl,su in zip(segments[:-1],segments[1:]):
                            aa,ab=raw(sl),raw(su)
                            if aa==ab:continue
                            for target in np.geomspace(min(aa,ab),max(aa,ab),max(2,int(abs(np.log10(ab/aa)))+2))[1:-1]:
                                extra.append(sl+(su-sl)*(target-aa)/(ab-aa))
                        try:expected=integrate(lambda t:(minimum/raw(t))**power,x,y,extra)
                        except IntegrationWarning as e:raise RuntimeError((row,index,cell,power,x,y,extra)) from e
                        measured=p[key][cell]*(S*minimum)**power/H
                        maxima[scale_key]=max(maxima[scale_key],abs(measured-expected)/max(expected,1e-300))
    assert maxima['volume_scaled']<2e-12 and maxima['centroid_scaled']<1e-9 and maxima['gradient_scaled']<1e-9,maxima
    assert maxima['inverse_area_relative']<1e-8 and maxima['inverse_cube_relative']<1e-8 and maxima['aperture_scaled']<1e-12,maxima
    checks=subprocess.run([str(binary)],capture_output=True,text=True);assert checks.returncode==0,checks.stderr
    log=out/'contour-native.log';log.write_text(checks.stdout)
    sources=['src/glottal_contour.cpp','include/ducklips/glottal_contour.hpp','tools/glottal_contour_probe.cpp','tools/glottal_contour_reference.py']
    receipt=dict(scope='piecewise-affine signed medial profile; numerical quadrature of storage, gradients and passage integrals; geometry only',
        valid_cases=len(fixtures),invalid_cases=len(invalid),motion_poses=len(poses),maxima=maxima,
        input=str(inp.relative_to(ROOT)),input_sha256=sha(inp),output=str(output.relative_to(ROOT)),output_sha256=sha(output),
        motion_input=str(motion.relative_to(ROOT)),motion_input_sha256=sha(motion),native_checks=str(log.relative_to(ROOT)),native_checks_sha256=sha(log),
        hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),binary]},versions=dict(numpy=np.__version__,scipy=scipy.__version__))
    (ROOT/'reports/glottal_contour_reference.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
