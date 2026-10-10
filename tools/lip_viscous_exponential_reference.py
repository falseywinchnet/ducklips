#!/usr/bin/env python3
"""Check the full-matrix audit propagator at representative lip stiffnesses."""
import hashlib,json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.linalg import expm
from planar_continuous_reference import decay_exponential
from slit_viscous_reference import spectrum
ROOT=Path(__file__).resolve().parents[1]
def main():
    w,p=spectrum();mp.mp.dps=50;records=[];paths=[]
    for case in ['closing','opening']:
        path=ROOT/f'out/planar-lip-viscous-16/{case}-1-512.json';paths.append(path)
        data=json.loads(path.read_text());items=[]
        for step in data['steps']:
            for site in step['viscous_sites']:
                time=step['dt']/2*(1.86e-5/1.14)/site['half_gap']**2
                items.append((time,site['segment_inertance']*site['branch_conductance']))
        items.sort()
        for label,index in [('least_stiff',0),('median',len(items)//2),('most_stiff',-1)]:
            time,ratio=items[index]
            matrix=np.zeros((25,25));matrix[0,0]=-ratio*(3+sum(w))
            matrix[0,1:]=matrix[1:,0]=np.sqrt(ratio*w*p);matrix[1:,1:]=-np.diag(p)
            matrix*=time
            reference=np.array(mp.expm(mp.matrix(matrix.tolist())).tolist(),dtype=float)
            ordinary=expm(matrix);increment=decay_exponential(matrix)
            error=float(np.max(abs(increment-reference)))
            assert error<2e-13,(case,label,error)
            records.append(dict(case=case,selection=label,time=time,inertia_fraction=ratio,
                matrix_norm_1=float(np.linalg.norm(matrix,1)),
                scipy_expm_max_absolute_error=float(np.max(abs(ordinary-reference))),
                increment_exponential_max_absolute_error=error))
    paths += [Path(__file__),ROOT/'tools/planar_continuous_reference.py',ROOT/'tools/slit_viscous_reference.py']
    report=dict(scope='Six selected actual matrices, compared with 50-digit mpmath full matrix exponential; not a bound on every trajectory or matrix.',precision_digits=50,cases=records,
        hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths})
    (ROOT/'reports/lip_viscous_exponential_reference.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(records,indent=2))
if __name__=='__main__':main()
