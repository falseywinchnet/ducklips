#!/usr/bin/env python3
"""Closed-port and reverse-flow checks with an attached kinetic column."""
import json
import subprocess
from pathlib import Path
import numpy as np
from graph_continuous_glottis_reference import ROOT, EDGES, C, G, NODAL, SERIES, PORT, INITIAL, XINITIAL
from planar_inertive_motion_reference import sha

out=ROOT/'out/continuous-outlet-column-controls';out.mkdir(exist_ok=True)
records=[]
for build in ['build','.build/asan']:
    binary=ROOT/build/'ducklips_graph_continuous_glottis_probe'
    for name in ['sealed','reverse']:
        y=INITIAL.copy();y[:2]=0.;y[7:]=0.
        if name=='sealed':y[2:4]=[3e-6,-1e-7];y[4:7]=0.
        else:y[4:7]=-1e-4
        lines=[f'3 4 1 {PORT}',*[f'{a} {b}' for a,b in EDGES],' '.join(map(str,y)), '1','80',' '.join(map(str,XINITIAL))]
        lines.append(' '.join(map(str,np.r_[1e-7,0.,2.,160. if name=='sealed' else 80.,C,C,G,G,NODAL,SERIES])))
        stem=out/(('asan-' if build.startswith('.') else '')+name)
        fixture='\n'.join(lines)+'\n';stem.with_suffix('.input').write_text(fixture)
        result=subprocess.run([str(binary),'--column-fixture'],input=fixture,capture_output=True,text=True)
        stem.with_suffix('.txt').write_text(result.stdout);stem.with_suffix('.stderr').write_text(result.stderr)
        assert result.returncode==0,result.stderr
        data=np.fromstring(result.stdout,sep=' ');assert len(data)==44
        if name=='sealed':assert data[36]==0 and data[33]==0 and data[42]==0
        else:assert data[36]<0 and data[30]>0
        assert abs(data[32])<1e-13 and abs(data[35])<1e-20
        records.append(dict(name=name,build=build,flow_m3_s=float(data[36]),outlet_loss_J=float(data[30]),geometry_work_J=float(data[42]),source_residual_J=float(data[32]),
            hashes={str(p.relative_to(ROOT)):sha(p) for p in [binary,*[stem.with_suffix(e) for e in ['.input','.txt','.stderr']]]}))
report=dict(scope='Native and ASAN/UBSAN sealed-column geometry change and reverse-flow controls',cases=records,
    hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'tools/graph_continuous_glottis_probe.cpp']})
(out/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(records,indent=2))
