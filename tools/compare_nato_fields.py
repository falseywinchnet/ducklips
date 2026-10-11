#!/usr/bin/env python3
"""Compare complete enhanced fields at identical native time/frequency spacing."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.io import wavfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).parent/'enhanced_specgram'))
from display import enhanced_rgba
p=argparse.ArgumentParser();p.add_argument('reference',type=Path);p.add_argument('synthetic',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--gain',type=float,default=1);p.add_argument('--label',default='Current sheet engine');a=p.parse_args()
r=np.load(a.reference);s=np.load(a.synthetic)
for key in ['sample_rate','hop','target_n_fft']:assert int(r[key])==int(s[key]),f'Cannot compare unlike {key}'
assert np.array_equal(r['hz'],s['hz']),'Frequency axes differ'
ref=.17808298734486477;rows=len(r['hz']);widths=[len(r['seconds']),len(s['seconds'])];W=sum(widths)+190;H=rows+105
fig=plt.figure(figsize=(W/100,H/100),dpi=100,facecolor='#111814')
for z,width,left,title,gain in [(r,widths[0],65,'Human reference',1),(s,widths[1],130+widths[0],a.label,a.gain)]:
 ax=fig.add_axes([left/W,55/H,width/W,rows/H]);data=enhanced_rgba(z['enhanced']*gain,ref,70)
 hop=float(z['hop'])/float(z['sample_rate']);step=z['hz'][1]-z['hz'][0]
 ax.imshow(data,interpolation='nearest',aspect='auto',extent=[float(z['seconds'][0])-hop/2,float(z['seconds'][-1])+hop/2,-step/2,float(z['hz'][-1])+step/2])
 ax.set_title(title,color='#e8e9df',fontsize=10);ax.tick_params(colors='#becab8',labelsize=8);ax.set_xlabel('Seconds · no time warping',color='#becab8',fontsize=8);ax.set_ylabel('Hz',color='#becab8',fontsize=8)
fig.text(.04,.013,f'684 native rows · 1 pixel per time cell · magma / 70 dB · synthetic display gain {a.gain:.6g}',color='#becab8',fontsize=8)
a.out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.out,dpi=100,facecolor=fig.get_facecolor());plt.close(fig)
a.out.with_suffix('.json').write_text(json.dumps(dict(reference=str(a.reference.resolve()),synthetic=str(a.synthetic.resolve()),synthetic_display_gain=a.gain,display_reference=ref,rows=rows,widths=widths,pixel_resampling=False,time_warping=False),indent=2))
