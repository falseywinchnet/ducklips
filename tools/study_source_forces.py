#!/usr/bin/env python3
"""Observe the rejected source's actual force law without changing its audio.

Uses the frozen current source-stage snapshot, not an older synthesis engine.
No model fitting, waveform prescription, or replacement acoustic analysis.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import numpy as np
from scipy.io import wavfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'out/source-stage-study'
OUT=ROOT/'out/source-force-study'
OUT.mkdir(exist_ok=True)
source=(BASE/'src/sources.cpp').read_text()
source=source.replace('#include <cmath>', '#include <cmath>\n#include <ostream>')
source=source.replace('namespace ducklips {', 'namespace ducklips {\nextern std::ostream* fold_force_trace;\nextern std::size_t fold_force_frame;',1)
marker='            lower_velocity_[side] += step * lower_force / lower_mass[side];'
assert source.count(marker)==1
observer='''            if(side==0 && fold_force_trace){
                *fold_force_trace << fold_force_frame << ',' << substep << ','
                    << lower_position_[side] << ',' << upper_position_[side] << ','
                    << lower_velocity_[side] << ',' << upper_velocity_[side] << ','
                    << lower_drive << ',' << upper_drive << ','
                    << -lower_stiffness[side]*lower_position_[side]-coupling[side]*(lower_position_[side]-upper_position_[side]) << ','
                    << -upper_stiffness[side]*upper_position_[side]-coupling[side]*(upper_position_[side]-lower_position_[side]) << ','
                    << -lower_contact*lower_overlap << ',' << -upper_contact*upper_overlap << ','
                    << -lower_damping*lower_velocity_[side] << ',' << -upper_damping*upper_velocity_[side] << ','
                    << lower_force << ',' << upper_force << ','
                    << strip_gap(rest_area,0,lower_sum) << ',' << strip_gap(rest_area,0,upper_sum) << ','
                    << lower_mass[side] << ',' << upper_mass[side] << '\\n';
            }
'''
(OUT/'sources.cpp').write_text(source.replace(marker,observer+marker))
command=['clang++','-std=c++20','-O3','-I'+str(BASE/'include'),str(ROOT/'tools/source_force_probe.cpp'),str(OUT/'sources.cpp'),str(BASE/'libducklips-base.a'),'-o',str(OUT/'probe')]
subprocess.run(command,check=True)
subprocess.run([str(OUT/'probe'),str(ROOT/'out/visual-round-018/body.body'),str(OUT)],check=True)
raw=OUT/'free-folds-raw.wav'
old=ROOT/'out/source-stage-listening/01-free-folds-raw.wav'
assert raw.read_bytes()==old.read_bytes(), 'Observation changed the source audio'
fs,audio=wavfile.read(raw)
gain=315.98259903736243
wavfile.write(OUT/'free-folds.wav',fs,np.rint(audio*gain*32767).astype(np.int16))
assert (OUT/'free-folds.wav').read_bytes()==(ROOT/'out/source-stage-listening/01-free-folds.wav').read_bytes()
force=np.genfromtxt(OUT/'forces.csv',delimiter=',',names=True)
flow=np.genfromtxt(OUT/'flow.csv',delimiter=',',names=True)
t=(force['frame']+force['substep']/2)/fs
fig,ax=plt.subplots(5,1,figsize=(13,12),sharex=True)
for name in ['lower_gap','upper_gap']:
    ax[0].plot(t,force[name]*1e6,label=name)
ax[0].set_ylabel('opening, mm²')
for name in ['lower_pressure','upper_pressure']:
    ax[1].plot(t,force[name],label=name)
ax[1].set_ylabel('surface Pa')
for name in ['lower_elastic','lower_contact','lower_damping','lower_force']:
    ax[2].plot(t,force[name]*1e3,label=name)
ax[2].plot(t,force['lower_pressure']*.014*.0025*1e3,label='lower aerodynamic')
ax[2].set_ylabel('force, mN')
ax[3].plot(flow['time'],flow['flow']*1e6,label='airflow')
ax[3].set_ylabel('cm³/s')
ax[4].plot(flow['time'],flow['derivative'],label='flow derivative')
ax[4].set_ylabel('m³/s²');ax[4].set_xlabel('s')
for axis in ax:axis.legend(loc='upper right',fontsize=8)
fig.suptitle('Actual forces in the rejected bare source; complete recording')
fig.tight_layout();fig.savefig(OUT/'forces-full.png',dpi=110)
for axis in ax:axis.set_xlim(.5,.526)
fig.suptitle('Supplementary cycle detail; complete enhanced field retained separately')
fig.savefig(OUT/'forces-cycle.png',dpi=110)
# Surface work is an observation, not a naturalness score. Each fold is
# symmetric; multiply one fold by two. Use actual symplectic displacement.
h=1/(2*fs)
dx1=h*(force['lower_v']+h*force['lower_force']/force['lower_mass'])
dx2=h*(force['upper_v']+h*force['upper_force']/force['upper_mass'])
selected=(t>=.2)&(t<1.2)
airwork=2*np.sum((force['lower_pressure']*.014*.0025*dx1+force['upper_pressure']*.014*.0005*dx2)[selected])
dampwork=2*np.sum((force['lower_damping']*dx1+force['upper_damping']*dx2)[selected])
record={'purpose':'Non-invasive source mechanics observation; no new voice candidate',
 'raw_audio_unchanged':True,'raw_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),
 'playback_sha256':hashlib.sha256((OUT/'free-folds.wav').read_bytes()).hexdigest(),
 'force_source_sha256':hashlib.sha256((OUT/'sources.cpp').read_bytes()).hexdigest(),
 'source_snapshot':str(BASE),'tissue_substeps':2,'pressure_pa':600,'tension':1,
 'window_s':[.2,1.2],'aerodynamic_tissue_work_j':float(airwork),'damping_work_j':float(dampwork),
 'upper_surface_pressure_max_pa':float(np.max(abs(force['upper_pressure']))),
 'full_enhanced_same_signal':str(ROOT/'out/source-stage-listening/priority-free-enhanced-cells.png'),
 'enhanced_inspection':'pending','limitations':['No empirical human glottal waveform supplied','Symmetric fixed controls; no noise model','Force work is not a full fluid-tissue energy audit']}
(OUT/'receipt.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
