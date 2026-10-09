#!/usr/bin/env python3
"""Isolate contact localization from stiff viscous decay in the current source."""
import argparse
import concurrent.futures
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
from continuous_phonation_reference import Reference, YSCALE, VSCALE

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'out/continuous-tract-contact'
RECEIPT = ROOT/'reports/continuous_tract_contact_reference.json'
CONFIGS = [('baseline-6400',6400,False,False),('local-6400',6400,True,False),
    ('local-12800',12800,True,False),('fitted-3200',3200,True,True),
    ('fitted-6400',6400,True,True),('fitted-12800',12800,True,True)]

def timescales(states):
    result = []
    for state in states:
        y = np.array(state['source'])/YSCALE
        pressure = state['tract_pressure_Pa']
        records = []
        for order,h in [(96,1e-6),(192,5e-7)]:
            ref = Reference(order)
            def rhs(z):
                physical = z*YSCALE
                c = ref.configuration(physical[2:4],physical[:2])
                derivative = ref.rhs(z,1400.)
                load = np.array([0.,0.,-pressure,0.,0.])
                derivative[4:] += np.linalg.solve(c['G']*VSCALE[:,None]*VSCALE[None,:],VSCALE*load)
                return derivative
            jacobian = np.column_stack([(rhs(y+h*e)-rhs(y-h*e))/(2*h) for e in np.eye(9)])
            modes = np.linalg.eigvals(jacobian)
            fastest = min(modes,key=lambda value:value.real)
            records.append(dict(quadrature_order=order,difference_step=h,fastest_real_s=float(fastest.real),
                fastest_imag_s=float(fastest.imag),maximum_frequency_Hz=float(max(abs(modes.imag))/(2*np.pi))))
        assert abs(records[0]['fastest_real_s']/records[1]['fastest_real_s']-1) < 1e-5
        result.append(dict(**state,refinements=records))
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--collect',action='store_true')
    args = parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if RECEIPT.exists():
        previous = json.loads(RECEIPT.read_text())
        seed, samples = previous['seed'],previous['timescale_inputs']
    else:
        seed = json.loads((OUT/'seed.json').read_text())
        rows = list(map(json.loads,(ROOT/'out/graph-continuous-glottis/tube-19200.jsonl').read_text().splitlines()))
        samples = []
        for t in [.008,.0083,.0084,.008425]:
            row = min(rows,key=lambda x:abs(x['time']-t))
            samples.append(dict(time=row['time'],source=row['source'],tract_pressure_Pa=row['charge'][0]/(2e-4*.175/8/139650.)))
    text = ' '.join(map(str,seed['source']+seed['charge']+seed['momentum']))+'\n'
    (OUT/'seed.input').write_text(text)
    current = ROOT/'build/ducklips_continuous_tract_probe'
    archived = OUT/'control-probe'
    controls = archived if archived.exists() else current
    commands = {}
    for name,n,localized,fitted in CONFIGS:
        commands[name] = [str(current if fitted else controls),str(n),'.016','--graph-state']
        if localized:
            commands[name] += ['--contact-step','7.8125001e-8']
        if fitted:
            commands[name] += ['--relaxation']
    def run(name):
        result = subprocess.run(commands[name],input=text,text=True,capture_output=True)
        (OUT/(name+'.jsonl')).write_text(result.stdout)
        (OUT/(name+'.stderr')).write_text(result.stderr)
        assert result.returncode == 0,result.stderr
    if not args.collect:
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            list(pool.map(run,commands))
    runs = {}
    for name,n,localized,fitted in CONFIGS:
        summary = json.loads((OUT/(name+'.stderr')).read_text())
        rows = list(map(json.loads,(OUT/(name+'.jsonl')).read_text().splitlines()))
        assert summary['complete'] and rows[-1]['step'] == n
        assert np.max(abs((np.array(rows[0]['source'])-seed['source'])/YSCALE)) < 1e-8
        assert rows[0]['charge'] == seed['charge'] and rows[0]['momentum'] == seed['momentum']
        assert summary['maximum_work_residual_J'] < 1e-12 and summary['maximum_mass_residual_m3'] < 1e-19
        runs[name] = dict(summary=summary,final=rows[-1],command=commands[name],localized_contact=localized,fitted_damping=fitted)
    differences = []
    for a,b in [('baseline-6400','local-12800'),('local-6400','local-12800'),
                ('fitted-3200','fitted-6400'),('fitted-6400','fitted-12800'),('local-12800','fitted-12800')]:
        x,y = runs[a]['final'],runs[b]['final']
        differences.append(dict(cases=[a,b],maximum_scaled_source_difference=float(np.max(abs((np.array(x['source'])-y['source'])/YSCALE))),
            maximum_tract_pressure_difference_Pa=float(np.max(abs(np.array(x['charge'])-y['charge']))/(2e-4*.175/8/139650.))))
    paths = [Path(__file__),ROOT/'tools/continuous_tract_probe.cpp',ROOT/'src/continuous_glottis.cpp',
        ROOT/'include/ducklips/continuous_glottis.hpp',ROOT/'src/glottal_relaxation.cpp',ROOT/'include/ducklips/glottal_relaxation.hpp',
        ROOT/'tools/continuous_phonation_reference.py',current,controls,OUT/'seed.input',
        *[OUT/(name+'.'+extension) for name,_,_,_ in CONFIGS for extension in ['jsonl','stderr']]]
    report = dict(scope='Common-state restart at 8 ms of the frozen-tube control; separates event-step localization and mode-fitted viscous decay. Neither alone establishes converged contact motion.',
        seed=seed,seed_origin='3a9eb8b tube-19200 at 0.008 s',elapsed_duration_s=.016,
        timescale_scope='Local open-source Jacobian with tract effort held at the recorded value; not equilibrium stability or a full coupled-tract eigensystem',
        timescale_inputs=samples,timescales=timescales(samples),native=runs,comparisons=differences,
        temporal_convergence_established=False,
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    RECEIPT.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(comparisons=differences,timescales=report['timescales']),indent=2))

if __name__ == '__main__':
    main()
