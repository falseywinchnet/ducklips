#!/usr/bin/env python3
"""Collect fresh outlet closure checks and verify every nested artifact digest."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def main():
    verified={}
    def verify(value):
        if isinstance(value,dict):
            for key,item in value.items():
                if key in ['hashes','files','input_poses'] and isinstance(item,dict):
                    for name,digest in item.items():
                        assert sha(ROOT/name)==digest,('stale artifact',name)
                        verified[name]=digest
                else:verify(item)
        elif isinstance(value,list):
            for item in value:verify(item)
    reports=['lip_outlet_viscous_reference','blocked_viscous_reference','viscous_forced_flow_reference','planar_outlet_viscous_reference',
        'planar_outlet_viscous_16_reference','planar_outlet_viscous_asan_reference','outlet_viscous_refinement_location']
    for name in reports:
        path=ROOT/f'reports/{name}.json';verify(json.loads(path.read_text()));verified[str(path.relative_to(ROOT))]=sha(path)
    regressions=ROOT/'out/outlet-viscous-regressions/receipt.json';reg=json.loads(regressions.read_text());verify(reg)
    assert len(reg['cases'])==27 and all(case['exact'] for case in reg['cases'])
    tests=[]
    for folder,log in [('build','out/outlet-viscous-ctest.log'),('.build/asan','out/outlet-viscous-asan-ctest.log')]:
        assert '100% tests passed, 0 tests failed out of 18' in (ROOT/log).read_text()
        tests.append(dict(build=folder,passed=18,failed=0));verified[log]=sha(ROOT/log)
    sanitizer=json.loads((ROOT/'reports/planar_outlet_viscous_asan_reference.json').read_text())
    assert len(sanitizer['runs'])==4 and all('.build/asan/' in r['command'][0] for r in sanitizer['runs'])
    orders=[]
    for case in ['closing','opening']:
        for n in [512,1024]:
            pair=[ROOT/f'out/{folder}/{case}-1-{n}.json' for folder in ['planar-outlet-viscous','planar-outlet-viscous-16']]
            a,b=[json.loads(path.read_text())['steps'][-1] for path in pair]
            assert a['volume']==b['volume'] and a['conductance']==b['conductance']
            delta=[y-x for x,y in zip(a['pressure'],b['pressure'])]
            orders.append(dict(case=case,subdivisions=n,pressure_difference_16_minus_8_Pa=delta,
                maximum_pressure_difference_Pa=max(map(abs,delta))))
    closure=[]
    for case in ['closing','opening']:
        closed_steps=retired=0;surviving_energy=0.;stiffness=[]
        for interval in [0,1]:
            data=json.loads((ROOT/f'out/planar-outlet-viscous/{case}-{interval}-1024.json').read_text())
            previous=data['initial']
            if case=='opening' and interval==0:
                core=[site for site in previous['viscous_sites'] if site['branch']==6][0]
                assert core['segment_inertance']==0 and all(z==0 for z in core['momentum'])
            for step in data['steps']:
                sites=[site for site in step['viscous_sites'] if site['branch']==6]
                prior=[site for site in previous['viscous_sites'] if site['branch']==6]
                assert len(sites)==9
                if sites[0]['segment_inertance']==0:
                    closed_steps+=1;assert sites[0]['half_gap']==0 and sites[0]['branch_conductance']==0 and all(z==0 for z in sites[0]['momentum'])
                    assert step['momentum'][6]==0
                    if prior[0]['segment_inertance']>0:retired+=1
                    surviving_energy=max(surviving_energy,sum(site['energy'] for site in sites[1:]))
                else:
                    resistance=sum(3*site['segment_inertance']*(1.86e-5/1.14)/site['half_gap']**2 for site in sites)
                    stiffness.append(step['dt']*step['conductance'][6]*resistance)
                previous=step
        assert (case=='closing' and retired==1 and closed_steps>0 and surviving_energy>0) or (case=='opening' and retired==0)
        closure.append(dict(case=case,closed_endpoint_steps=closed_steps,core_retirements=retired,maximum_surviving_blocked_profile_energy_J=surviving_energy,
            open_endpoint_steps=len(stiffness),steps_with_endpoint_DC_decay_times_step_above_one=sum(value>1 for value in stiffness),maximum_endpoint_DC_decay_times_step=max(stiffness)))
    paths=[Path(__file__),regressions,ROOT/'README.md',ROOT/'RESEARCH.md',ROOT/'reports/CONTINUOUS_OUTLET_VISCOSITY.md',ROOT/'CMakeLists.txt',
        ROOT/'out/nato-spatial-pressure.wav',ROOT/'include/ducklips/graph_continuous_glottis.hpp',ROOT/'tools/graph_continuous_glottis_probe.cpp']
    for path in paths:verified[str(path.relative_to(ROOT))]=sha(path)
    assert verified['out/nato-spatial-pressure.wav']=='8284f366f2df5f3496d3139f7cbbaf4a8947ad839737d6a7e40b78635083b307'
    result=dict(scope='Outlet core closure and surviving shoulder profile integration, with a verified stiff forced-flow splitting limitation. No new speech or performance claim.',
        ctest=tests,moving_sanitizer_intervals=4,historical_exact_json_regressions=27,closure_state_checks=closure,
        spatial_order_comparison=orders,verified_artifact_count=len(verified),hashes=dict(sorted(verified.items())))
    (ROOT/'reports/continuous_outlet_viscous_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(verified_artifact_count=len(verified),spatial_order_comparison=orders,closure_state_checks=closure),indent=2))
if __name__=='__main__':main()
