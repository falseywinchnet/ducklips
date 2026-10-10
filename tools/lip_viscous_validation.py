#!/usr/bin/env python3
"""Collect fresh lip-path checks and verify every nested artifact digest."""
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
    reports=['lip_viscous_reference','lip_viscous_exponential_reference','planar_lip_viscous_reference',
        'planar_lip_viscous_16_reference','planar_lip_viscous_asan_reference','lip_viscous_refinement_location']
    for name in reports:
        path=ROOT/f'reports/{name}.json';verify(json.loads(path.read_text()));verified[str(path.relative_to(ROOT))]=sha(path)
    regressions=ROOT/'out/lip-viscous-regressions/receipt.json';reg=json.loads(regressions.read_text());verify(reg)
    assert len(reg['cases'])==23 and all(case['exact'] for case in reg['cases'])
    tests=[]
    for folder,log in [('build','out/lip-viscous-ctest.log'),('.build/asan','out/lip-viscous-asan-ctest.log')]:
        assert '100% tests passed, 0 tests failed out of 18' in (ROOT/log).read_text()
        tests.append(dict(build=folder,passed=18,failed=0));verified[log]=sha(ROOT/log)
    sanitizer=json.loads((ROOT/'reports/planar_lip_viscous_asan_reference.json').read_text())
    assert len(sanitizer['runs'])==4 and all('.build/asan/' in r['command'][0] for r in sanitizer['runs'])
    orders=[]
    for case in ['closing','opening']:
        for n in [512,1024]:
            pair=[ROOT/f'out/{folder}/{case}-1-{n}.json' for folder in ['planar-lip-viscous','planar-lip-viscous-16']]
            a,b=[json.loads(path.read_text())['steps'][-1] for path in pair]
            assert a['volume']==b['volume'] and a['conductance']==b['conductance']
            delta=[y-x for x,y in zip(a['pressure'],b['pressure'])]
            orders.append(dict(case=case,subdivisions=n,pressure_difference_16_minus_8_Pa=delta,
                maximum_pressure_difference_Pa=max(map(abs,delta))))
    paths=[Path(__file__),regressions,ROOT/'README.md',ROOT/'reports/CONTINUOUS_LIP_VISCOUS.md',ROOT/'CMakeLists.txt',
        ROOT/'out/nato-spatial-pressure.wav']
    for path in paths:verified[str(path.relative_to(ROOT))]=sha(path)
    assert verified['out/nato-spatial-pressure.wav']=='8284f366f2df5f3496d3139f7cbbaf4a8947ad839737d6a7e40b78635083b307'
    result=dict(scope='Actual open internal lip-path drag integration and independent spatial/matrix audits. No new speech or performance claim.',
        ctest=tests,moving_sanitizer_intervals=4,historical_exact_json_regressions=23,
        spatial_order_comparison=orders,verified_artifact_count=len(verified),hashes=dict(sorted(verified.items())))
    (ROOT/'reports/continuous_lip_viscous_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(verified_artifact_count=len(verified),spatial_order_comparison=orders),indent=2))
if __name__=='__main__':main()
