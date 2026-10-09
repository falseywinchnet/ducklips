#!/usr/bin/env python3
"""Couple persistent-clear region correspondence to native stored-air state.

A moving sealed wall is a controlled two-chamber test, not speech anatomy.
Its geometric regions deliberately change output order during motion.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import shapely
from shapely.geometry import box
from planar_region_motion_reference import fixture,fixture_regions,correspondence

ROOT=Path(__file__).resolve().parents[1]
BULK=1.14*350**2


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def pose(shift,time):
    p=fixture(box(-.02,-.01,.02,.01),[box(shift-.00005,-.02,shift+.00005,.02)],time)
    p.update(depth=.022,x_edges=[-.02,.02],y_edges=[-.01,.01])
    return p


def execute(divisions,binary,out):
    shifts=[0.,.0002,-.0001,0.]
    poses=[pose(0,0)]
    for step,(a,b) in enumerate(zip(shifts,shifts[1:])):
        for sub in range(1,divisions+1):
            fraction=sub/divisions;poses.append(pose((1-fraction)*a+fraction*b,(step+fraction)*.001))
    all_geometry=[];all_cells=[];volumes=[];sides=[]
    for index,p in enumerate(poses):
        geometry,cells=fixture_regions(p);assert len(geometry)==2
        if index%2:geometry=geometry[::-1]
        all_geometry.append(geometry);all_cells.append(cells)
        volumes.append(np.asarray([g.area*p['depth'] for g in geometry]))
        sides.append([0 if g.centroid.x<0 else 1 for g in geometry])
    initial_pressure=np.asarray([600.,0.]);initial_mass=volumes[0]*(1+initial_pressure/BULK)
    lines=[f'2 {len(poses)-1}']
    lines.extend(f'{v:.17g} {v/BULK:.17g} {pressure:.17g}' for v,pressure in zip(volumes[0],initial_pressure))
    correspondence_records=[];maximum_rejected_overlap=0
    for step in range(1,len(poses)):
        old,new=poses[step-1],poses[step]
        result=correspondence(old,new,all_geometry[step-1],all_cells[step-1],all_geometry[step],all_cells[step])
        assert len(result['one_to_one'])==2 and not result['split_candidates'] and not result['merge_candidates']
        origin={j:i for i,j in result['one_to_one']}
        for i,j in result['one_to_one']:assert sides[step-1][i]==sides[step][j]
        maximum_rejected_overlap=max(maximum_rejected_overlap,max((r[2] for r in result['endpoint_overlaps_removed']),default=0))
        lines.append(f'2 0 {new["time"]-old["time"]:.17g}')
        lines.extend(f'{origin[j]} {v:.17g} {v/BULK:.17g} 0' for j,v in enumerate(volumes[step]))
        correspondence_records.append(result)
    source=out/f'sealed-wall-{divisions}.txt';destination=source.with_suffix('.csv');log=source.with_suffix('.log')
    source.write_text('\n'.join(lines)+'\n')
    with source.open() as f,destination.open('w') as g:
        run=subprocess.run([str(binary)],stdin=f,stdout=g,stderr=subprocess.PIPE,text=True)
    log.write_text(run.stderr);assert run.returncode==0,run.stderr
    data=np.loadtxt(destination,delimiter=',',skiprows=1);pressure=data[:,2].reshape(len(poses),2)
    assert np.array_equal(data[:,0],np.repeat(np.arange(len(poses)),2))
    pressure_error=mass_error=0;density=[]
    for step,v in enumerate(volumes):
        expected_mass=initial_mass[sides[step]]
        expected_pressure=(expected_mass-v)/(v/BULK)
        pressure_error=max(pressure_error,float(max(abs(expected_pressure-pressure[step]))))
        mass_error=max(mass_error,float(max(abs(v*(1+pressure[step]/BULK)-expected_mass))))
        density.extend((1+pressure[step]/BULK).tolist())
    assert pressure_error<1e-7 and mass_error<1e-18
    restored=pressure[-1][np.argsort(sides[-1])]
    assert max(abs(restored-initial_pressure))<1e-7
    correspondence_path=out/f'sealed-wall-{divisions}-correspondence.json'
    correspondence_path.write_text(json.dumps(correspondence_records,indent=2)+'\n')
    result=dict(subdivisions=divisions,steps=len(poses)-1,maximum_pressure_error_Pa=pressure_error,
                maximum_individual_chamber_mass_error_m3=mass_error,
                final_pressure_by_chamber_Pa=restored.tolist(),relative_density_range=[min(density),max(density)],
                maximum_endpoint_overlap_excluded_m3=maximum_rejected_overlap,
                command=[str(binary)],input=str(source.relative_to(ROOT)),input_sha256=sha(source),
                output=str(destination.relative_to(ROOT)),output_sha256=sha(destination),log=str(log.relative_to(ROOT)),log_sha256=sha(log),
                correspondence=str(correspondence_path.relative_to(ROOT)),correspondence_sha256=sha(correspondence_path))
    print(json.dumps(result),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_cavity_region_motion_probe')
    parser.add_argument('--output',type=Path,default=ROOT/'out/cavity-region-motion')
    parser.add_argument('--report',type=Path,default=ROOT/'reports/cavity_region_motion_reference.json')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    results=[execute(n,args.binary,args.output) for n in [1,2,4,8]]
    sources=['include/ducklips/cavity_network.hpp','src/cavity_network.cpp','tests/cavity_network_tests.cpp',
             'tools/cavity_region_motion_probe.cpp','tools/cavity_region_motion_reference.py','tools/planar_region_motion_reference.py']
    report=dict(use='native variable-graph stored-air state with independently established geometric continuation; sealed-wall control, not speech',
                bulk_modulus_Pa=BULK,boundary='two sealed chambers separated by a moving wall; no interchamber passage',
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],runs=results,
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),args.binary]},
                versions=dict(numpy=np.__version__,shapely=shapely.__version__,geos=shapely.geos_version_string))
    args.report.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
