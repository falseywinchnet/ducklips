#!/usr/bin/env python3
"""Independent component/port audit of the native planar air-space geometry.

Checks every pose, region, port and whole-fluid component against GEOS.
This is spatial geometry, not temporal transport or speech rendering.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import shapely
import cavity_graph_reference as pressure_reference
from shapely.geometry import Polygon, box, Point, LineString
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'out/planar-regions'
BINARY = ROOT / 'build/ducklips_planar_regions_probe'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def components(regions, links):
    parent = list(range(len(regions)))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j, _, _, area in links:
        assert area > 0
        parent[root(i)] = root(j)
    grouped = {}
    for i, row in enumerate(regions):
        grouped.setdefault(root(i), []).append(i)
    return list(grouped.values())


def audit(path):
    p = json.loads(path.read_text())
    nx, ny = len(p['x_edges']) - 1, len(p['y_edges']) - 1
    lines = [f'{nx} {ny} {len(p["solids"])} {p["depth"]}', ' '.join(map(str, p['x_edges'])), ' '.join(map(str, p['y_edges']))]
    polygons = [p['domain'], *p['solids']]
    for polygon in polygons:
        lines.append(str(len(polygon)))
        lines.extend(' '.join(map(str, point)) for point in polygon)
    source = OUT / f'pose-{p["time"]:.3f}.txt'
    source.write_text('\n'.join(lines) + '\n')
    destination = source.with_suffix('.json')
    with source.open() as input_file, destination.open('w') as output_file:
        subprocess.run([str(BINARY)], stdin=input_file, stdout=output_file, check=True)
    native = json.loads(destination.read_text())
    result = check_geometry(p,native)
    result.update(geometry=str(path.relative_to(ROOT)), geometry_sha256=sha(path),
                  command=[str(BINARY)], input=str(source.relative_to(ROOT)), input_sha256=sha(source),
                  output=str(destination.relative_to(ROOT)), output_sha256=sha(destination))
    print(json.dumps({k:result[k] for k in ['time','regions','links','split_cells','components']}), flush=True)
    return result


def check_geometry(p,native):
    nx, ny = len(p['x_edges']) - 1, len(p['y_edges']) - 1
    polygons = [Polygon(points) for points in [p['domain'], *p['solids']]]
    assert all(g.is_valid for g in polygons)
    fluid = polygons[0].difference(unary_union(polygons[1:]))
    exact_components = [g for g in shapely.get_parts(fluid) if g.area > 0]
    squares = [box(p['x_edges'][col], p['y_edges'][row], p['x_edges'][col+1], p['y_edges'][row+1])
               for row in range(ny) for col in range(nx)]
    clipped = shapely.intersection(squares, fluid)
    cells = {}
    for i, (cell, volume, x, y) in enumerate(native['regions']):
        assert volume > 0
        cells.setdefault(cell, []).append(i)
    exact_region = {}
    volume_error = 0.
    split = []
    for cell, geometry in enumerate(clipped):
        parts = [g for g in shapely.get_parts(geometry) if g.area > 0]
        ids = cells.get(cell, [])
        assert len(parts) == len(ids), (p['time'], cell, 'component count', len(parts), len(ids))
        if len(ids) > 1:
            split.append(cell)
        used = set()
        for i in ids:
            _, volume, x, y = native['regions'][i]
            matches = [k for k, g in enumerate(parts) if g.covers(Point(x, y))]
            assert len(matches) == 1 and matches[0] not in used, (p['time'], cell, 'representative', i, matches)
            used.add(matches[0])
            exact_region[i] = parts[matches[0]]
            volume_error = max(volume_error, abs(volume - exact_region[i].area * p['depth']))
    expected = {}
    for cell, ids in cells.items():
        row, col = divmod(cell, nx)
        faces = []
        if col + 1 < nx:
            faces.append((cell+1, row*(nx+1)+col+1, 1, LineString([(p['x_edges'][col+1], p['y_edges'][row]), (p['x_edges'][col+1], p['y_edges'][row+1])])))
        if row + 1 < ny:
            faces.append((cell+nx, (row+1)*nx+col, 0, LineString([(p['x_edges'][col], p['y_edges'][row+1]), (p['x_edges'][col+1], p['y_edges'][row+1])])))
        for other, face, vertical, line in faces:
            for i in ids:
                for j in cells.get(other, []):
                    area = exact_region[i].intersection(exact_region[j]).intersection(line).length * p['depth']
                    if area > 0:
                        expected[(min(i,j), max(i,j), face, vertical)] = area
    actual = {}
    for i, j, face, vertical, area in native['links']:
        key = (min(i,j), max(i,j), face, vertical)
        assert key not in actual
        actual[key] = area
    assert actual.keys() == expected.keys(), (p['time'], 'port mismatch', set(actual)-set(expected), set(expected)-set(actual))
    area_error = max((abs(actual[k] - expected[k]) for k in expected), default=0.)
    grouped = components(native['regions'], native['links'])
    exact_volumes = sorted(g.area*p['depth'] for g in exact_components)
    volumes = sorted(sum(native['regions'][i][1] for i in ids) for ids in grouped)
    assert len(volumes) == len(exact_volumes), (p['time'], 'whole-domain components')
    component_error = max(abs(a-b) for a,b in zip(volumes, exact_volumes))
    # Matching counts/volumes alone could hide an incorrect exchange of parts.
    # Every native component must map to exactly one Boolean fluid component.
    mapped = []
    for ids in grouped:
        owners = set()
        for i in ids:
            point = exact_region[i].representative_point()
            matches = [k for k,g in enumerate(exact_components) if g.covers(point)]
            assert len(matches) == 1
            owners.add(matches[0])
        assert len(owners) == 1
        mapped.extend(owners)
    assert len(set(mapped)) == len(exact_components)
    assert volume_error < 1e-18 and area_error < 1e-15 and component_error < 1e-16
    return dict(time=p['time'], regions=len(native['regions']), links=len(native['links']), split_cells=split,
                components=len(grouped), component_volumes_m3=volumes, maximum_region_volume_error_m3=volume_error,
                maximum_port_area_error_m2=area_error, maximum_component_volume_error_m3=component_error)


def pressure_control():
    # Only geometry differs: both models use the same declared linear mobility.
    # Closed outer boundary; each actual chamber starts at uniform pressure.
    p = json.loads((ROOT/'out/planar-airspace/motion/pose-5.032.json').read_text())
    native = json.loads((OUT/'pose-5.032.json').read_text())
    fluid = Polygon(p['domain']).difference(unary_union([Polygon(s) for s in p['solids']]))
    chambers = [g for g in shapely.get_parts(fluid) if g.area > 0]
    assert len(chambers) == 2
    rear = min(range(2), key=lambda i: chambers[i].centroid.x)
    owner = np.array([next(k for k,g in enumerate(chambers) if g.covers(Point(x,y))) for _,_,x,y in native['regions']])
    volume = np.array([r[1] for r in native['regions']])
    bulk, mobility, dt = 1.14*350**2, 1e-7, .001
    initial = np.where(owner == rear, 600., 0.)
    links = []
    for i,j,face,vertical,area in native['links']:
        if vertical:
            column = face % len(p['x_edges'])
            distance = .5*(p['x_edges'][column+1]-p['x_edges'][column-1])
        else:
            row = face // (len(p['x_edges'])-1)
            distance = .5*(p['y_edges'][row+1]-p['y_edges'][row-1])
        links.append((i,j,area,distance/(mobility*area),0.,0.,0.))
    pressure_reference.OUT = OUT
    result, flows, receipt = pressure_reference.run('pressure-separated',volume,volume/bulk,initial,np.zeros(len(volume)),links,dt)
    assert max(abs(result-initial)) == 0 and max(abs(flows)) == 0
    # Reproduce the old merged-cell representation, including its initial
    # pressure mixing. Each cell retains exactly the same total stored mass.
    cells = sorted({int(r[0]) for r in native['regions']})
    index = {cell:i for i,cell in enumerate(cells)}
    merged = np.array([index[int(r[0])] for r in native['regions']])
    merged_volume = np.bincount(merged,weights=volume)
    merged_initial = np.bincount(merged,weights=volume*initial)/merged_volume
    routes = {}
    for i,j,area,resistance,*_ in links:
        key = tuple(sorted((int(merged[i]),int(merged[j]))))
        prior_area, prior_g = routes.get(key,(0.,0.))
        routes[key] = prior_area+area, prior_g+1/resistance
    merged_links = [(i,j,area,1/g,0.,0.,0.) for (i,j),(area,g) in routes.items()]
    old, oldflows, old_receipt = pressure_reference.run('pressure-merged',merged_volume,merged_volume/bulk,
                    merged_initial,np.zeros(len(cells)),merged_links,dt)
    front = owner != rear
    front_gain = float(np.dot(volume[front]/bulk,old[merged[front]]-initial[front]))
    initial_mix = float(np.dot(volume[front]/bulk,merged_initial[merged[front]]-initial[front]))
    mean_front_pressure = float(np.dot(volume[front],old[merged[front]])/sum(volume[front]))
    assert front_gain > 0 and initial_mix > 0 and mean_front_pressure > 0
    return dict(time=5.032, boundary='sealed', initial_rear_Pa=600., initial_front_Pa=0., dt=dt,
                mobility_m2_per_Pa_s=mobility, separated_maximum_pressure_change_Pa=float(max(abs(result-initial))),
                separated_maximum_flow_m3_per_s=float(max(abs(flows))),
                merged_initial_front_gain_m3=initial_mix, merged_final_front_gain_m3=front_gain,
                merged_mean_front_pressure_Pa=mean_front_pressure, separated=receipt, merged=old_receipt)


def main():
    global BINARY, OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=BINARY)
    parser.add_argument('--output',type=Path,default=OUT)
    parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_regions_reference.json')
    args = parser.parse_args()
    BINARY, OUT = args.binary.resolve(), args.output.resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    results = [audit(p) for p in sorted((ROOT/'out/planar-airspace/motion').glob('pose-*.json'))]
    assert len(results) == 66
    report = dict(use='spatial component geometry; no moving-state identity or acoustic speech integration',
                  versions=dict(shapely=shapely.__version__, geos=shapely.geos_version_string, numpy=np.__version__),
                  audit_command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]], poses=results, pressure_control=pressure_control())
    sources = ['include/ducklips/planar_airspace.hpp','src/planar_airspace.cpp','tests/planar_airspace_tests.cpp',
               'tools/planar_regions_probe.cpp','tools/planar_regions_reference.py',str(BINARY.relative_to(ROOT)),
               'tools/cavity_graph_reference.py','build/ducklips_cavity_graph_probe']
    report['hashes'] = {name:sha(ROOT/name) for name in sources}
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print('All 66 spatial decompositions and all their geometric openings agree with the independent geometry.',flush=True)


if __name__ == '__main__':
    main()
