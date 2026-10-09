#!/usr/bin/env python3
"""Conservative zero-storage pressure tests on audited moving Foxtrot poses.

Uses an explicit linear mobility fixture, not calibrated vocal-tract friction.
Rejects every pose containing a disconnected fluid fragment inside a grid cell.
No speech or acoustic state is generated. Temporal refinement interpolates the
saved volumes/apertures; it does not reconstruct true intermediate anatomy.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import shapely
import scipy
from shapely.geometry import Polygon, box
from shapely.ops import unary_union
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'out/cavity-transport'
BINARY = ROOT / 'build/ducklips_cavity_transport_probe'
BULK = 1.14 * 350 ** 2
MOBILITY = 1e-7  # m^2/(Pa s), a prescribed numerical fixture, not anatomy


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def topology(p):
    nx, ny = len(p['x_edges']) - 1, len(p['y_edges']) - 1
    links, length, aperture = [], [], []
    for row in range(ny):
        for column in range(nx):
            cell = row * nx + column
            if column + 1 < nx:
                links.append((cell, cell + 1))
                length.append(.5 * (p['x_edges'][column + 2] - p['x_edges'][column]))
                aperture.append(p['vertical_face'][row * (nx + 1) + column + 1])
            if row + 1 < ny:
                links.append((cell, cell + nx))
                length.append(.5 * (p['y_edges'][row + 2] - p['y_edges'][row]))
                aperture.append(p['horizontal_face'][(row + 1) * nx + column])
    return np.array(links), np.array(length), np.array(aperture)


def audit_geometry(path):
    p = json.loads(path.read_text())
    polygons = [Polygon(p['domain']), *[Polygon(s) for s in p['solids']]]
    assert all(g.is_valid for g in polygons), f'invalid polygon in {path}'
    fluid = polygons[0].difference(unary_union(polygons[1:]))
    regions = [g for g in shapely.get_parts(fluid) if g.area > 0]
    nx, ny = len(p['x_edges']) - 1, len(p['y_edges']) - 1
    squares = [box(p['x_edges'][column], p['y_edges'][row], p['x_edges'][column + 1], p['y_edges'][row + 1])
               for row in range(ny) for column in range(nx)]
    clipped = shapely.intersection(squares, fluid)
    exact_volume = shapely.area(clipped) * p['depth']
    error = float(max(abs(exact_volume - p['volume'])))
    assert error < 1e-18
    split = []
    for cell, g in enumerate(clipped):
        pieces = [part for part in shapely.get_parts(g) if part.area > 0]
        if len(pieces) > 1:
            split.append(dict(cell=cell, volume_m3=[part.area * p['depth'] for part in pieces],
                              fluid_regions=[next(i for i, region in enumerate(regions) if region.covers(part.representative_point())) for part in pieces]))
    links, length, aperture = topology(p)
    active = np.array(p['volume']) > 0
    graph = coo_matrix((aperture, (links[:, 0], links[:, 1])), shape=(len(active), len(active))).tocsr()
    graph.eliminate_zeros()
    count = connected_components(graph[active][:, active], directed=False, return_labels=False)
    if not split:
        assert count == len(regions)
    return p, dict(time=p['time'], sha256=sha(path), split_cells=split, exact_fluid_components=len(regions),
                   merged_cell_graph_components=int(count), maximum_cell_volume_error_m3=error)


def execute(name, poses, subdivisions):
    links, length, initial_aperture = topology(poses[0])
    n = len(poses[0]['volume'])
    volumes = [np.array(poses[0]['volume'])]
    apertures = [initial_aperture]
    times = [poses[0]['time']]
    for old, new in zip(poses, poses[1:]):
        oldv, newv = np.array(old['volume']), np.array(new['volume'])
        olda, newa = topology(old)[2], topology(new)[2]
        for substep in range(1, subdivisions + 1):
            f = substep / subdivisions
            volumes.append((1 - f) * oldv + f * newv)
            apertures.append((1 - f) * olda + f * newa)
            times.append((1 - f) * old['time'] + f * new['time'])
    source = OUT / f'{name}-{subdivisions}.txt'
    with source.open('w') as output:
        output.write(f'{n} {len(links)} {len(volumes)}\n')
        for i, j in links:
            output.write(f'{i} {j}\n')
        for v in volumes[0]:
            output.write(f'{v:.17g} {v / BULK:.17g} 0\n')
        for step in range(1, len(volumes)):
            output.write(f'{times[step] - times[step - 1]:.17g}\n')
            for v in volumes[step]:
                output.write(f'{v:.17g} {v / BULK:.17g} 0\n')
            area = .5 * (apertures[step - 1] + apertures[step])
            for a, distance in zip(area, length):
                resistance = distance / (MOBILITY * a) if a > 0 else 0.
                output.write(f'{a:.17g} {resistance:.17g} 0 0 0\n')
    destination = OUT / f'{name}-{subdivisions}.csv'
    with source.open() as input_file, destination.open('w') as output:
        result = subprocess.run([str(BINARY)], stdin=input_file, stdout=output, stderr=subprocess.PIPE, text=True)
    (OUT / f'{name}-{subdivisions}.log').write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(result.stderr)
    data = np.loadtxt(destination, delimiter=',', skiprows=1)
    assert len(data) == len(volumes) * n
    assert np.array_equal(data[:, 0].astype(int), np.repeat(np.arange(len(volumes)), n))
    assert np.array_equal(data[:, 1].astype(int), np.tile(np.arange(n), len(volumes)))
    pressures = data[:, 2].reshape(len(volumes), n)
    mass0 = sum(volumes[0])
    error = mass_error = local_error = 0.
    min_density = np.inf
    max_density = 0.
    born = vanished = 0
    oracle = np.zeros(n)
    for step in range(1, len(volumes)):
        dt = times[step] - times[step - 1]
        v0, v1 = volumes[step - 1], volumes[step]
        c0, c1 = v0 / BULK, v1 / BULK
        area = .5 * (apertures[step - 1] + apertures[step])
        conductance = MOBILITY * area / length
        i, j = links.T
        values = np.concatenate([conductance, conductance, -conductance, -conductance])
        laplacian = coo_matrix((values, (np.concatenate([i, j, i, j]), np.concatenate([i, j, j, i]))), shape=(n, n)).tocsc()
        matrix = diags(c1 / dt) + laplacian
        active = matrix.diagonal() > 0
        right = (c0 * oracle - (v1 - v0)) / dt
        assert np.all(right[~active] == 0)
        oracle = np.zeros(n)
        oracle[active] = spsolve(matrix[active][:, active], right[active])
        assert np.all(np.isfinite(oracle))
        p = pressures[step]
        error = max(error, float(max(abs(oracle - p))))
        charge = v1 + c1 * p
        mass_error = max(mass_error, float(abs(sum(charge) - mass0)))
        old_charge = v0 + c0 * pressures[step - 1]
        residual = charge - old_charge + dt * (laplacian @ p)
        local_error = max(local_error, float(max(abs(residual))))
        min_density = min(min_density, float(min(charge[v1 > 0] / v1[v1 > 0])))
        max_density = max(max_density, float(max(charge[v1 > 0] / v1[v1 > 0])))
        born += np.count_nonzero((v0 == 0) & (v1 > 0))
        vanished += np.count_nonzero((v0 > 0) & (v1 == 0))
    assert error < 1e-4 and mass_error < 1e-13 and local_error < 1e-13 and min_density > 0
    receipt = dict(name=name, start=times[0], end=times[-1], subdivisions=subdivisions, steps=len(volumes)-1,
                   born=int(born), vanished=int(vanished), maximum_pressure_error_Pa=error,
                   maximum_total_mass_error_m3=mass_error, maximum_local_mass_error_m3=local_error,
                   minimum_relative_density=min_density, maximum_relative_density=max_density, native_report=result.stderr.strip(),
                   command=[str(BINARY)], input=str(source.relative_to(ROOT)), input_sha256=sha(source),
                   output=str(destination.relative_to(ROOT)), output_sha256=sha(destination))
    print(json.dumps(receipt), flush=True)
    return receipt, pressures[-1]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    poses, audit = zip(*[audit_geometry(path) for path in sorted((ROOT / 'out/planar-airspace/motion').glob('pose-*.json'))])
    assert len(poses) == 66
    invalid = {r['time'] for r in audit if r['split_cells']}
    flags = [np.array(p['volume']) > 0 for p in poses]
    born = sum(int(np.count_nonzero(b & ~a)) for a, b in zip(flags, flags[1:]))
    vanished = sum(int(np.count_nonzero(a & ~b)) for a, b in zip(flags, flags[1:]))
    runs, refinement = [], []
    for name, start, end in [('preclosure', 5.020, 5.031), ('release', 5.038, 5.078)]:
        selected = [p for p in poses if start <= p['time'] <= end]
        assert len(selected) == round((end - start) * 1000) + 1
        assert not any(p['time'] in invalid for p in selected), 'split cells cannot be pressure nodes'
        previous, differences, weighted, worst = None, [], [], []
        for divisions in (1, 2, 4):
            receipt, final = execute(name, selected, divisions)
            runs.append(receipt)
            if previous is not None:
                # Algebraic disappearing-cell multipliers are not physical
                # endpoint pressures; compare only positive final air volumes.
                mask = np.array(selected[-1]['volume']) > 0
                delta = np.where(mask, final - previous, 0)
                volume = np.array(selected[-1]['volume'])
                cell = int(np.argmax(abs(delta)))
                differences.append(float(max(abs(delta))))
                weighted.append(float(np.sqrt(np.dot(volume, delta * delta) / sum(volume))))
                worst.append(dict(cell=cell, volume_m3=float(volume[cell]), previous_pressure_Pa=float(previous[cell]), pressure_Pa=float(final[cell])))
            previous = final
        refinement.append(dict(name=name, endpoint_pressure_differences_Pa=differences, volume_weighted_RMS_differences_Pa=weighted, worst_cells=worst))
    report = dict(use='mean-pressure mass transport reference; no acoustic integration or new speech',
                  whole_transition_valid=False, pressure_refinement_established=False,
                  decision='Do not promote this trajectory into speech: unresolved split cells, temporal geometry and large density changes',
                  mobility_m2_per_Pa_s=MOBILITY, boundaries='sealed outer boundary',
                  versions=dict(numpy=np.__version__, scipy=scipy.__version__, shapely=shapely.__version__, geos=shapely.geos_version_string),
                  all_pose_cell_events=dict(born=born, vanished=vanished),
                  geometry_command=['./build/ducklips_airspace_geometry_probe', 'out/nato-airspace-geometry.csv', '5.020', '5.085', '512', 'out/planar-airspace/motion.csv', '--frames', 'out/planar-airspace/motion'],
                  temporal_geometry='linear interpolation of saved volumes and apertures; trapezoid interval aperture',
                  audit_command=[sys.executable, str(Path(__file__).resolve())], geometry_audit=audit,
                  rejected_times=sorted(invalid), runs=runs, refinement=refinement)
    sources = ['include/ducklips/cavity_network.hpp', 'src/cavity_network.cpp', 'tests/cavity_network_tests.cpp',
               'tools/airspace_geometry_probe.cpp', 'tools/cavity_transport_probe.cpp', 'tools/cavity_transport_reference.py',
               'build/ducklips_cavity_transport_probe', 'build/ducklips_airspace_geometry_probe', 'out/nato-airspace-geometry.csv',
               'out/cavity-transport-native.log', 'out/cavity-transport-native-asan.log']
    report['hashes'] = {name:sha(ROOT/name) for name in sources}
    (ROOT/'reports/cavity_transport_reference.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'rejected_times': sorted(invalid), 'refinement': refinement}, indent=2), flush=True)


if __name__ == '__main__':
    main()
