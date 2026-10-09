#!/usr/bin/env python3
"""Independent fixed-geometry pressure audit; does not synthesize speech.

The Foxtrot graphs retain every positive native cell and aperture. Their
circular-passage resistances are a numerical pressure fixture, not a model of
viscous flow in a planar mouth: artificial grid faces are not duct walls.
"""
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import scipy
from scipy.optimize import root as find_root
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'out/cavity-graph'
BINARY = ROOT / 'build/ducklips_cavity_graph_probe'
RHO, SOUND, MU = 1.14, 350., 1.86e-5


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(name, volume, compliance, pressure, inflow, links, dt):
    lines = [f'{len(volume)} {len(links)} 1 {dt:.17g}']
    lines += [' '.join(format(float(x), '.17g') for x in row)
              for row in zip(volume, compliance, pressure, inflow)]
    lines += [f'{i} {j} ' + ' '.join(format(float(x), '.17g') for x in parameters)
              for i, j, *parameters in links]
    source = OUT / f'{name}.txt'
    source.write_text('\n'.join(lines) + '\n')
    result = subprocess.run([str(BINARY)], input=source.read_text(), text=True, capture_output=True, check=True)
    destination = OUT / f'{name}.csv'
    destination.write_text(result.stdout)
    (OUT / f'{name}.log').write_text(result.stderr)
    values = list(csv.DictReader(io.StringIO(result.stdout)))
    p = np.array([float(row['value']) for row in values if row['kind'] == 'pressure'])
    q = np.array([float(row['value']) for row in values if row['kind'] == 'flow'])
    assert len(p) == len(volume) and len(q) == len(links)
    return p, q, {'command': [str(BINARY)], 'input': str(source.relative_to(ROOT)), 'input_sha256': sha(source),
                  'output': str(destination.relative_to(ROOT)), 'output_sha256': sha(destination),
                  'native_report': result.stderr.strip()}


def flow_and_slope(drop, area, resistance, forward, reverse):
    if area == 0:
        return 0., 0.
    loss = forward if drop >= 0 else reverse
    k = loss * RHO / (2 * area * area)
    q = 2 * drop / (resistance + np.sqrt(resistance ** 2 + 4 * k * abs(drop)))
    return q, 1 / (resistance + 2 * k * abs(q))


def balance(p, initial, c, inflow, links, dt):
    residual = c / dt * (p - initial) - inflow
    jacobian = np.diag(c / dt)
    for i, j, area, resistance, reservoir, forward, reverse in links:
        q, g = flow_and_slope(p[i] - (p[j] if j >= 0 else reservoir), area, resistance, forward, reverse)
        residual[i] += q
        jacobian[i, i] += g
        if j >= 0:
            residual[j] -= q
            jacobian[j, j] += g
            jacobian[i, j] -= g
            jacobian[j, i] -= g
    return residual, jacobian


def mass_error(p, initial, c, inflow, links, q, dt):
    external = sum(flow for link, flow in zip(links, q) if link[1] == -1)
    return abs(c @ (p - initial) + dt * (external - sum(inflow)))


def nonlinear_audit():
    random = np.random.default_rng(7319)
    results = []
    for case in range(24):
        n, dt = 12, .0002
        volume = np.full(n, 1e-5)
        c = 10 ** random.uniform(-10, -8, n)
        initial = random.uniform(-500, 500, n)
        inflow = random.uniform(-1e-6, 1e-6, n)
        pairs = [(i, (i + 1) % n) for i in range(n)] + [(i, i + 4) for i in range(8)] + [(0, -1), (n - 1, -1)]
        links = []
        for e, (i, j) in enumerate(pairs):
            area = 10 ** random.uniform(-6, -4)
            if e < 20 and (case + e) % 7 == 0:
                area = 0.
            links.append((i, j, area, 10 ** random.uniform(3, 6), (-1 if case % 2 else 1) * 600 if j == -1 else 0.,
                          float(e % 3), float((e + 1) % 3)))
        p, q, receipt = run(f'nonlinear-{case:02}', volume, c, initial, inflow, links, dt)
        oracle = find_root(lambda x: balance(x, initial, c, inflow, links, dt)[0], initial,
                           jac=lambda x: balance(x, initial, c, inflow, links, dt)[1], tol=1e-11)
        oracle_residual = max(abs(balance(oracle.x, initial, c, inflow, links, dt)[0]))
        error = max(abs(p - oracle.x))
        mass = mass_error(p, initial, c, inflow, links, q, dt)
        assert oracle_residual < 1e-12 and error < 1e-5 and mass < 1e-15, (case, error, mass, oracle_residual)
        results.append(dict(case=case, pressure_error_Pa=error, mass_error_m3=mass,
                            oracle_residual_m3_per_s=oracle_residual, **receipt))
    return results


def snapshot_audit(time):
    source = ROOT / f'out/planar-airspace/pose-{time:.3f}.json'
    pose = json.loads(source.read_text())
    nx, ny = len(pose['x_edges']) - 1, len(pose['y_edges']) - 1
    cells = np.flatnonzero(np.array(pose['volume']) > 0)
    index = {int(cell): i for i, cell in enumerate(cells)}
    volume = np.array(pose['volume'])[cells]
    c = volume / (RHO * SOUND ** 2)
    links = []
    for cell, i in index.items():
        row, column = divmod(cell, nx)
        candidates = []
        if column + 1 < nx:
            candidates.append((cell + 1, pose['vertical_face'][row * (nx + 1) + column + 1],
                               .5 * (pose['x_edges'][column + 2] - pose['x_edges'][column])))
        if row + 1 < ny:
            candidates.append((cell + nx, pose['horizontal_face'][(row + 1) * nx + column],
                               .5 * (pose['y_edges'][row + 2] - pose['y_edges'][row])))
        for other, area, length in candidates:
            if area > 0:
                assert other in index
                links.append((i, index[other], area, 8 * np.pi * MU * length / area ** 2, 0., 0., 0.))
    rows, columns, values = [], [], []
    for i, j, area, resistance, *_ in links:
        rows.extend([i, j, i, j]); columns.extend([i, j, j, i]); values.extend([1 / resistance, 1 / resistance, -1 / resistance, -1 / resistance])
    laplacian = coo_matrix((values, (rows, columns)), shape=(len(c), len(c))).tocsc()
    count, labels = connected_components(laplacian, directed=False)
    initial = 400 + 200 * np.sin(cells / nx)
    rear = None
    if time == 5.040:
        assert count == 2
        x = np.array([.5 * (pose['x_edges'][cell % nx] + pose['x_edges'][cell % nx + 1]) for cell in cells])
        rear = int(labels[np.argmin(x)])
        initial[labels == rear] = 600.
        initial[labels != rear] -= 300.
    else:
        assert count == 1
    dt = .0001
    inflow = np.zeros(len(c))
    p, q, receipt = run(f'foxtrot-{time:.3f}', volume, c, initial, inflow, links, dt)
    oracle = spsolve(diags(c / dt) + laplacian, c / dt * initial)
    error = max(abs(p - oracle))
    component_mass = [abs(c[labels == k] @ (p - initial)[labels == k]) for k in range(count)]
    # Backward Euler energy identity for a closed, fixed, linear network.
    change = .5 * np.dot(c, (p - initial) * (p + initial))
    temporal_loss = .5 * np.dot(c, (p - initial) ** 2)
    viscous_loss = dt * sum(link[3] * flow ** 2 for link, flow in zip(links, q))
    energy_error = abs(change + temporal_loss + viscous_loss)
    isolated_error = max(abs(p[labels == rear] - 600)) if rear is not None else None
    assert error < 1e-5 and max(component_mass) < 1e-15 and energy_error < 1e-12
    assert isolated_error is None or isolated_error < 1e-8
    return dict(time=time, geometry=str(source.relative_to(ROOT)), geometry_sha256=sha(source),
                nodes=len(c), edges=len(links), components=count, cycle_rank=len(links)-len(c)+count,
                pressure_error_Pa=error, component_mass_error_m3=component_mass,
                energy_balance_error_J=energy_error, isolated_rear_pressure_error_Pa=isolated_error, **receipt)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = dict(use='fixed-geometry mean-pressure graph reference; not speech, moving-state transport, or calibrated planar viscosity',
                  audit_command=[sys.executable, str(Path(__file__).resolve())],
                  scipy_version=scipy.__version__, numpy_version=np.__version__,
                  nonlinear=nonlinear_audit(), foxtrot=[snapshot_audit(t) for t in (5.025, 5.040, 5.071, 5.073)])
    paths = ['include/ducklips/cavity_network.hpp', 'src/cavity_network.cpp', 'tools/cavity_graph_probe.cpp',
             'tests/cavity_network_tests.cpp', 'tools/cavity_graph_reference.py', 'build/ducklips_cavity_graph_probe']
    report['hashes'] = {name: sha(ROOT / name) for name in paths}
    (ROOT / 'reports/cavity_graph_reference.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'nonlinear_max_pressure_error_Pa': max(r['pressure_error_Pa'] for r in report['nonlinear']),
                      'nonlinear_max_mass_error_m3': max(r['mass_error_m3'] for r in report['nonlinear']),
                      'foxtrot': report['foxtrot']}, indent=2))


if __name__ == '__main__':
    main()
