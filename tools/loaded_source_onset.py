#!/usr/bin/env python3
"""Open-source equilibrium and local onset with the declared airway impedance.

Uses independent continuous source equations and eliminates passive memories
analytically in the Laplace domain. This is a local mode diagnostic, not a
through-contact simulation or an anatomical validation.
"""
import json
from pathlib import Path
import re

import numpy as np
from scipy.optimize import least_squares, root
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

from continuous_outlet_column_reference import ColumnReference
from continuous_phonation_reference import YSCALE, VSCALE
from planar_source_reference import EDGES
from slit_thermal_reference import spectrum as thermal_spectrum
from slit_viscous_reference import spectrum as viscous_spectrum
from coupled_phonation import sha

ROOT = Path(__file__).resolve().parents[1]


class Load:
    def __init__(self, row, edges, port, spectra):
        self.c = np.array(row['compliance'])
        self.g = np.array(row['conductance'])
        self.port = port
        self.edges = edges
        self.b = np.zeros((len(self.c), len(self.g)))
        for k, (i, j) in enumerate(edges):
            self.b[i, k] = 1
            if j >= 0:self.b[j, k] = -1
        tw, tp, vw, vp = spectra
        self.thermal = [(site['node'], .4 * tw * site['compliance'],
                         2.3e-5 / site['half_gap']**2 * tp) for site in row['thermal_sites']]
        self.viscous = []
        self.dc = np.zeros(len(self.g))
        for site in row['viscous_sites']:
            l, gap, branch = site['segment_inertance'], site['half_gap'], site['branch']
            if l == 0:continue
            rate = (1.86e-5 / 1.14) / gap**2
            self.dc[branch] += 3 * l * rate
            self.viscous.append((branch, l * rate * vw, rate * vp))
        self.rad_branch = next(k for k, (_, j) in enumerate(edges) if j < 0)
        h = row['radiation_state']['inverse_inertance']
        self.rad_r = 1 / ((1.4 * 1.14 / (4 * np.pi * 350.)) * h*h)
        self.rad_rate = self.rad_r*h

    def impedance(self, s):
        compliance = self.c.astype(complex)
        z = s / self.g + self.dc
        for node, storage, rates in self.thermal:
            compliance[node] += np.sum(storage * rates / (s + rates))
        for branch, resistance, rates in self.viscous:
            z[branch] += np.sum(resistance * s / (s + rates))
        z[self.rad_branch] += self.rad_r * s / (s + self.rad_rate)
        admittance = np.diag(s * compliance) + (self.b / z) @ self.b.T
        unit = np.zeros(len(self.c));unit[self.port] = 1.
        return np.linalg.solve(admittance, unit)[self.port]

    def dc_impedance(self):
        # Contract zero-resistance edges instead of assigning a tiny invented
        # resistance or choosing an arbitrary circulation in a lossless loop.
        n = len(self.c)
        parent = list(range(n + 1))
        def find(i):
            while parent[i] != i:i = parent[i]
            return i
        for r, (i, j) in zip(self.dc, self.edges):
            if r == 0:parent[find(i)] = find(n if j < 0 else j)
        ground = find(n)
        groups = sorted({find(i) for i in range(n)} - {ground})
        index = {group: i for i, group in enumerate(groups)}
        matrix = np.zeros((len(groups), len(groups)))
        for r, (i, j) in zip(self.dc, self.edges):
            if r == 0:continue
            a, b = find(i), find(n if j < 0 else j)
            vector = np.zeros(len(groups))
            if a != ground:vector[index[a]] += 1
            if b != ground:vector[index[b]] -= 1
            matrix += np.outer(vector, vector) / r
        unit = np.zeros(len(groups));unit[index[find(self.port)]] = 1.
        return float(np.linalg.solve(matrix, unit)[index[find(self.port)]])

    def full_memory_impedance(self, s):
        # Separate frequency-domain solve with every physical memory explicit.
        # Unknowns are node pressure, branch flow, then thermal/viscous/radiation
        # memories. This checks the rational elimination and port sign.
        n, m = len(self.c), len(self.g)
        count = n + m + sum(len(a) for _, a, _ in self.thermal + self.viscous) + 1
        matrix = lil_matrix((count, count), dtype=complex)
        matrix[:n, :n] = np.diag(s*self.c)
        matrix[:n, n:n+m] = self.b
        matrix[n:n+m, :n] = -self.b.T
        matrix[n:n+m, n:n+m] = np.diag(s/self.g + self.dc)
        index = n + m
        for node, storage, rates in self.thermal:
            for k, rate in zip(storage, rates):
                matrix[node, index] = s
                matrix[index, node] = -k*rate
                matrix[index, index] = s+rate
                index += 1
        for branch, resistance, rates in self.viscous:
            for r, rate in zip(resistance, rates):
                matrix[n+branch, index] = s
                matrix[index, n+branch] = -r
                matrix[index, index] = s+rate
                index += 1
        matrix[n+self.rad_branch, index] = s
        matrix[index, n+self.rad_branch] = -self.rad_r
        matrix[index, index] = s+self.rad_rate
        rhs = np.zeros(count);rhs[self.port] = 1e-4
        solution = spsolve(matrix.tocsc(), rhs)
        residual = float(np.max(abs(matrix @ solution - rhs)))
        return solution[self.port]/1e-4, count, residual


def mode(load, column, order, difference, guess, pressure=1400.):
    ref = ColumnReference(order);ref.inertance = column
    z0 = load.dc_impedance()
    def stationary(x):
        force = ref.stationary(x, pressure)
        force[2] += z0*x[4]*1e-4/1000.
        return force
    sol = least_squares(stationary, guess, bounds=([-100, -100, .001, .001, 0.], [100, 100, 100, 100, 100]),
                        xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=100)
    residual = float(np.max(abs(stationary(sol.x))))
    if residual > 1e-8:raise RuntimeError(('equilibrium residual', residual))
    y = np.r_[sol.x[:4], np.full(3, sol.x[4]), 0., 0.]
    load_pressure = z0*sol.x[4]*1e-4
    def coupling(y):
        state = ref.configuration(y[2:4]*1e-6, y[:2]*1000.)
        force = np.array([0., 0., -1., 0., 0.])
        return np.r_[np.zeros(4), np.linalg.solve(state['G'] * VSCALE[:, None]*VSCALE[None, :], VSCALE*force)]
    def rhs(y):return ref.rhs(y, pressure) + coupling(y)*load_pressure
    a = np.column_stack([(rhs(y + np.eye(9)[j]*difference) - rhs(y - np.eye(9)[j]*difference))/(2*difference)
                         for j in range(9)])
    b = coupling(y)
    c = np.zeros(9);c[6] = YSCALE[6]
    dc_modes = sorted([dict(growth_s=float(s.real), frequency_Hz=float(s.imag/(2*np.pi)))
                       for s in np.linalg.eigvals(a + z0*np.outer(b, c)) if s.imag > 0],
                      key=lambda m:m['growth_s'], reverse=True)
    def secular(value):
        s = complex(*value)
        feedback = load.impedance(s) * (c @ np.linalg.solve(s*np.eye(9) - a, b))
        value = 1 - feedback
        return np.array([value.real, value.imag])
    solutions = []
    for initial in [[-30., 2*np.pi*320], [30., 2*np.pi*200], [-100., 2*np.pi*120]]:
        found = root(secular, initial, tol=1e-10)
        error = float(np.max(abs(secular(found.x))))
        s = complex(*found.x)
        if error < 1e-7 and s.imag > 0 and not any(abs(s-complex(*r['pole'])) < 1e-3 for r in solutions):
            solutions.append(dict(pole=list(found.x), growth_s=float(s.real), frequency_Hz=float(s.imag/(2*np.pi)),
                                  secular_residual=error))
    return dict(pressure_Pa=pressure, quadrature_order=order, derivative_step=difference, dc_load_Pa_s_m3=z0,
                equilibrium_scaled=y.tolist(), equilibrium_residual=residual,
                equilibrium_rhs_residual=float(np.max(abs(rhs(y)))), constant_load_pressure_Pa=load_pressure,
                instantaneous_resistance_modes=dc_modes, modes=solutions)


def main():
    path = ROOT / 'out/actual-coupled-phonation-v2/final-hold-017-100.json'
    row = json.loads(path.read_text())['steps'][-1]
    source = row['continuous_source']
    expected_body = dict(lower_mass=1.25e-4, upper_mass=2.5e-5, lower_stiffness=80.,
                         upper_stiffness=8., coupling_stiffness=25., lower_damping_ratio=.1,
                         upper_damping_ratio=.6, length=.014, lower_thickness=.0025,
                         upper_thickness=.0005, asymmetry=1.)
    body = dict(line.split('=', 1) for line in (ROOT/'bodies/reference.body').read_text().splitlines()
                if '=' in line and not line.lstrip().startswith('#'))
    assert body['preset'] == 'baritone' and source['rest_area'] == 3e-6
    header = (ROOT/'include/ducklips/sources.hpp').read_text()
    for name, expected in expected_body.items():
        default = float(re.search(r'double\s+' + name + r'\{([^}]+)', header).group(1))
        assert float(body.get('folds.'+name, default)) == expected, name
    for key in ['lower_position', 'upper_position', 'lower_velocity', 'upper_velocity']:
        assert abs(source[key][0] - source[key][1]) < 1e-14
    tw, tp = thermal_spectrum();vw, vp = viscous_spectrum();spectra = tw, tp, vw, vp
    actual = Load(row, EDGES, 13, spectra)
    # Same eight-section tube and loss layout as joint_contact_probe.cpp.
    n, area, dx, rho = 8, 2e-4, .175/8, 1.14
    c = np.full(n, area*dx/(rho*350**2));g = np.full(n, area/(rho*dx));g[-1] *= 2
    h = np.sqrt(np.pi*area)/(rho*.7)
    tube = dict(compliance=c, conductance=g,
                thermal_sites=[dict(node=i, compliance=c[i], half_gap=area/.06) for i in range(n)],
                viscous_sites=[dict(branch=i, segment_inertance=1/g[i], half_gap=area/.06) for i in range(n)],
                radiation_state=dict(inverse_inertance=h))
    tube_load = Load(tube, [(i, i+1 if i+1<n else -1) for i in range(n)], 0, spectra)
    guess = np.r_[np.array(source['effort'])/1000, np.array(source['area'])/1e-6, source['flow'][2]/1e-4]
    cases = []
    for name, load, column in [('actual-mouth', actual, source['outlet_inertance']),
                               ('fixed-tube', tube_load, rho*dx/(2*area))]:
        records, impedance_checks = [], []
        for frequency in [20., 150., 320., 1000., 8000.]:
            s = 2j*np.pi*frequency
            rational = load.impedance(s)
            full, count, residual = load.full_memory_impedance(s)
            error = float(abs(full-rational)/abs(full))
            assert error < 1e-9 and rational.real > 0
            impedance_checks.append(dict(frequency_Hz=frequency, real_impedance=float(rational.real),
                                         imaginary_impedance=float(rational.imag), explicit_unknowns=count,
                                         relative_discrepancy=error, full_equation_residual=residual))
        for order, step in [(32, 1e-3), (48, 5e-4), (48, 2.5e-4)]:
            result = mode(load, column, order, step, guess)
            print(json.dumps(dict(name=name, **result)), flush=True)
            records.append(result)
        assert all(r['equilibrium_rhs_residual'] < 1e-6 and len(r['modes']) >= 1 for r in records)
        pressure_scan = []
        for pressure in [700., 1050., 1750., 2000.]:
            result = mode(load, column, 48, 2.5e-4, guess, pressure)
            assert result['equilibrium_rhs_residual'] < 1e-6 and result['modes']
            pressure_scan.append(result)
            print(json.dumps(dict(name=name, scan_pressure_Pa=pressure, modes=result['modes'])), flush=True)
        cases.append(dict(name=name, source_outlet_inertance=column,
                          impedance_checks=impedance_checks, records=records, pressure_scan=pressure_scan))
    paths = [Path(__file__), path, ROOT/'tools/continuous_phonation_reference.py',
             ROOT/'tools/continuous_outlet_column_reference.py', ROOT/'tools/slit_thermal_reference.py',
             ROOT/'tools/slit_viscous_reference.py', ROOT/'tools/planar_source_reference.py',
             ROOT/'tools/joint_contact_probe.cpp', ROOT/'bodies/reference.body',
             ROOT/'include/ducklips/sources.hpp', ROOT/'src/voice.cpp', ROOT/'src/body_config.cpp']
    report = dict(scope='Independent smooth symmetric source equations and frozen passive airway impedance; refinement at 1400 Pa and a declared pressure scan. '
                        'Local modes found from declared seeds, not an exhaustive stability proof or contact trajectory.',
                  source_parameters=expected_body, tension=1., rest_area_m2=3e-6,
                  cases=cases, hashes={str(p.relative_to(ROOT)): sha(p) for p in paths})
    (ROOT/'reports/loaded_source_onset.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
