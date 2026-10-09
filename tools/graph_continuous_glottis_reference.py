#!/usr/bin/env python3
"""Independent continuous and full-collocation checks of source/tract feedback."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from continuous_phonation_reference import Reference, YSCALE, VSCALE

ROOT = Path(__file__).resolve().parents[1]
EDGES = [(0, 1), (1, 2), (2, 0), (2, -1)]
B = np.array([[1., 0., -1., 0.], [-1., 1., 0., 0.], [0., -1., 1., 1.]])
C = np.array([1e-11, 1.4e-11, 8e-12])
G = np.array([.002, .003, .001, .0025])
NODAL = np.array([0., 1e-6, -2e-6])
SERIES = np.array([0., 30., 0., -10.])
PORT = 1
XSCALE = np.r_[np.full(3, 1e-8), np.full(4, .1)]
INITIAL = np.array([112.18465038080178, -3.138675468892969,
    3.6868231667933133e-6, 4.149026342131548e-6,
    .00017466261521028343, .00017466261521028343, .00017466261521028343, 0., 0.])
XINITIAL = np.r_[C * [120., -80., 40.], [.01, -.02, .015, .005]]
GAUSS = np.array([[.25, .25 - np.sqrt(3)/6], [.25 + np.sqrt(3)/6, .25]])
STAGES = GAUSS.sum(axis=1)

def metrics(t, moving):
    if not moving:
        return C, G
    return (C * (1 + .15*np.sin(2*np.pi*400*t + np.arange(3)*.5)),
            G * (1 + .2*np.sin(2*np.pi*300*t + np.arange(4)*.7)))

def matrix(c, g):
    a = np.zeros((7, 7))
    a[:3, 3:] = -B*g
    a[3:, :3] = B.T/c
    return a*XSCALE[None, :]/XSCALE[:, None]

def continuous(moving, order=32, tolerance=2e-9, feedback=True):
    ref = Reference(order)
    def rhs(t, y):
        physical = y[:9]*YSCALE
        x = y[9:]*XSCALE
        c, g = metrics(t, moving)
        derivative = ref.rhs(y[:9], 1400.)
        if feedback:
            state = ref.configuration(physical[2:4], physical[:2])
            load = np.array([0., 0., -x[PORT]/c[PORT], 0., 0.])
            derivative[4:] += np.linalg.solve(state['G']*VSCALE[:, None]*VSCALE[None, :], VSCALE*load)
        qdot = -B@(g*x[3:]) + NODAL
        qdot[PORT] += physical[6]
        pdot = B.T@(x[:3]/c) + SERIES
        return np.r_[derivative, np.r_[qdot, pdot]/XSCALE]
    sol = solve_ivp(rhs, (0., .002), np.r_[INITIAL/YSCALE, XINITIAL/XSCALE],
                    method='Radau', rtol=tolerance, atol=tolerance*.01, max_step=1e-5, dense_output=True)
    assert sol.success, sol.message
    return sol

def collocation(x, c, g, dt, flow):
    matrices = [matrix(c[s], g[s]) for s in range(2)]
    drives = NODAL.copy()
    drives[PORT] += flow
    source = np.r_[drives, SERIES]/XSCALE
    blocks = np.block([[np.eye(7)*(r == s) - dt*GAUSS[r, s]*matrices[s]
                        for s in range(2)] for r in range(2)])
    right = np.array([x/XSCALE + dt*STAGES[r]*source for r in range(2)])
    stage = np.linalg.solve(blocks, right.ravel()).reshape(2, 7)
    derivatives = np.array([matrices[s]@stage[s]+source for s in range(2)])
    endpoint = x + dt*.5*derivatives.sum(axis=0)*XSCALE
    physical = stage*XSCALE
    pressure = .5*(physical[:, :3]/c).sum(axis=0)
    flow = .5*(physical[:, 3:]*g).sum(axis=0)
    return endpoint, pressure, flow

def run_native(binary, out, moving, steps, reference):
    dt = .002/steps
    lines = [f'3 4 {steps} {PORT}', *[f'{i} {j}' for i, j in EDGES],
             ' '.join(map(str, INITIAL)), '1', ' '.join(map(str, XINITIAL))]
    for k in range(steps):
        c, g = zip(*(metrics((k+s)*dt, moving) for s in STAGES))
        row = np.r_[dt, 1400., 0., c[0], c[1], g[0], g[1], NODAL, SERIES]
        lines.append(' '.join(map(str, row)))
    name = f'{"moving" if moving else "fixed"}-{steps}'
    input_path = out/(name+'.input')
    input_path.write_text('\n'.join(lines)+'\n')
    result = subprocess.run([str(binary), '--fixture'], input=input_path.read_text(),
                            capture_output=True, text=True)
    (out/(name+'.stderr')).write_text(result.stderr)
    (out/(name+'.txt')).write_text(result.stdout)
    assert result.returncode == 0, result.stderr
    data = np.loadtxt(out/(name+'.txt'))
    assert data.shape == (steps, 42), data.shape
    graph = data[:, 9:16]
    ledger = data[:, 23:]
    oracle = XINITIAL.copy()
    maximum = dict(charge_equivalent_Pa=0., momentum_Pa_s=0., average_pressure_Pa=0., average_flow_m3_s=0.)
    max_mass = 0.
    source_mass = 0.
    net_transfer = 0.
    # Initial source mass from independent field quadrature.
    source_mass = Reference(48).configuration(INITIAL[2:4], INITIAL[:2])['mass'].sum()
    for k in range(steps):
        c, g = map(np.array, zip(*(metrics((k+s)*dt, moving) for s in STAGES)))
        oracle, pressure, flow = collocation(oracle, c, g, dt, ledger[k, 13])
        maximum['charge_equivalent_Pa'] = max(maximum['charge_equivalent_Pa'], float(np.max(np.abs(oracle[:3]-graph[k, :3])/C)))
        maximum['momentum_Pa_s'] = max(maximum['momentum_Pa_s'], float(np.max(np.abs(oracle[3:]-graph[k, 3:]))))
        maximum['average_pressure_Pa'] = max(maximum['average_pressure_Pa'], float(np.max(np.abs(pressure-data[k, 16:19]))))
        maximum['average_flow_m3_s'] = max(maximum['average_flow_m3_s'], float(np.max(np.abs(flow-data[k, 19:23]))))
        net_transfer += dt*(ledger[k, 17] + NODAL.sum() - data[k, 22])
        max_mass = max(max_mass, abs(ledger[k, 1:3].sum()+graph[k, :3].sum()-source_mass-XINITIAL[:3].sum()-net_transfer))
    assert maximum['charge_equivalent_Pa'] < 1e-8 and maximum['average_pressure_Pa'] < 1e-8, maximum
    assert maximum['momentum_Pa_s'] < 1e-12 and maximum['average_flow_m3_s'] < 1e-14, maximum
    assert max_mass < 1e-19, max_mass
    assert np.max(np.abs(ledger[:, 12])) < 1e-20
    assert np.max(np.abs(ledger[:, 18])) < 1e-18
    assert np.min(ledger[:, 7]) >= 0
    native = np.c_[data[:, :9]/YSCALE, graph/XSCALE]
    times = np.arange(1, steps+1)*dt
    error = np.abs(native-reference.sol(times).T)
    receipt = dict(steps=steps, dt_s=dt, endpoint_scaled_error=float(error[-1].max()),
        maximum_trajectory_scaled_error=float(error.max()),
        maximum_port_residual_J=float(np.abs(ledger[:, 12]).max()),
        maximum_source_step_residual_J=float(np.abs(ledger[:, 9]).max()),
        maximum_mass_residual_m3=max_mass,
        maximum_kinetic_loss_split_residual_J=float(np.abs(ledger[:, 18]).max()),
        outlet_loss_J=float(ledger[:, 7].sum()), full_collocation=maximum,
        final_state_scaled=native[-1].tolist(), pressure_range_Pa=[float(data[:, 16:19].min()), float(data[:, 16:19].max())])
    if not moving:
        # Independent source energy expression in the symmetric open subspace.
        from continuous_phonation_reference import STIFF, REST
        s = Reference(48).configuration(INITIAL[2:4], INITIAL[:2])
        p = Reference(48).fields(Reference(48).t, INITIAL[2:4], INITIAL[:2])[3]
        source_energy = INITIAL[:2]@s['mass'] - np.sum(s['w']*s['a']*p)
        source_energy += .5*INITIAL[4:]@s['G']@INITIAL[4:] + .5*(INITIAL[2:4]-REST)@STIFF@(INITIAL[2:4]-REST)
        e0 = source_energy + .5*np.sum(XINITIAL[:3]**2/C) + .5*np.sum(G*XINITIAL[3:]**2)
        energies = ledger[:, 0] + .5*np.sum(graph[:, :3]**2/C, axis=1) + .5*np.sum(G*graph[:, 3:]**2, axis=1)
        losses = np.cumsum(ledger[:, 5:9].sum(axis=1))
        work = np.cumsum(ledger[:, 3:5].sum(axis=1)+ledger[:, 11])
        receipt['maximum_joint_work_residual_J'] = float(np.abs(energies-e0+losses-work).max())
        assert receipt['maximum_joint_work_residual_J'] < 1e-13, receipt
    return receipt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', default='build/ducklips_graph_continuous_glottis_probe')
    parser.add_argument('--output', default='out/graph-continuous-glottis')
    parser.add_argument('--receipt', default='reports/graph_continuous_glottis_reference.json')
    parser.add_argument('--steps', nargs='+', type=int, default=[200, 400, 800, 1600, 3200, 6400])
    args = parser.parse_args()
    out = ROOT/args.output
    out.mkdir(parents=True, exist_ok=True)
    cases = []
    for moving in [False, True]:
        fine = continuous(moving, 48, 2e-10)
        coarse = continuous(moving, 32, 2e-9)
        sample_times = np.linspace(0., .002, 12801)
        np.savez_compressed(out/f'continuous-{moving}.npz', times=sample_times, fine=fine.sol(sample_times), coarse=coarse.sol(sample_times))
        reference_error = float(np.max(np.abs(fine.sol(np.linspace(0, .002, 201))-coarse.sol(np.linspace(0, .002, 201)))))
        assert reference_error < 2e-7, reference_error
        runs = []
        for steps in args.steps:
            run = run_native(ROOT/args.binary, out, moving, steps, fine)
            runs.append(run)
            (out/f'{"moving" if moving else "fixed"}-{steps}.receipt.json').write_text(json.dumps(run, indent=2)+'\n')
            print(json.dumps(dict(moving=moving, steps=steps, error=run['maximum_trajectory_scaled_error'], endpoint=run['endpoint_scaled_error'])), flush=True)
        ratios = [a['maximum_trajectory_scaled_error']/b['maximum_trajectory_scaled_error'] for a, b in zip(runs, runs[1:])]
        endpoint_ratios = [a['endpoint_scaled_error']/b['endpoint_scaled_error'] for a, b in zip(runs, runs[1:])]
        # Arbitrarily preloaded tract pressure excites a fast source transient.
        # Preserve coarse failures; require the final refinement to reach the
        # asymptotic regime rather than declaring conservation to be accuracy.
        assert ratios and 3.7 < ratios[-1] < 4.3, ratios
        assert endpoint_ratios and 3.7 < endpoint_ratios[-1] < 4.3, endpoint_ratios
        uncoupled = continuous(moving, 32, 2e-9, False)
        difference = np.abs(uncoupled.y[:9, -1]-fine.y[:9, -1])
        assert difference.max() > .01, difference
        cases.append(dict(moving_metrics=moving, reference_refinement_scaled_error=reference_error,
            native=runs, trajectory_refinement_ratios=ratios,
            endpoint_refinement_ratios=endpoint_ratios,
            removing_feedback_source_endpoint_scaled_difference=float(difference.max()),
            reference_final_state_scaled=fine.y[:, -1].tolist()))
    paths = [Path(__file__), ROOT/args.binary, ROOT/'tools/graph_continuous_glottis_probe.cpp',
             ROOT/'tools/continuous_phonation_reference.py', ROOT/'src/graph_continuous_glottis.cpp',
             ROOT/'include/ducklips/graph_continuous_glottis.hpp', ROOT/'src/continuous_glottis.cpp',
             ROOT/'include/ducklips/continuous_glottis.hpp', ROOT/'src/moving_graph.cpp',
             *out.glob('fixed-*.input'), *out.glob('fixed-*.txt'), *out.glob('fixed-*.stderr'),
             *out.glob('moving-*.input'), *out.glob('moving-*.txt'), *out.glob('moving-*.stderr'), *out.glob('continuous-*.npz')]
    report = dict(scope='Open full-span continuous source with fixed-topology linear acoustic graph; 2 ms nonlinear feedback, no speech or physiological claim',
        graph_edges=EDGES, source_node=PORT, reservoir_pressure_Pa=1400.,
        source_initial=INITIAL.tolist(), graph_initial=XINITIAL.tolist(),
        source_scale=YSCALE.tolist(), graph_scale=XSCALE.tolist(),
        reference='Separate continuous-field ODE plus charge/momentum graph ODE, Radau; independent full-state graph collocation', cases=cases,
        hashes={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/args.receipt).write_text(json.dumps(report, indent=2)+'\n')
    print('PASS implicit feedback, full graph collocation, mass/work balance, and second-order coupled convergence', flush=True)

if __name__ == '__main__':
    main()
