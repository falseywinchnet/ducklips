#!/usr/bin/env python3
"""Exact frozen surface path elimination; compare complex responses to impulses.

This is a frequency-domain reference, not a new speech renderer. No paths,
cavities, transverse channels, or passive filters are omitted.
"""
import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


class FrozenPorts:
    def __init__(self, snapshot):
        self.snapshot = snapshot
        self.fs = snapshot['sample_rate']
        paths = snapshot['paths']
        self.nodes = np.array([[p[0], p[1]] for p in paths], dtype=int)
        self.y = np.array([p[2] for p in paths])
        self.delay = np.array([p[3] for p in paths])
        self.b = np.array([p[4] for p in paths])
        self.a = np.array([p[5] for p in paths])
        if not np.array_equal(self.nodes[::2], self.nodes[1::2, ::-1]):
            raise ValueError('Expected adjacent reciprocal paths')
        if not np.allclose(self.y[::2], self.y[1::2], rtol=1e-14, atol=0):
            raise ValueError('Unequal reciprocal admittances')
        self.ports = snapshot['ports']
        self.port_nodes = np.array([p[0] for p in self.ports], dtype=int)
        self.port_y = np.array([p[1] for p in self.ports])
        self.port_b = np.array([p[3] for p in self.ports])
        self.port_a = np.array([p[4] for p in self.ports])
        active = np.zeros(snapshot['pieces'], dtype=bool)
        active[self.nodes[self.y > 0].ravel()] = True
        active[self.port_nodes[self.port_y > 0]] = True
        self.active = np.flatnonzero(active)
        self.remap = np.full(snapshot['pieces'], -1, dtype=int)
        self.remap[self.active] = np.arange(len(self.active))
        self.pair_mask = self.y[::2] > 0
        self.pair_nodes = self.nodes[::2][self.pair_mask]
        self.pair_y = self.y[::2][self.pair_mask]

    @staticmethod
    def transfer(b, a, q):
        powers = np.array([1, q, q*q])
        return (b @ powers) / (a @ powers)

    def response(self, frequency, damping=0.0):
        q = np.exp(-(damping + 2j*np.pi*frequency) / self.fs)
        h = q**self.delay * self.transfer(self.b, self.a, q)
        forward, backward = h[::2][self.pair_mask], h[1::2][self.pair_mask]
        product = forward * backward
        divisor = 1 - product
        if np.min(np.abs(divisor)) < 1e-12:
            raise ValueError('Path-elimination pole; use a positive Laplace damping contour')
        diagonal = self.pair_y * (1 + product) / divisor
        ab = -2*self.pair_y*backward/divisor
        ba = -2*self.pair_y*forward/divisor
        a, b = self.pair_nodes.T
        rows = np.concatenate([a, b, a, b, self.port_nodes])
        cols = np.concatenate([a, b, b, a, self.port_nodes])
        reflection = self.transfer(self.port_b, self.port_a, q)
        shunt = self.port_y * (1-reflection)/(1+reflection)
        values = np.concatenate([diagonal, diagonal, ab, ba, shunt])
        matrix = coo_matrix((values, (self.remap[rows], self.remap[cols])),
                            shape=(len(self.active), len(self.active))).tocsc()
        rhs = np.zeros(len(self.active), dtype=complex)
        source = self.remap[self.snapshot['source_piece']]
        if source < 0:
            raise ValueError('Source is on an unconnected piece')
        rhs[source] = 1
        p = spsolve(matrix, rhs)
        if not np.all(np.isfinite(p)):
            raise RuntimeError('Nonfinite frozen pressure solve')
        pressure = p[self.remap[self.snapshot['pressure_pieces']]]
        radiation = sum(shunt[i]*p[self.remap[port[0]]] for i, port in enumerate(self.ports) if port[2])
        microphone = self.fs*(1-q)*radiation
        # Every two-terminal path and each one-terminal port must dissipate
        # nonnegative real power. Summation then preserves this for the mesh.
        off_hermitian = (ab + np.conj(ba))/2
        least_path_power = diagonal.real - np.abs(off_hermitian)
        scale = max(float(np.max(np.abs(diagonal))), float(np.max(np.abs(shunt))), 1e-30)
        power_margin = min(float(least_path_power.min()), float(shunt.real.min())) / scale
        residual = np.linalg.norm(matrix @ p-rhs) / max(np.linalg.norm(rhs), 1e-30)
        if not (np.isfinite(microphone) and np.isfinite(power_margin) and np.isfinite(residual)):
            raise RuntimeError('Nonfinite frozen response or validation metric')
        return np.append(pressure, microphone), power_margin, residual


def analyze(prefix, damping=40.0):
    receipt_path = prefix.with_suffix('.json')
    snapshot = json.loads(receipt_path.read_text())
    data = np.fromfile(prefix.with_suffix('.f64'), dtype=np.float64).reshape(snapshot['samples'], snapshot['channels'])
    if not np.all(np.isfinite(data)):
        raise ValueError('Nonfinite time-domain reference')
    model = FrozenPorts(snapshot)
    # Exponential window evaluates exactly the same transfer function at
    # z=exp((alpha+j*w)/fs). It makes finite-tail truncation explicit/bounded.
    window = np.exp(-damping*np.arange(len(data))/model.fs)
    spectrum = np.fft.rfft(data*window[:,None], axis=0) / snapshot['impulse_m3_per_s']
    bins = np.unique(np.round(np.geomspace(50, 12000, 220)*len(data)/model.fs).astype(int))
    frequencies = bins*model.fs/len(data)
    started = time.perf_counter()
    direct, margins, residuals = [], [], []
    for frequency in frequencies:
        result, margin, residual = model.response(frequency, damping)
        direct.append(result); margins.append(margin); residuals.append(residual)
    elapsed = time.perf_counter()-started
    direct, observed = np.array(direct), spectrum[bins]
    relative = np.abs(direct-observed)/np.maximum(np.abs(direct), np.max(np.abs(direct),axis=0)*1e-6)
    channel_reports = []
    for channel in range(snapshot['channels']):
        label = ('pressure_piece_'+str(snapshot['pressure_pieces'][channel])
                 if channel < len(snapshot['pressure_pieces']) else 'radiated_flow_derivative')
        channel_reports.append({'channel':label, 'max_relative_complex_error':float(relative[:,channel].max()),
                                'rms_relative_complex_error':float(np.sqrt(np.mean(relative[:,channel]**2)))})
    report = {'fixture':snapshot['fixture'], 'source':snapshot['source'],
              'snapshot_sha256':hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
              'response_sha256':hashlib.sha256(prefix.with_suffix('.f64').read_bytes()).hexdigest(),
              'pieces':snapshot['pieces'], 'active_pressure_unknowns':len(model.active),
              'allocated_delay_samples':snapshot['allocated_delay_samples'],
              'directed_paths':len(model.y), 'laplace_damping_per_second':damping,
              'frequency_count':len(frequencies), 'frequency_range_hz':[float(frequencies.min()),float(frequencies.max())],
              'frequency_solve_seconds':elapsed, 'impulse_render_seconds':snapshot['render_seconds'],
              'impulse_duration':len(data)/model.fs, 'instrumented_render_speed':len(data)/model.fs/snapshot['render_seconds'],
              'gather_seconds':snapshot['gather_seconds'], 'scatter_seconds':snapshot['scatter_seconds'],
              'minimum_scaled_passive_power_margin':min(margins), 'max_linear_residual':max(residuals),
              'channels':channel_reports}
    np.savez_compressed(str(prefix)+'-ports.npz', hz=frequencies, direct=direct, impulse=observed,
                        relative_complex_error=relative)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix', type=Path, nargs='+')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--binary', type=Path, default=Path('build/surface/surface_response'))
    args = parser.parse_args()
    reports = [analyze(prefix) for prefix in args.prefix]
    reciprocity = []
    # Independent drive locations must give the same reversed transfer in a
    # reciprocal frozen geometry. Only compare identical exported networks.
    for index, first in enumerate(args.prefix):
        one = json.loads(first.with_suffix('.json').read_text())
        for second in args.prefix[index+1:]:
            two = json.loads(second.with_suffix('.json').read_text())
            if (one['source_piece'] == two['source_piece'] or one['paths'] != two['paths']
                    or one['ports'] != two['ports'] or one['sample_rate'] != two['sample_rate']
                    or two['source_piece'] not in one['pressure_pieces']
                    or one['source_piece'] not in two['pressure_pieces']):
                continue
            with np.load(str(first)+'-ports.npz') as a, np.load(str(second)+'-ports.npz') as b:
                if not np.array_equal(a['hz'], b['hz']):
                    continue
                forward = a['direct'][:, one['pressure_pieces'].index(two['source_piece'])]
                reverse = b['direct'][:, two['pressure_pieces'].index(one['source_piece'])]
                error = float(np.max(np.abs(forward-reverse)/np.maximum(np.abs(forward),1e-30)))
                reciprocity.append({'first':str(first), 'second':str(second),
                                    'max_relative_complex_error':error})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    receipt = {'schema': 1, 'platform': platform.platform(),
               'processor': platform.machine(), 'numpy_version': np.__version__,
               'probe_binary': str(args.binary),
               'probe_binary_sha256': hashlib.sha256(args.binary.read_bytes()).hexdigest(),
               'git_head_at_measurement': subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
               'working_tree_modified': bool(subprocess.check_output(['git','status','--porcelain'],text=True).strip()),
               'observations': reports, 'reciprocity_checks': reciprocity}
    args.out.write_text(json.dumps(receipt, indent=2)+'\n')
    for report in reports:
        worst = max(c['max_relative_complex_error'] for c in report['channels'])
        print(f"{report['fixture']} / {report['source']}: worst complex error {worst:.3g}, "
              f"{report['instrumented_render_speed']:.3f}x impulse render, "
              f"{report['frequency_solve_seconds']:.3f} s frozen frequency solve")
    if any(r['max_relative_complex_error'] > 1e-9 for r in reciprocity) or any(r['minimum_scaled_passive_power_margin'] < -1e-10
           or r['max_linear_residual'] > 1e-9
           or any(c['max_relative_complex_error'] > 1e-5 for c in r['channels']) for r in reports):
        raise SystemExit('Frozen-response verification gate failed; inspect errors, do not accept reduction')

if __name__ == '__main__':
    main()
