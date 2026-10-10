#!/usr/bin/env python3
"""Export completed actual-mouth segments from their physical radiation port."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import wavfile

from coupled_phonation import power_wave, save, sha
from joint_accuracy_moving_reference import ROOT, audit

# The actual-mouth probe uses AirProperties{} and NetworkLosses{} at this
# radiation boundary. These are declared model constants, not microphone SPL.
RHO, SOUND_SPEED, RADIATION_FACTOR = 1.14, 350., 1.4
COEFFICIENT = RADIATION_FACTOR * RHO / (4 * np.pi * SOUND_SPEED)
OUTLET_BRANCH = 6  # Declared two-region / four-lip / eight-throat adapter.


def join(before, after):
    for key in ['volume', 'compliance', 'conductance', 'charge', 'momentum',
                'radiation_state', 'thermal_sites', 'viscous_sites']:
        if before[key] != after[key]:
            raise ValueError(f'physical continuation resets {key}')
    a, b = before['continuous_source'], after['continuous_source']
    mass = float(np.max(abs(np.array(a['mass']) - b['mass'])))
    momentum = float(np.max(abs(np.array(a['canonical_momentum']) - b['canonical_momentum'])))
    energy = abs(a['energy'] - b['energy'])
    if not (mass < 1e-20 and momentum < 1e-14 and energy < 1e-13):
        raise ValueError('source changes at continuation')
    return dict(mass_m3=mass, canonical_momentum=momentum, energy_J=energy)


def observation(row, time):
    if len(row['charge']) != 14 or len(row['flow']) != 15:
        raise ValueError('unsupported actual-mouth topology')
    radiation = row['radiation_state']
    h = radiation['inverse_inertance']
    q = row['flow'][OUTLET_BRANCH]
    memory = h * radiation['momentum']
    if not np.isfinite(h) or h < 0:
        raise ValueError('invalid radiation inertia')
    if h == 0 and (q != 0 or radiation['momentum'] != 0):
        raise ValueError('flow or exterior memory at a sealed outlet')
    r = 1 / (COEFFICIENT * h * h) if h > 0 else 0.
    source = row['continuous_source']
    if min(source['mass']) <= 0:
        raise ValueError('nonpositive source gas')
    for j, port in [(0, 0), (1, 2)]:
        if source['area'][j] <= 0 and source['flow'][port] != 0:
            raise ValueError('flow through a sealed source port')
    return dict(time=time, outlet_flow=q, radiation_flow=memory, radiation_resistance=r,
                source_area=source['area'], source_flow=source['flow'],
                lower_position=source['lower_position'], upper_position=source['upper_position'])


def collect(run_path):
    manifest_bytes = run_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    for name in ['include/ducklips/physics.hpp', 'include/ducklips/acoustic_network.hpp']:
        if sha(ROOT / name) != manifest['source_hashes'][name]:
            raise ValueError(f'radiation defaults differ from the recorded engine: {name}')
    binary = run_path.parent / 'probe'
    if sha(binary) != manifest['binary_sha256']:
        raise ValueError('executable differs from trajectory provenance')
    if manifest.get('seed_run'):
        seed = manifest['seed_run']
        seed_path = Path(seed['path'])
        if sha(seed_path) != seed['sha256']:
            raise ValueError('completed seed manifest changed')
        rows, segments, joins, previous, manifests = collect(seed_path)
        offset = rows[-1]['time']
    else:
        rows, segments, joins, previous, manifests, offset = [], [], [], None, [], 0.
    for receipt in manifest['segments']:
        paths = {Path(name).suffix: ROOT / name for name in receipt['files']}
        for name, expected in receipt['files'].items():
            if sha(ROOT / name) != expected:
                raise ValueError(f'completed segment changed: {name}')
        data = json.loads(paths['.json'].read_text())
        checks = audit(data, 1e-4)
        if previous is not None:
            joins.append(join(previous, data['initial']))
        else:
            rows.append(observation(data['initial'], 0.))
        local = 0.
        for row in data['steps']:
            local += row['dt']
            rows.append(observation(row, offset + local))
        offset += local
        previous = data['steps'][-1]
        segments.append(dict(path=str(paths['.json']), sha256=sha(paths['.json']),
                             duration_s=local, checks=checks,
                             source_closures=sum(r['accuracy']['source_closures'] for r in data['steps']),
                             source_reopenings=sum(r['accuracy']['source_reopenings'] for r in data['steps'])))
    if not rows or abs(offset - manifest['completed_duration_s']) > 1e-12:
        raise ValueError('completed duration and native segments disagree')
    # Save the exact manifest bytes read, even if its producer is still active.
    import hashlib
    manifests.append(dict(path=str(run_path.resolve()), sha256=hashlib.sha256(manifest_bytes).hexdigest(),
                          status_when_read=manifest['status'], binary_sha256=manifest['binary_sha256'],
                          completed_duration_s=manifest['completed_duration_s']))
    return rows, segments, joins, previous, manifests


def export(run_path, prefix, gain=1., peak=None):
    if not np.isfinite(gain) or gain <= 0:
        raise ValueError('positive fixed playback gain required')
    manifest = json.loads(run_path.read_bytes())
    word = manifest.get('word')
    if word:
        target = manifest.get('target_duration_s', float('nan'))
        completed = manifest.get('completed_duration_s', float('nan'))
        if (manifest.get('status') != 'complete' or not np.isfinite(target)
                or target <= 0 or not np.isfinite(completed)
                or abs(completed-target) > 1e-12):
            raise ValueError('word audio requires the complete requested utterance')
        if any(prefix.with_suffix(s).exists() for s in ['.npz', '.wav', '.json']):
            raise ValueError('use a new word audio prefix; preserve previous renders')
    if peak is not None and (not np.isfinite(peak) or not 0 < peak < 1 or gain != 1.):
        raise ValueError('peak must be between zero and one, without an explicit gain')
    rows, segments, joins, _, manifests = collect(run_path)
    times, signal, audio, native_rate = power_wave(rows, 48000)
    if word and abs(times[-1]-target) > 1e-12:
        raise ValueError('word trajectory does not cover the requested utterance')
    if peak is not None:
        amplitude = float(np.max(abs(audio)))
        if not np.isfinite(amplitude) or amplitude <= 0:
            raise ValueError('cannot set playback level for silent or invalid audio')
        gain = peak / amplitude
    played = gain * audio
    if not np.all(np.isfinite(played)) or np.max(abs(played)) >= 1:
        raise ValueError('fixed playback gain clips')
    prefix.parent.mkdir(parents=True, exist_ok=True)
    arrays = {key: np.array([r[key] for r in rows]) for key in
              ['source_area', 'source_flow', 'lower_position', 'upper_position',
               'outlet_flow', 'radiation_flow', 'radiation_resistance']}
    np.savez(prefix.with_suffix('.npz'), seconds=times, radiation_power_wave_sqrt_W=signal,
             resampled_power_wave_sqrt_W=audio, sample_rate=48000, **arrays)
    wavfile.write(prefix.with_suffix('.wav'), 48000, np.rint(32767 * played).astype(np.int16))
    report = dict(scope='Current continuous-source actual-mouth diagnostic; '
                        'two oral pressure regions, not a NATO word or calibrated microphone recording.',
                  duration_s=float(times[-1]), audio_duration_s=len(audio) / 48000,
                  sample_rate_hz=48000, observation_rate_hz=native_rate,
                  formula='sqrt(R)*(Q_outlet-H*p_radiation); R=1/(K*H^2)',
                  radiation_coefficient=COEFFICIENT, outlet_branch=OUTLET_BRANCH,
                  signal_units='sqrt(W)', fixed_playback_gain_per_sqrt_W=gain,
                  peak_before_gain=float(np.max(abs(audio))),
                  resampling='polyphase Kaiser beta 8.6, zero exterior padding; no time stretching or looping',
                  sampled_radiation_energy_J=float(np.trapezoid(signal**2, times)),
                  energy_scope='Observation integral, not the native integrated radiation-loss receipt.',
                  manifests=manifests, segments=segments, joins=joins,
                  files={str(p.relative_to(ROOT)): sha(p) for p in
                         [Path(__file__), ROOT / 'tools/coupled_phonation.py',
                          ROOT / 'include/ducklips/physics.hpp', ROOT / 'include/ducklips/acoustic_network.hpp',
                          prefix.with_suffix('.npz'), prefix.with_suffix('.wav')]})
    if word:
        report.update(scope=manifest['scope'], word=word,
                      input_start_s=manifest['input_start_s'], input_end_s=manifest['input_end_s'],
                      applied_controls=manifest['applied_controls'],
                      unapplied_controls=manifest['unapplied_controls'],
                      body_sha256=manifest['body_sha256'],
                      perceptual_status='Newly rendered utterance; intelligibility unassessed.')
    report['playback_level_selection'] = ('single gain from complete waveform peak'
                                        if peak is not None else 'explicit fixed gain')
    save(prefix.with_suffix('.json'), report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    levels = parser.add_mutually_exclusive_group()
    levels.add_argument('--gain', type=float, default=1.)
    levels.add_argument('--peak', type=float, help='set one fixed gain from the complete waveform peak')
    args = parser.parse_args()
    report = export(args.run.resolve(), args.out.resolve(), args.gain, args.peak)
    print(json.dumps({k: report[k] for k in ['duration_s', 'peak_before_gain', 'sampled_radiation_energy_J']}))


if __name__ == '__main__':
    main()
