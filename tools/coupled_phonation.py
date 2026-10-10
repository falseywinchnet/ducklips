#!/usr/bin/env python3
"""Continue the fixed-tube coupled voice and export its physical radiation port.

This is a research trajectory driver, not the articulated speech renderer.
Every segment keeps the native gas, tissue, tract and loss-memory state.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.io import wavfile
from scipy.signal import resample_poly

ROOT = Path(__file__).resolve().parents[1]
SOURCE_SCALE = np.array([1e-6, 1000, 1000, *([1e-4] * 4),
                         *([.03571428571428571] * 4), *([1e-4] * 3)])


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def read_trace(path):
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    if len(rows) < 2 or rows[0]['time'] != 0 or 'restart' not in rows[-1]:
        raise ValueError(f'{path}: not a completed native segment')
    t = np.array([r['time'] for r in rows])
    if not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError(f'{path}: invalid sample times')
    for row in rows:
        if not (0 <= row['work_residual'] < 1e-12 and
                0 <= row['mass_residual'] < 1e-18 and min(row['gas_mass']) > 0):
            raise ValueError(f'{path}: invalid work, air or gas receipt')
        for j, port in [(0, 0), (1, 2)]:
            if row['area'][j] <= 0 and row['flow'][port] != 0:
                raise ValueError(f'{path}: flow through closed source port')
    for key in ['physical_loss', 'numerical_loss', 'closures', 'reopenings']:
        values = np.array([r[key] for r in rows])
        if not np.all(np.isfinite(values)) or np.any(values < 0) or np.any(np.diff(values) < 0):
            raise ValueError(f'{path}: invalid cumulative {key}')
    return rows


def restart_text(row):
    state = row['restart']
    if state['version'] != 1:
        raise ValueError('unknown native restart format')
    values = []
    for key, count in [('source', 14), ('charge', 8), ('momentum', 8),
                       ('radiation', 1), ('thermal', 192), ('viscous', 192)]:
        part = [state[key]] if key == 'radiation' else state[key]
        if len(part) != count or not np.all(np.isfinite(part)):
            raise ValueError(f'invalid restart field {key}')
        values.extend(part)
    return ' '.join(format(v, '.17g') for v in values) + '\n'


def check_join(before, after):
    a, b = before['restart'], after['restart']
    restart_text(before)
    restart_text(after)
    for key in ['charge', 'momentum', 'radiation', 'thermal', 'viscous']:
        if a[key] != b[key]:
            raise ValueError(f'continuation resets {key}')
    error = float(np.max(np.abs(np.array(a['source']) - b['source']) / SOURCE_SCALE))
    energy_error = abs(before['energy'] - after['energy'])
    mass_error = float(np.max(np.abs(np.array(before['gas_mass']) - after['gas_mass'])))
    if not (error < 1e-10 and energy_error < 1e-14 and mass_error < 1e-20):
        raise ValueError('source reconstruction changes the physical state')
    for key in ['outlet_flow', 'radiation_flow', 'radiation_resistance']:
        if before[key] != after[key]:
            raise ValueError(f'continuation changes radiation observation {key}')
    return dict(maximum_scaled_source_difference=error,
                energy_difference_J=energy_error, gas_mass_difference_m3=mass_error)


def collect(paths):
    rows, segments, joins = [], [], []
    offset = 0.
    for path in paths:
        data = read_trace(path)
        if rows:
            joins.append(check_join(previous, data[0]))
        # Adjacent segments share an endpoint: retain it only once.
        for row in data[1:] if rows else data:
            rows.append(dict(row, time=offset + row['time']))
        previous = data[-1]
        offset += previous['time']
        segments.append(dict(path=str(path.resolve()), sha256=sha(path),
                             duration_s=previous['time'], closures=previous['closures'],
                             reopenings=previous['reopenings'],
                             work_residual_J=previous['work_residual'],
                             mass_residual_m3=previous['mass_residual']))
    if not rows:
        raise ValueError('no completed trajectories')
    return rows, segments, joins


def power_wave(rows, rate):
    times = np.array([r['time'] for r in rows])
    native_rate = round((len(times) - 1) / times[-1])
    if native_rate <= 0 or not np.allclose(np.diff(times), 1 / native_rate, rtol=1e-8, atol=1e-14):
        raise ValueError('audio export requires a uniform physical observation clock')
    resistance = np.array([r['radiation_resistance'] for r in rows])
    difference = np.array([r['outlet_flow'] - r['radiation_flow'] for r in rows])
    if not np.all(np.isfinite(resistance)) or np.any(resistance < 0) or not np.all(np.isfinite(difference)):
        raise ValueError('invalid radiation state')
    # Boundary pressure P=R(Q-Hr), radiated power=P^2/R.
    # Signed power wave sqrt(R)(Q-Hr) has units sqrt(W), not microphone Pa.
    signal = np.sqrt(resistance) * difference
    divisor = math.gcd(rate, native_rate)
    # Keep the final endpoint as filter context but not an extra audio interval.
    audio = resample_poly(signal, rate // divisor, native_rate // divisor,
                          window=('kaiser', 8.6), padtype='constant')
    # Do not invent a final fractional playback interval beyond simulated time.
    count = math.floor(times[-1] * rate + 1e-8)
    if count < 1:
        raise ValueError('trajectory shorter than one output sample')
    audio = audio[:count]
    return times, signal, audio, native_rate


def export(paths, prefix, rate=48000, gain=1.):
    if rate <= 0 or not np.isfinite(gain) or gain <= 0:
        raise ValueError('positive output rate and fixed gain required')
    rows, segments, joins = collect(paths)
    times, signal, audio, native_rate = power_wave(rows, rate)
    played = audio * gain
    if not np.all(np.isfinite(played)) or np.max(np.abs(played)) >= 1:
        raise ValueError('playback gain would clip; choose a smaller fixed gain')
    prefix.parent.mkdir(parents=True, exist_ok=True)
    np.savez(prefix.with_suffix('.npz'), seconds=times, radiation_power_wave_sqrt_W=signal,
             resampled_power_wave_sqrt_W=audio, sample_rate=rate)
    wavfile.write(prefix.with_suffix('.wav'), rate, np.rint(played * 32767).astype(np.int16))
    report = dict(scope='Fixed-tube coupled voice diagnostic; not articulated speech or microphone pressure.',
                  duration_s=float(times[-1]), audio_duration_s=len(audio) / rate,
                  native_rate_hz=native_rate, sample_rate_hz=rate,
                  formula='sqrt(R) * (outlet_flow - radiation_flow)', units='sqrt(W)',
                  playback_gain_per_sqrt_W=gain, peak_before_gain=float(np.max(abs(audio))),
                  resampling='scipy.signal.resample_poly, Kaiser beta=8.6, zero exterior padding',
                  endpoint_policy='retain shared endpoint once; final endpoint is filter context only',
                  loops=0, time_stretch=False, peak_normalization=False,
                  sampled_radiation_energy_J=float(np.trapezoid(signal * signal, times)),
                  energy_scope='Trapezoidal observation integral, not the native integrated loss receipt.',
                  segments=segments, joins=joins, tool_sha256=sha(__file__),
                  outputs={suffix: sha(prefix.with_suffix(suffix)) for suffix in ['.npz', '.wav']})
    save(prefix.with_suffix('.json'), report)
    return report


def run(args):
    binary = args.binary.resolve()
    build = json.loads(args.build.read_text())
    binary_hash = sha(binary)
    if binary_hash != build['binary_sha256']:
        raise ValueError('binary differs from its build receipt')
    paths = [p.resolve() for p in args.seed]
    rows, segments, joins = collect(paths)
    if not np.isfinite(args.duration) or args.duration <= rows[-1]['time']:
        raise ValueError('target duration must extend the completed seed')
    for value in [args.segment, args.step, args.tolerance, args.contact_step]:
        if not np.isfinite(value) or value <= 0:
            raise ValueError('positive finite integration settings required')
    # A distinct directory prevents accidental overwrite of another live run.
    args.out.mkdir(parents=True, exist_ok=False)
    receipt = dict(scope='Fixed tube, 1400 Pa lung pressure, native physical continuation.',
                   binary=str(binary), binary_sha256=binary_hash, build=build,
                   driver_sha256=sha(__file__), target_duration_s=args.duration,
                   completed_duration_s=rows[-1]['time'], seeds=segments, segments=[], status='running')
    save(args.out / 'run.json', receipt)
    try:
        last = read_trace(paths[-1])[-1]
        while receipt['completed_duration_s'] < args.duration - 1e-12:
            # The native CLI uses ceil(duration / macro_step). Decimalize the
            # subtraction so a 20 us remainder does not become three samples
            # merely because the binary subtraction lands above 20 us.
            duration = float(format(min(args.segment, args.duration - receipt['completed_duration_s']), '.12g'))
            count = math.ceil(duration / args.step)
            if not math.isclose(duration / count, args.step, rel_tol=1e-8, abs_tol=1e-14):
                raise ValueError('segment duration must preserve the native observation clock')
            stem = args.out / f'segment-{len(receipt["segments"]):04d}'
            stem.with_suffix('.input').write_text(restart_text(last))
            command = [str(binary), 'phonation', format(duration, '.17g'), str(args.step),
                       str(args.tolerance), str(args.contact_step), '--initial-state']
            if sha(binary) != binary_hash:
                raise ValueError('executable changed during continuation')
            receipt['active_command'] = command
            receipt['active_trace'] = str(stem.with_suffix('.jsonl').resolve())
            save(args.out / 'run.json', receipt)
            start = time.monotonic()
            with stem.with_suffix('.input').open() as stdin, stem.with_suffix('.jsonl').open('w') as stdout, stem.with_suffix('.log').open('w') as stderr:
                child = subprocess.Popen(command, stdin=stdin, stdout=stdout, stderr=stderr)
                receipt['active_pid'] = child.pid
                save(args.out / 'run.json', receipt)
                code = child.wait()
            elapsed = time.monotonic() - start
            if code != 0:
                raise RuntimeError(f'native segment failed ({code}); retained {stem}')
            data = read_trace(stem.with_suffix('.jsonl'))
            join = check_join(last, data[0])
            if not math.isclose(data[-1]['time'], duration, rel_tol=0, abs_tol=1e-12):
                raise ValueError('native segment did not reach requested end')
            paths.append(stem.with_suffix('.jsonl').resolve())
            last = data[-1]
            receipt['completed_duration_s'] += duration
            receipt['segments'].append(dict(command=command, wall_seconds=elapsed, join=join,
                output_sha256=sha(paths[-1]), input_sha256=sha(stem.with_suffix('.input')),
                log_sha256=sha(stem.with_suffix('.log')), closures=last['closures'], reopenings=last['reopenings']))
            save(args.out / 'run.json', receipt)
            export(paths, args.out / 'radiation', args.rate, args.gain)
            print(json.dumps(dict(completed_s=receipt['completed_duration_s'], wall_s=elapsed)), flush=True)
        receipt['status'] = 'complete'
    except Exception as error:
        receipt['status'] = 'failed'
        receipt['error'] = str(error)
        raise
    finally:
        save(args.out / 'run.json', receipt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    audio = commands.add_parser('export')
    audio.add_argument('traces', nargs='+', type=Path)
    audio.add_argument('--out', required=True, type=Path)
    continuation = commands.add_parser('run')
    continuation.add_argument('--binary', required=True, type=Path)
    continuation.add_argument('--build', required=True, type=Path)
    continuation.add_argument('--seed', required=True, action='append', type=Path)
    continuation.add_argument('--out', required=True, type=Path)
    continuation.add_argument('--duration', required=True, type=float, help='total physical seconds including seeds')
    continuation.add_argument('--segment', type=float, default=.005)
    continuation.add_argument('--step', type=float, default=1e-5)
    continuation.add_argument('--tolerance', type=float, default=1e-4)
    continuation.add_argument('--contact-step', type=float, default=1e-7)
    for subparser in [audio, continuation]:
        subparser.add_argument('--rate', type=int, default=48000)
        subparser.add_argument('--gain', type=float, default=1., help='fixed PCM gain per sqrt(W), no peak normalization')
    args = parser.parse_args()
    if args.command == 'run':
        run(args)
    else:
        print(json.dumps(export(args.traces, args.out, args.rate, args.gain), indent=2))


if __name__ == '__main__':
    main()
