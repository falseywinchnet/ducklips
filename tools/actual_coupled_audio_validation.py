#!/usr/bin/env python3
"""Check actual-mouth export, native field pixels and measured mode decay."""
import copy
import json
from pathlib import Path
import sys

import numpy as np
from PIL import Image
from scipy.io import wavfile
from scipy.optimize import least_squares

from actual_coupled_audio import ROOT, collect, join, observation
from coupled_phonation import sha
from native_fields import enhanced_rgba


def must_reject(function):
    try:function()
    except ValueError:return
    raise AssertionError('invalid physical observation accepted')


def main():
    out = ROOT / 'out/actual-coupled-phonation-50ms'
    run_path = out / 'run.json'
    run = json.loads(run_path.read_text())
    assert run['status'] == 'complete' and abs(run['completed_duration_s']-.05) < 1e-12
    rows, segments, joins, last, manifests = collect(run_path)
    assert len(rows) == 5001 and len(segments) == 50 and len(joins) == 49
    arrays = np.load(out / 'radiation.npz')
    audio_report = json.loads((out / 'radiation.json').read_text())
    assert audio_report['manifests'] == manifests
    for name, expected in audio_report['files'].items():assert sha(ROOT / name) == expected, name
    times = arrays['seconds'];signal = arrays['radiation_power_wave_sqrt_W']
    expected = np.array([np.sqrt(r['radiation_resistance'])*(r['outlet_flow']-r['radiation_flow']) for r in rows])
    assert np.array_equal(signal, expected)
    fs, pcm = wavfile.read(out / 'radiation.wav')
    assert fs == 48000 and len(pcm) == 2400
    assert np.array_equal(pcm, np.rint(32767*arrays['resampled_power_wave_sqrt_W']).astype(np.int16))
    assert not any(s['source_closures'] or s['source_reopenings'] for s in segments)
    bad = copy.deepcopy(last);bad['thermal_sites'][0]['charge'][0] += 1e-15
    must_reject(lambda: join(last, bad))
    sealed = copy.deepcopy(last)
    sealed['radiation_state'].update(inverse_inertance=0., momentum=0.)
    sealed['flow'][6] = 0.
    assert observation(sealed, 0.)['radiation_resistance'] == 0
    sealed['flow'][6] = 1e-6
    must_reject(lambda: observation(sealed, 0.))

    visual = json.loads((out / 'enhanced.json').read_text())
    field = np.load(out / 'enhanced.npz')
    png = np.asarray(Image.open(out / 'enhanced.png').convert('RGBA'))
    pixels = np.repeat(enhanced_rgba(field['enhanced'], visual['display_reference'], visual['db_range']), 4, axis=1)
    assert np.array_equal(png, pixels)
    assert field['enhanced'].shape == (684, 51) and png.shape == (684, 204, 4)
    assert visual['start_sample'] == 0 and visual['end_sample'] == 2400 and visual['hop'] == 48
    assert visual['source_sha256'] == sha(out / 'radiation.wav')

    # Fit the native tissue observations independently of the onset solver.
    # Two damped sinusoidal components plus an offset; repeat the fit windows
    # rather than using a chosen peak pair to manufacture an envelope slope.
    fits = []
    areas = arrays['source_area']/1e-6
    for start in [.015, .020, .025, .030]:
        mask = times >= start;clock = times[mask]-start;target = areas[mask]
        def residual(x, details=False):
            a, f, b, g = x
            basis = np.column_stack([np.ones(len(clock)), np.exp(a*clock)*np.cos(2*np.pi*f*clock),
                                     np.exp(a*clock)*np.sin(2*np.pi*f*clock),
                                     np.exp(b*clock)*np.cos(2*np.pi*g*clock),
                                     np.exp(b*clock)*np.sin(2*np.pi*g*clock)])
            coefficients = np.linalg.lstsq(basis, target, rcond=None)[0]
            error = basis@coefficients-target
            return (error, coefficients) if details else error.ravel()
        solved = least_squares(residual, [-25., 320., -60., 155.],
                               bounds=([-150, 250, -300, 100], [50, 400, 0, 200]),
                               xtol=1e-12, ftol=1e-12, gtol=1e-12)
        error, coefficients = residual(solved.x, True)
        fits.append(dict(first_time_s=start, last_time_s=float(times[-1]),
                         growth_s=float(solved.x[0]), frequency_Hz=float(solved.x[1]),
                         secondary_growth_s=float(solved.x[2]), secondary_frequency_Hz=float(solved.x[3]),
                         RMS_area_residual_m2=float(np.sqrt(np.mean(error**2))*1e-6),
                         fitted_mean_area_m2=(coefficients[0]*1e-6).tolist()))
    onset_path = ROOT / 'reports/loaded_source_onset.json'
    onset = json.loads(onset_path.read_text())
    for name, expected in onset['hashes'].items():assert sha(ROOT / name) == expected, name
    reference = max(onset['cases'][0]['records'][-1]['modes'], key=lambda r:r['frequency_Hz'])
    assert all(abs(r['growth_s']-reference['growth_s']) < 1 and
               abs(r['frequency_Hz']-reference['frequency_Hz']) < 1 for r in fits)
    paths = [Path(__file__), run_path, Path(run['seed_run']['path']), onset_path,
             ROOT/'tools/actual_coupled_phonation.py', ROOT/'tools/enhanced_specgram/render.py',
             ROOT/'tools/enhanced_specgram/README.md',
             *[out/f'radiation.{ext}' for ext in ['json','npz','wav']],
             *[out/f'enhanced.{ext}' for ext in ['json','npz','png']]]
    report = dict(scope='Completed 50 ms actual-mouth diagnostic and local decay agreement; no sustained speech claim.',
                  segments=50, unique_observations=5001, sample_rate_hz=48000, audio_samples=2400,
                  maximum_work_residual_J=max(s['checks']['work_residual_J'] for s in segments),
                  maximum_air_residual_m3=max(s['checks']['mass_residual_m3'] for s in segments),
                  source_closures=0, source_reopenings=0, enhanced_pixels_exact=True,
                  enhanced_shape=list(png.shape), native_decay_fits=fits, independent_local_mode=reference,
                  hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})
    (ROOT/'reports/actual_coupled_audio_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS 50 ms physical continuation, radiation audio, full enhanced pixels and native mode decay')


if __name__ == '__main__':main()
