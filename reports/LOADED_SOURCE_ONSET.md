# Why the actual mouth does not sustain the tube's voice

The current source completes 50 ms in the actual-mouth control at 1400 Pa,
including the initial held pose, one millisecond of actual body/lip motion,
and physical continuation of the final pose. Its oscillation decays rather
than reaching a sustained voice. No source opening closes in this trajectory.
The fixed-tube reference at the same pressure has already completed repeated
closures. A successful tube source is therefore not sufficient evidence of
phonation in the actual mouth.

An independent local calculation explains the observed difference. It uses
the previously checked continuous gas/tissue equations in the symmetric open
configuration, the declared tissue parameters, the actual glottis-to-throat
half-column inertia, and the airway's frequency-dependent pressure response.
It does not change the native engine's forces, damping or source waveform.

## Pressure feedback and the local modes

The frozen passive airway gives a source-port impedance `P(s) = Z(s) Q(s)`.
Its response includes branch air inertia, each thermal storage mode, each
active viscous profile mode and the radiation boundary. Eliminating these
memories analytically leaves a nodal pressure solve. A separate sparse solve
retains every active memory explicitly: 1,158 unknowns for the actual mouth
and 401 for the fixed tube. At five frequencies from 20 Hz to 8 kHz, maximum
relative impedance discrepancies are 6.23e-15 and 3.67e-14 respectively.
The checked real-frequency impedances have positive real parts.

After finding a steady open source against the airway's DC resistance, the
independent source equations are linearized about that state. If `A` is the
source matrix, `b` applies port pressure and `c` observes outlet flow, coupled
modes satisfy

```
1 - Z(s) c (s I - A)^(-1) b = 0
```

The attached half-column is included in the source inertia once, as in the
native engine. The DC solve contracts zero-resistance edges; it does not add
an artificial small resistance to fix a circulation ambiguity.

| Configuration at 1400 Pa | Mode frequency | Amplitude growth rate |
|---|---:|---:|
| Actual mouth, independent local mode | 319.551 Hz | −25.725 s⁻¹ |
| Actual mouth, native fits over four late-time windows | 319.742–319.932 Hz | −25.968 to −25.641 s⁻¹ |
| Fixed tube, independent local mode | 281.846 Hz | +297.598 s⁻¹ |

Negative growth means decay; positive growth means a small perturbation grows.
The native fits use both fold-opening traces, two damped sinusoidal components
and a fitted offset. Fit windows begin at 15, 20, 25 and 30 ms and all end at
50 ms. This agreement is evidence for the observed local mode, not a proof of
an exhaustive spectrum, global stability, contact behavior or voice realism.
Changing quadrature from 32 to 48 points and halving the source derivative
step retains the mode's sign and changes its reported frequency by less than
0.000022 Hz.

Mean airway load alone does not explain the difference. The computed steady
source-port pressures are only about 1.825 Pa for the actual mouth and
0.719 Pa for the tube. Replacing the dynamic airway by its DC resistance gives
different local modes: in both configurations the growing mode is near 149 Hz.
That approximation misses the actual simulation's observed 320 Hz decay.
The phase and frequency dependence of the pressure feedback matter.

The [onset receipt](loaded_source_onset.json) includes a declared pressure
scan. In the actual-mouth model, the tracked higher-frequency mode changes
from decay at 1400 Pa to growth at 1750 Pa: +43.014 s⁻¹ near 347.035 Hz.
An 80 ms native case at 1750 Pa has been started to test that prediction.
It is not yet claimed complete or voiced. These are model pressures and
threshold behavior, not a measured human physiological calibration.

## Radiation audio from the actual mouth

`tools/actual_coupled_audio.py` now exports completed actual-mouth segments,
including a chain of continued runs. It verifies segment hashes, executable
identity, unchanged geometry and memory values at joins, source mass and
momentum reconstruction, work/air receipts, and the recorded radiation
defaults. A duplicated shared endpoint is retained once.

The output signal is the signed radiation power wave
`sqrt(R) * (Q_outlet - H * p_radiation)`, with units `sqrt(W)`.
It comes from the actual mouth's outlet flow and radiation state. The WAV
uses a declared fixed gain, 48 kHz anti-alias resampling and PCM16 storage.
There is no peak normalization, looped cycle or stretched time. Microphone
distance, directivity and calibrated sound-pressure level remain unspecified.

The completed diagnostic is 50 ms: 5,001 unique native observations become
2,400 playback samples. It is a short decaying onset, not a word or a sustained
voice recording. Across its 50 completed segments, the largest checked work
residual is 4.07e-16 J and air residual is 3.03e-20 m³. The
[audio and trajectory validation](actual_coupled_audio_validation.json) binds
the native records, export, mode fits and visual artifact.

The full enhanced field uses the original apertures, frequency rows and magma
colors, with explicitly denser 48-sample time centers. Every computed field
cell survives: 684 frequency rows and 51 time centers produce a 204×684 RGBA
PNG, with four image columns per time center and no interpolation. Pixel
equality is checked. No ordinary-STFT comparison panel is presented. Because
the longest analysis window is longer than this 50 ms diagnostic, padding
affects the field substantially; it cannot establish steady harmonic quality.

```sh
python tools/actual_coupled_audio.py \
  out/actual-coupled-phonation-50ms/run.json \
  --out out/actual-coupled-phonation-50ms/radiation
python tools/enhanced_specgram/render.py \
  out/actual-coupled-phonation-50ms/radiation.wav --max-hz 8000 --hop 48 \
  --out out/actual-coupled-phonation-50ms/enhanced
python tools/loaded_source_onset.py
python tools/actual_coupled_audio_validation.py
```

Use the scientific Python environment. The actual-mouth geometry remains a
research discretization with two oral pressure regions. Speech articulation,
realism, body calibration and the original efficiency requirement remain open.
