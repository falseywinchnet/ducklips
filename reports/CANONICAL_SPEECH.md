# Canonical motion connected to speech

The `--canonical-motion` path now couples the full canonical moving-network
state to the existing physical vocal folds, aerodynamic cavity, articulators,
body configuration, nudge queue and microphone path. It does not run the old
modal resonator, modal eigensolver or geometry reprojection. It also does not
use the experimental reduced `advance_gauss` kernel. The full canonical
network is the source-coupled physical reference from which that reduction
can subsequently be assessed. Performance optimization is deferred at the
user's explicit request; the 10x final requirement remains unmet.

## Source coupling and the instability it exposed

The physical states are nodal acoustic charge `q` and branch momentum `lambda`:

```
q'      = -B G(t) lambda + nodal_flow
lambda' =  B^T C(t)^-1 q + series_pressure
```

Every branch retains its momentum when geometry or source location changes.
No gradient-only projection discards the source history. `MovingTreeReference`
now accepts held nodal flow and series pressure in its two-stage Gauss solve.
An independent forced-oscillator check verifies the affine update.

The renderer evolves two solutions for each audio sample: the free response
including aspiration/frication and a zero-initial-state response to unit
volume flow at the glottis. The orifice equation determines the actual flow;
superposition injects it into the acoustic and loss-memory states. Pressure
feedback uses the Gauss integral of pressure throughout the sample. With a
held flow, its acoustic work is exactly `dt * flow * average_pressure`.

The first integration reused the old endpoint-pressure trapezoid. That was
not a passive coupling for the full network, which contains modes above
output Nyquist. It produced runaway pressure and nonfinite output. Replacing
that approximation with the work-conjugate pressure integral stabilized Alfa
without state suppression or output normalization. Tests check the source
work identity at angular frequencies through 400,000 rad/s.

## Causal losses

A resistance-only intermediate render stayed finite but rang after excitation;
the full alphabet also clipped briefly. The current path includes physical
wall resistance, viscous series impedance, thermal shunt admittance, and a
passive one-pole radiation load.

Viscous/thermal channels approximate `sqrt(2 s)` by positive terms
`sum weight_j * s/(s+pole_j)`. A 24-point logarithmic quadrature spans poles
from `2 pi 1e-4` to `2 pi 1e10` per second. Its continuous transfer differs
from `sqrt(2 s)` by at most 0.200% on the measured 50 Hz–20 kHz grid.
Multiplying by the network's physical coefficient gives a real component
approximating `coefficient * sqrt(omega)` and retains the associated causal
reactive component. This is not a fitted human-speaker loss model.

For a branch, each channel stores canonical memory momentum `psi`, with
`psi' = weight * flow - pole * psi`. Its drop subtracts from the main branch
momentum derivative. Shunt channels use charge and pressure analogously.
The memory energy is `pole * psi^2 / (2 weight)` for positive weight. A joint
implicit midpoint update of the primary state and its memory branches is
passive at fixed geometry; tests check decay of their combined energy.

For radiation, let `L` be the existing end-correction inertance and `K` the
network's low-frequency radiation-resistance coefficient. The positive load
is `Z(s) = (L^2/K) s/(s+L/K)`. Its low-frequency expansion supplies reactive
inertance `L` and real resistance `K omega^2`. The end-correction inertance
is removed from the explicit acoustic branch, so it is not counted twice.
This one-pole continuation is an approximation, not a full radiation solution.

Loss and acoustic evolution use symmetric splitting with eight substeps per
48 kHz output sample by default. The source pressure integral is accumulated
through these substeps. `--canonical-substeps` changes temporal resolution
for convergence investigation, not output gain.

## Motion, output, and current limits

Positive compliance and nonnegative inverse inertance interpolate continuously
between control poses. This eliminates discrete state reprojection. Geometry
slopes and turbulence source placement are still updated at control boundaries;
a complete anatomical closure law, smooth source motion and a moving-wall
work audit including the loss memory remain unfinished. A vanishing loss
weight is not yet a verified physical closure limit.

Radiation uses the change in branch volume flow over the audio sample, so it
includes changing geometry and pressure-source contributions. The original
body geometry is retained; this integration does not add the detailed surface
model's transverse propagation or tongue-floor anatomy. Audio-rate held
sources and radiation averaging still need bandwidth/convergence assessment.
The result is a working research renderer, not a claim of correctly pronounced
or human-realistic NATO words.

Both CLIs explicitly print which acoustic path is running. The WAV CLI now
returns failure if its pre-clamp audit finds nonfinite or clipped samples,
while retaining the diagnostic file. The stream API accepts a synthesis
configuration; its canonical path passes no-allocation, host-buffer identity,
nonzero/finite/unclipped output and live additional-nudge checks.

```
./build/ducklips_say --canonical-motion --body bodies/reference.body \
  --dictionary samples/nato.dict --text 'Alfa.' --out out/canonical-alfa.wav
./build/ducklips_stream --canonical-motion --body bodies/reference.body \
  < examples/diphthongs.nudges > out/canonical-diphthongs.f32
```

The enhanced Alfa image preserves all 684 rows and 123 columns, the original
magma colors, the existing synthetic amplitude reference and 70 dB range.
It contains no ordinary-STFT companion, crop, grayscale conversion or reduced
resolution preview. Numerical speech receipts are recorded separately.

The first causal-loss Alfa render has no nonfinite or clipped samples.
Changing from eight to sixteen substeps gives 0.997% normalized waveform
RMS difference over the complete 1.31-second file, with the same 45 counted
fold cycles. This is a sensitivity measurement, not a converged-human-voice
claim. The late tail from 1.1 seconds has RMS 4.997e-5, versus 1.454e-3
in the resistance-only render (29.3 dB lower). These numbers concern the
new canonical path in both cases.

The full 26-codeword causal-loss run completed with zero nonfinite and zero
clipped samples, peak 0.666 and RMS 0.0122 without normalization. It used zero
modal eigensolves. Its 23.24 seconds of audio took about 309 seconds to render
while other validation jobs also ran; this is an offline integration receipt,
not an accepted performance result. Release and ASan/UBSan each passed all
seven registered tests, including the new streaming and energy checks.

The source/audio hashes, complete render logs, convergence comparison and
known limitations are in [canonical_speech_reference.json](canonical_speech_reference.json).
The final Alfa build produces exactly the same PCM as the causal-loss render
used to generate the enhanced field. No human word-recognition score has
been assigned to any of these files.
