# Flat tube experiment

This is a fresh implementation of the requested side-cutaway tube. It does not
import the previous Ducklips synthesizer, source, gesture planner or recordings.
The first objective is an editable, inspectable physical instrument. Natural
human phonation and intelligible speech are not established by this version.

## Direct manipulation

The two walls define an equivalent circular tube: A = pi ((upper + lower)/2)^2.
The cutaway retains separate upper/lower motion, but it is not a full 3D fluid
or shell discretization. Six bands contain four axial samples each. Band 1 is
the return chamber; band 2 is the initially tight elastic band. There is no
special waveform source in that band. All bands use the same engine.

Volume blocks specify added enclosed air volume, in mL. Their compact spatial
weights are normalized on the actual axial grid, so integral(delta A dx)
equals their stated volume. They expand the walls; they are not air obstacles
or additional inertial tissue. Dragging a block relocates that volume.

Each surface has a support handle. Vertical movement changes its rest position;
horizontal movement relocates the support spring within that material patch.
The support slider changes stiffness. At present each segment's vertical rest
offset is uniform; this is a coarse patch model, not an arbitrary spline wall.
Undo, reset and JSON save/open retain the physical design. The browser saves
the current design locally. Nothing is published or uploaded to an external site.

## Receptive materials

Each surface patch has a diagonal mass matrix and a symmetric stiffness matrix
assembled from membrane tension, a bulk support term, and the positioned
support. Wall arc length affects prestretch. A small eigendecomposition of
M^(-1/2) K M^(-1/2) produces four modes. The displayed frequencies are results
of that calculation, not pitch targets. Material density, modulus, thickness,
damping and three FIR coefficients are explicit in engine.mjs. These are
prototype material settings, not measurements of a particular human body.

Modal displacement/velocity retains physical state when a support changes:
old nodal displacement and velocity are projected into the new basis. Moving
a support or inserting volume changes rest geometry and can do external work;
a complete ledger for that editing work has not been implemented.

The linked example [Computing vibrations of a moose drumskin, Alexander
Gustafsson](https://www.youtube.com/watch?v=OEsBoatyHvU) demonstrates modes tied
to the supported domain. This implementation uses that principle on a much
coarser membrane patch. It does not claim to reproduce the video's arbitrary
2D/3D domain calculation.

## The two communication paths

Air has two representative lanes with axial, transverse and crisscross links.
Each link has geometric propagation delay and the same three-tap attenuation
FIR. Fractional-delay interpolation is folded into the attenuation taps, so
there is one combined FIR per directed link and no iterative propagation pass.
At a junction, admittance-weighted incoming waves determine pressure and the
outgoing waves. These relations follow pressure continuity and volume balance;
see [Julius O. Smith, Digital Waveguide Theory](https://dsprelated.com/freebooks/pasp/Digital_Waveguide_Theory.html).

Internally, delay buffers store sqrt(admittance) times pressure. This energy
normalization avoids inventing wave energy simply by changing a port area.
The junction receives pressure, and wall motion exchanges volume with it.
The same surface S converts p to force and wall velocity to volume velocity.
A closed-form local midpoint exchange accounts for the air-port impedance;
it requires no global per-sample linear solve. Supported material responses
advance as modal states; air and solid connections each advance once per step.

Solid connections carry traveling velocity amplitudes between matching nodes
of adjacent material patches on each wall. Their FIR is the geometric mean
of the adjacent materials' attenuation taps. It transmits only what survives
both materials, including their different losses. Finite delays use a speed
estimate from the lowest material mode and patch size. Modal receptance is
local to each material, while the FIR governs propagation between materials;
these are distinct operations. The mechanical port impedance is currently a
fixed prototype value. This path is a representative network, not a validated
continuum elastic propagator. Tissue modes themselves use recursive states,
not a truncated sampled impulse-response convolution.

## Pressure, contact and observation

Only the inlet reservoir supplies sustained external energy. Its target
pressure approaches through a 35 ms reservoir response. Seeded low-rate
stochastic pressure wander provides the requested slight variability. It is
not an audio-rate noise source or a prescribed oscillation. The first chamber
remains bidirectionally connected to the second band and the rest of the tube.

Opposing walls have a nonadhesive compression/recoil contact force. Minimum
hydraulic radius is 25 micrometres: complete airtight sealing is not yet
represented. The present air model is a moving acoustic network with reservoir
flow; it does not yet resolve a separated glottal jet, shear-layer turbulence,
or a validated viscous constriction law. Consequently, any observed sustained
oscillation is a candidate behavior of this approximation, not established
physiological phonation. In particular, regular spectral lines do not prove
that its energy transfer is anatomically correct.

The outlet sample is a two-tap pressure difference, a simple radiation
observation, with a fixed conversion gain. There is no source sine, imposed
pitch, phoneme scheduling, sample playback, output EQ or timbre correction.
Live monitor level and a final monitor clip guard are separate from the raw
float WAV capture. Faults stop the model visibly rather than substituting a
fallback sound. Stopping pressure releases the reservoir and lets the same material and
air states ring down. Reset clears state.

## Time and representation

The engine advances at 24 kHz, with causal linear reconstruction to a requested
48 kHz Web Audio context. Short propagation paths have a one-engine-sample
minimum delay; mesh dispersion and this delay floor require calibration before
physical formant claims. Four modes per surface are retained. This is a
bandlimited representative model, not 48 kHz fluid mechanics.

The first browser measurement was below 10x. Combining fractional delay and
attenuation into one FIR and storing only the required path history improved
the representation without reducing the modal count. The recorded final
0.8-second render measured 10.43x including trace collection; a subsequent
browser core measurement should be reported separately. Compilation, UI,
record encoding and offline enhanced analysis are outside the core benchmark.
No universal speed guarantee is implied.

## Audio -> depiction

Every recorded take automatically submits its complete raw WAV to the local
server for the existing enhanced field renderer. It returns the full magma
field at native computed cell size, with metadata and NPZ data. It does not
add an ordinary-STFT panel, crop the take or downsample the image. The viewer
scrolls instead of shrinking. The waveform uses min/max bins so narrow spikes
remain visible. Enhanced analysis is offline and takes longer than synthesis.

Run from the repository root:

    node web/tube/server.mjs
    node --test web/tube/engine.test.mjs
    node web/tube/render.mjs out/flat-tube-new 0.8

The server binds only 127.0.0.1:8765. No npm packages or build step are required.
The enhanced renderer uses ~/miniforge3/bin/python unless DUCKLIPS_PYTHON is set.
