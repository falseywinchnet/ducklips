# Connected sheet tube

The engine's material is a two-dimensional periodic sheet. Length runs along
one coordinate and circumference along the other; the circumferential edges
join. The two lines in the main editor are a side cutaway of this tube, not
two independently simulated strips.

## Geometry and controls

The below-mouth preset has nine sections and 60 axial positions, with eight
angular samples per ring: **480 radial degrees of freedom and 480 computed
modes**. The three shorter added sections have four axial samples each; the
six inherited sections retain eight each. Elastic connections cross every
section boundary and the periodic seam. Section materials differ in density,
thickness, modulus, damping and tissue-transport FIR coefficients.

The sequence is inlet expansion, light adapter, firm tissue, heavy soft tissue,
cord, heavy soft tissue, firm tissue, light adapter, and throat expansion.
The added upstream sequence reflects the downstream material order at half
its lengths. Main-tube diameter is about 5.373 mm. Both end expansions have
10.746 mm diameter: **2x diameter, 4x area**, not an inserted volume block.
Both end sections have unit mass, elasticity and support. Adapters inherit the
user's mass 0.2 and support 0.2, with elasticity 3; the cord has elasticity 1.5.
Heavy tissue retains mass 3.5, elasticity 1.5 and support 0.2. Firm tissue retains
mass 4.7, elasticity 0.2 and support 3.2; its base material modulus is higher.

The user's exact six-section configuration is preserved separately in
`presets/user-vocal-tube-2026-10-10.json` and selectable in the editor. It is
not silently migrated into the new geometry. The extended preset is exported
as `presets/below-mouth.json`.

There is one handle per section. It controls a ring support: vertical motion
changes the section's resting radius, and horizontal motion moves its support
footprint along that section. Both visible cutaway walls follow the same
ring. Mass, elasticity and support stiffness in the sidebar belong to the
selected section. Supports are circumferential rings, not isolated points.

The volume insertion tool is removed. Legacy JSON volume blocks still load
with their enclosed volume preserved. New designs use wall geometry alone.

Version 3 adds variable section sampling, whole-body elasticity (1–4), bounded
elasticity variation and a noise seed. Versions 1 and 2 still load. Version 1
migrates its two wall supports by averaging the two offsets and positions.
Input files are not modified. New edits use `ducklips-below-mouth-v3` in browser
storage; previous keys remain untouched.

## Supported sheet modes

`sheet.mjs` assembles a lumped mass M and symmetric stiffness K from axial
membrane tension, circumferential tension, hoop restoring stiffness, an
elastic foundation representing surrounding material, and the ring supports.
Axial distances include the radial slope between samples. The foundation's
15 mm support depth and the prestress/material parameters are explicit
prototype choices, not measured human anatomy.

Because rest geometry, materials and supports are uniform around each ring,
the periodic coordinate separates into an exact discrete Fourier basis.
Each angular family has a 60 × 60 mass-normalized eigenproblem. The full
set of real angular functions consists of a constant, three sine/cosine
pairs, and the alternating angular mode. All 480 eigenpairs are retained and
sorted by frequency; none are discarded to meet the runtime target.

The default dry supported-sheet frequencies span about **52.3–1062.4 Hz**.
There are **60 axisymmetric modes** and 420 circumferential modes. These are
material modes before air and delayed tissue-port loading, not acoustic
formants. The UI lists nominal modes at the set elasticity; live material
variation changes the supported response around that state.

The test suite applies the complete periodic-sheet K and M to every returned
mode, including the wraparound seam, and checks K phi = omega² M phi.
Changing whole-body mass by 4 divides frequencies by 2; changing elasticity
by 4 multiplies them by 2. Moving a support or reshaping volume recomputes
the supported response. Physical displacement and velocity survive edits;
editing the rest configuration can do work, for which there is not yet a
complete energy ledger.

## Efficient sum of modes

The runtime caches the sum of modes as a spatial receptance rather than
performing dense modal projections at each sample. For each angular family,
with mass-normalized coordinates, step h and eigenvector matrix Phi:

    R = Phi diag(1 / (1 + h² omega² / 4)) Phiᵀ
    z = h * velocity
    response = R (2 * displacement + z)
    next_displacement = response - displacement
    next_z = 2 * response - 4 * displacement - z

This is the implicit-midpoint response of the complete elastic modal system,
compiled when the design changes. At audio time it is a matrix application,
not a global solve. Its rapidly decaying spatial tails allow sparse storage:
terms smaller than 1e-14 are omitted, with the maximum omitted absolute row
sum measured and required to remain below 1e-12. Tests independently compare the cached update against explicit sums of all
modal responses, including the original 384-mode regression fixture.

With ring supports, a uniform inlet and a resting tube, all eight angular
sectors remain identical by symmetry. The normal playback path computes one
representative sector and its 60 axisymmetric modes. This removes duplicate
work; it does not replace the sheet with two walls or force an asymmetric
state to be symmetric. The full eight-sector path handles nonuniform initial
states, including diagnostic displacement of circumferential modes. Tests
compare its pressure-driven output with the compressed path through contact
and separately verify circumferential motion and air coupling. The diagnostic
`displaceMode` is not used as a running voice source or connected to the mode
preview control.

## Continual elasticity variation

Each section has its own seeded slow state: smoothed random targets, a
smoothed sawtooth relaxation, and signed irregular bubble envelopes. Their
bounded mixture changes elasticity by at most the selected amount (0–10%).
This is a phenomenological material-control law, not a validated neural or
muscle model. It does not inject noise into the pressure waveform or PCM.

The control runs at 200 Hz. Within each angular family the physical stiffness
matrix is tridiagonal. An LDL recurrence updates the full inverse response
`(I + h² K / 4)^-1` at control rate, independently checked against recomputing
the complete eigensystem at ±10% section elasticity. Audio applies the cached
response with sample-by-sample coefficient interpolation. All sheet modes
remain represented; there is no per-sample matrix solve. Each update verifies
that omitted response tails have absolute row sum below 1e-12.

The same varying elasticity updates nodal viscous damping and tissue-port
impedance and delay. Tissue histories reserve the longest delay at -10%
elasticity. Those local port coefficients update at control rate; only the
sheet-response coefficients interpolate per sample. This remains a numerical
approximation to a continuously changing material, with energy supplied or
removed by the changing constitutive parameters. A full energy ledger is
still absent.

Whole-body elasticity multiplies each section's own setting. In the dry sheet,
frequency scales as its square root, so 4x gives 2x mode frequency. The coupled
oscillation may change regime; this control is not a prescribed pitch track.
The browser and tests exercise the full 1–4 range without a model-domain fault
for the default preset. That does not guarantee stability for every edited
geometry or arbitrary abrupt movement.

## Air and tissue paths

Air connections run axially, diagonally across neighboring angular sectors,
and around each ring. Each directed connection has a geometric delay and the
same attenuation FIR. Fractional delay and attenuation share a combined FIR.
A delay has a one-sample minimum, so short-path dispersion remains a material
limitation at the 24 kHz physics rate. This is a representative pressure-wave
network, not a resolved fluid cross-section or a validated glottal jet model.

At each junction, incoming waves are combined by acoustic admittance. Delay
histories store pressure multiplied by sqrt(admittance), avoiding a simple
wave-energy jump when areas change. The same projected surface converts
pressure to radial force and wall velocity to air volume velocity. A scalar
midpoint exchange includes the local air impedance, tissue ports and viscous
damping. It does not perform a global per-sample solve.

Separate delayed connections in the material carry velocity wave amplitudes
axially and around the periodic seam. Their FIR taps are the geometric means
of the connected materials' taps. Their impedance derives from the local
membrane parameters with an explicit prototype coupling factor. These ports
provide additional delayed, lossy loading of the conservative sheet response;
their combined continuum/anatomical calibration is not established.

Both connection banks store delay histories contiguously. When all air delays
are exactly one sample, a common time-major history plane removes repeated
per-path circular addressing. Every path and direction remains independent;
fractional or longer delays use the original general representation. Edits
transfer the entire available history between representations. Tests require
sample-exact parity before and after such edits. This changed their
representation without changing the two-second reference PCM at all in the
measured comparison. Waveguide scattering follows the pressure/flow framework
in [Julius O. Smith's Digital Waveguide Theory](https://www.dsprelated.com/freebooks/pasp/Digital_Waveguide_Theory.html).
The supported-domain inspiration remains the user's
[drumskin mode example](https://www.youtube.com/watch?v=OEsBoatyHvU).

## Excitation, contact and limitations

The inlet reservoir is the only sustained external energy source. It approaches
the selected pressure through a 35 ms response with seeded slow pressure
wander. The expanded lower sections, central vocal cord and upper sections remain
bidirectionally connected. There is no pitch oscillator, prescribed phonation
wave, phoneme sequence, output timbre correction or voice recording in the
engine.

Contact currently imposes a 25 micrometre radial floor and removes inward
velocity at contact. This is an inelastic nonpenetration approximation, not
an airtight collision/jet model. The outlet is observed with a two-tap pressure
difference and fixed conversion gain. The monitor has a separate level and
clip guard; the recorded float WAV contains the raw output.

The sheet describes radial membrane motion on a cylindrical topology. It is
not a full shell model with tangential motion, a multilayer vocal fold, or a
3D fluid solver. Material values, external foundation, contact, port loading,
short-path delay floor and flow laws need physical calibration. Implementing
the requested topology and modes does not establish natural human phonation.
The new below-mouth field contains sustained, moving bands. It has not been
accepted as a realistic voice. Mouth, sinus, under-tongue, teeth and lip systems
are future connected stages, not hidden components of this preset.

## Running, timing and paired depiction

    node web/tube/server.mjs
    node --test web/tube/*.test.mjs
    node web/tube/render.mjs out/sheet-tube-check 2

Physics advances at 24 kHz and is linearly reconstructed to 48 kHz playback.
The initial dense modal path measured about 4x realtime. Cached receptance,
rotation symmetry and contiguous connection histories improved the measured
normal playback path to about 11x for the expanded, continually varying preset on this host. This is a pressure-on core
measurement, not a universal promise for every device, design or fully
asymmetric state. Compilation, drawing, audio encoding and offline enhanced
analysis are outside that measurement. The benchmark runs the loaded design
at page startup; subsequent geometry changes do not rerun it automatically.

Each recording automatically submits its entire raw WAV to the local enhanced
renderer. Its complete magma field preserves every computed row and column.
The viewer scrolls rather than shrinking it; no ordinary STFT companion or
cropped substitute is added. The receipt freezes the engine, sheet compiler,
audio renderer, both browser transports, app source and recording provenance.

The server binds 0.0.0.0:8765 for LAN access and prints interface addresses.
Use HOST=127.0.0.1 to restrict it to this computer. No build step or npm packages
are required. Localhost/HTTPS uses AudioWorklet; HTTP LAN uses the same renderer
in a Worker feeding scheduled Web Audio buffers, adding about 0.2 seconds of
playback buffering. Block-size parity tests cover the two paths. The host must
remain awake with the server running. The enhanced renderer uses
~/miniforge3/bin/python unless DUCKLIPS_PYTHON is set.

## Connected speech tubes (version 4)

The accepted nine-section below-mouth configuration remains intact as the source
stage and remains separately selectable. The default **Open** arrangement connects
it to a four-section mouth and a three-section teeth/lips tube. Two shorter tubes
branch from mouth sections 1 and 3: a sinus pocket and an under-tongue pocket.
Both pockets terminate closed in this experiment. The only driven reservoir is
the below-mouth inlet; the sampled open end is at the lips.

Each stage is its own periodic 2D material sheet. Its complete nominal modes are
computed on support/material edits. The network uses these stage responses and
couples their pressure and material velocities in both directions; it does not
cascade rendered audio through filters. A side connection is distributed around
its section's ring in this first model, preserving the rotational symmetry used
by the existing engine. This is an experimental tract topology, not anatomical
sinus/nasal airflow or a calibrated oral cavity.

Every inter-tube link has separate two-way air and solid histories. Its acoustic
admittance uses the smaller instantaneous opening. Air propagation uses the
existing three-tap attenuation FIR and a causal fractional delay based on the
adjacent half-cell lengths. Tissue propagation uses the endpoint material FIRs,
masses and elasticity. The link enters each endpoint's existing collocated
pressure/velocity junction. Incoming waves for **all** tubes are read before any
tube advances, and outgoing waves are committed afterward. Thus child evaluation
order cannot create a one-way connection or zero-delay loop. Normalized link
histories preserve wave-power scaling as the junction's aperture/impedance moves.
This does not prove global passivity under time-varying geometry and elasticity;
the existing active material changes can exchange energy with the system.

All stages receive the single 1–4 system elasticity multiplier and independent
seeded, bounded material variation. Tube-level source controls inside a version 4
file are overridden by the system controls: downstream tubes never gain a second
reservoir. Local section elasticity is retained. Saving exports all tubes,
connections, and controls. Older single-tube JSON still loads. The new editor uses
a separate local-storage key and does not overwrite the accepted below-mouth save.

Three inspectable starting arrangements:

| Arrangement | Mouth | Teeth/lips | Total computed modes |
| --- | --- | --- | ---: |
| Open | Four sections, wider tongue/body openings | Soft vestibule, firm teeth gap, elastic lips | 832 |
| Balanced | Four sections, narrower front opening | Firm teeth gap, elastic lips | 800 |
| Rounded | Three sections; long, heavy, weakly supported soft front chamber | Firm teeth gap, longer/narrower elastic lips | 768 |

The names describe geometry, not validated vowels. A mode preview is the selected
stage's nominal sheet mode, not an eigenmode of the entire coupled air/tissue
network. At the current ring-symmetric drive, the angularly asymmetric modes have
zero excitation; they remain computed and available in the full-sheet diagnostic.

### Compiled sample kernel

The connected JS engine crossed below the 10× core target. `native/core.c` carries
out the same sample operations with double precision in contiguous WebAssembly
memory, using the same histories, receptance coefficients, contact rule and
junctions. Control-rate material variation and geometry compilation remain in JS.
There is no new source waveform, reduced mode count, downsampled physics rate, or
output filter. Native/reference tests include exact PCM equality through changing
elasticity, closure/contact, fractional transport, support edits and network
feedback. The full asymmetric diagnostic retains the JS path.

`native/build.sh` regenerates `core-wasm.mjs`; the generated module is committed,
so serving and rendering require no compiler. Its embedded byte array also runs
inside AudioWorkletGlobalScope, where `atob` is unavailable. The current build uses
Homebrew Clang and the existing Rust wasm linker; environment overrides are
`DUCKLIPS_CLANG` and `DUCKLIPS_WASM_LD`. No fast-math or fused floating-point
reassociation is enabled.
