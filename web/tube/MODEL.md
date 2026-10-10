# Connected sheet tube

The engine's material is a two-dimensional periodic sheet. Length runs along
one coordinate and circumference along the other; the circumferential edges
join. The two lines in the main editor are a side cutaway of this tube, not
two independently simulated strips.

## Geometry and controls

Six sections each contain eight axial positions. Eight samples around each
ring make a 48 × 8 sheet with **384 radial degrees of freedom and 384 computed
modes**. Elastic connections cross every section boundary and the periodic
seam. Section materials may differ in density, thickness, modulus, damping
and tissue-transport FIR coefficients.

There is one handle per section. It controls a ring support: vertical motion
changes the section's resting radius, and horizontal motion moves its support
footprint along that section. Both visible cutaway walls follow the same
ring. Mass, elasticity and support stiffness in the sidebar belong to the
selected section. Supports are circumferential rings, not isolated points.

Volume blocks add enclosed air volume, not obstacles. Weights on the axial
grid are normalized so sum(delta A * dx) equals the requested mL. For the
resting circular geometry A = pi r². In the full moving state, angular
quadrature averages pi r(theta)². The editor retains its side cutaway while
the optional mode view displays the complete unrolled material sheet.

Version 2 designs have one `offset` and `anchor` per section. Version 1 saves
are migrated by averaging their upper/lower offsets and anchor positions;
the original saved JSON is not modified. The old browser storage key remains
available. New edits save under a separate version 2 key.

## Supported sheet modes

`sheet.mjs` assembles a lumped mass M and symmetric stiffness K from axial
membrane tension, circumferential tension, hoop restoring stiffness, an
elastic foundation representing surrounding material, and the ring supports.
Axial distances include the radial slope between samples. The foundation's
15 mm support depth and the prestress/material parameters are explicit
prototype choices, not measured human anatomy.

Because rest geometry, materials and supports are uniform around each ring,
the periodic coordinate separates into an exact discrete Fourier basis.
Each angular family has a 48 × 48 mass-normalized eigenproblem. The full
set of real angular functions consists of a constant, three sine/cosine
pairs, and the alternating angular mode. All 384 eigenpairs are retained and
sorted by frequency; none are discarded to meet the runtime target.

The default dry supported-sheet frequencies span about **72.2–1234.2 Hz**.
There are **48 axisymmetric modes** as well as 336 circumferential modes.
These are material modes before air and delayed tissue-port loading, not
acoustic formants. Their frequency range follows the chosen geometry and
constitutive parameters. A finite sheet has a finite spatial bandwidth.

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
sum measured and required to remain below 1e-12. The default axisymmetric
kernel stores 262 of 2304 coefficients. Tests independently compare the
cached update against explicit sums of all 384 modal responses.

With ring supports, a uniform inlet and a resting tube, all eight angular
sectors remain identical by symmetry. The normal playback path computes one
representative sector and its 48 axisymmetric modes. This removes duplicate
work; it does not replace the sheet with two walls or force an asymmetric
state to be symmetric. The full eight-sector path handles nonuniform initial
states, including diagnostic displacement of circumferential modes. Tests
compare its pressure-driven output with the compressed path through contact
and separately verify circumferential motion and air coupling. The diagnostic
`displaceMode` is not used as a running voice source or connected to the mode
preview control.

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

Both connection banks store delay histories contiguously. This changed their
representation without changing the two-second reference PCM at all in the
measured comparison. Waveguide scattering follows the pressure/flow framework
in [Julius O. Smith's Digital Waveguide Theory](https://www.dsprelated.com/freebooks/pasp/Digital_Waveguide_Theory.html).
The supported-domain inspiration remains the user's
[drumskin mode example](https://www.youtube.com/watch?v=OEsBoatyHvU).

## Excitation, contact and limitations

The inlet reservoir is the only sustained external energy source. It approaches
the selected pressure through a 35 ms response with seeded slow pressure
wander. The return chamber, second vocal section and rest of the tube remain
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
The current fresh field contains broadband energy and moving bands; it has
not been accepted as a realistic voice.

## Running, timing and paired depiction

    node web/tube/server.mjs
    node --test web/tube/*.test.mjs
    node web/tube/render.mjs out/sheet-tube-check 2

Physics advances at 24 kHz and is linearly reconstructed to 48 kHz playback.
The initial dense modal path measured about 4x realtime. Cached receptance,
rotation symmetry and contiguous connection histories improved the measured
normal playback path to about 11–13x on this host. This is a pressure-on core
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
