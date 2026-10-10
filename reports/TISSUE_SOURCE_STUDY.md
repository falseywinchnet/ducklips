# Pressure-driven deforming tissue source

The user rejected the isolated two-mass source as nonhuman, then explicitly
directed development toward deformation and recoil governed by pressure and
elasticity. This experiment is source work, not another pronunciation or speed
tuning round. No naturalness result is accepted.

## Current elastic volume representation (round 006 onward)

The scalar support experiment below has been superseded in `TissueSource`.
`tools/compile_tissue_body.py` constructs a two-layer elastic volume from
`data/tissue_body.json`, with all three displacement components. The backing
and anterior/posterior attachments are fixed; the other faces are free.
Selective integration of bulk strain permits nearly incompressible material
without locking the coarse hexahedral elements. Deviatoric strain is fully
integrated. Longitudinal prestress contributes geometric stiffness.

The compiler projects 735 free volume coordinates onto 16 mass-normalized
elastic deformation modes. Projected viscosity remains a coupled matrix.
The generated body contains the stiffness, damping, surface observation and
pressure-force projection operators. The runtime evolves only the 16
coordinates, driven by actual geometry-dependent pressure and contact. Mode
frequencies are body properties, not prescribed voice frequencies or phases.
The full tissue solve occurs only during body compilation (about 0.2 s here).

This removes the scalar relation previously assumed between normal surface
motion and interior shear. The internal displacement modes include vertical
motion, bulging and rotation. A static volume solution supplies the posturing
influence; changing its equilibrium coordinates preserves actual position and
velocity. The 16-mode surface posture approximation has an 8.4% relative
spatial error, which is a known truncation limit, not an acoustic quality score.

The current volume model is linear-strain elasticity about a reference body.
The aerodynamic approximation is still quasi-steady, with small-slope lateral
tractions and fixed axial sampling. It does not yet account for displaced air
volume or two-way airway acoustics. No claim of complete physiological accuracy
or accepted voice quality follows from the body-operator checks.

The tightly adducted round 006 did not sustain phonation. Its source-only
computation took 0.0695 s per generated second (about 14.4x), excluding trace
writes and posture initialization. This is a source observation, not a
full-engine performance qualification.

Round 007 uses the same elastic body with the superior edges just apposed at
rest, rather than carrying the scalar model's -0.6 mm compression into a
different body. Under constant 600 Pa breath, oscillation grows from rest,
settles into recurring pulses and includes upper-edge contact and zero-flow
intervals. The first 0.2 s is the onset transient; release starts at 0.8 s.
The full enhanced field has sustained harmonic structure without the broad
intermittent streaks of the rejected early surface runs. This is a sustained
source result, not an accepted natural-human timbre result.

Its measured source computation is 0.07084 s for one second of output
(14.12x), excluding diagnostic writes and initialization. Playback is one
constant gain, shared with the prior isolated-source diagnostic; there is no
waveform shaping, EQ, added noise or imposed F0. The latest raw/source WAV,
enhanced field, complete physical traces, snapshots and executable are at
`out/tissue-round-007`. Inspection and provenance are summarized in
`reports/tissue_source_reference.json`.

The focused source checks now distinguish sustained phonation from static
airflow: the initially apposed body must develop recurrent open/closed flow
under constant pressure. Zero-pressure rest, unforced energy decay and
adduction continuity are also checked. All pass. These checks cannot certify
perceptual realism. The experiment remains separate from the default speech
engine until source behavior and coupling are established.

## Initial surface-only representation (rounds 001–005)

`TissueSource` represents lateral displacement over two elastic layers. Each
layer has a 13 by 13 surface lattice, with fixed anterior/posterior attachment
points and free inferior/superior edges. The cover moves relative to a deeper
compliant layer, whose outer boundary is attached. Length, thickness, density,
prestress, shear modulus and viscosity enter the displacement equations in SI
units. Cover/backing coupling includes a positive quartic stored-energy term
for strain stiffening. Contact pressure grows continuously with local
penetration; its dissipative term vanishes at first contact.

For a cover displacement u and backing displacement b, the elastic surface
forces are the divergence of longitudinal and vertical stresses. The cover
inertia per unit area is rho*d_cover. A linear displacement profile through
the attached backing gives rho*d_body/3. Local pressure, layer shear and
contact supply the remaining forces. No oscillation frequency, temporal sine,
pulse shape, cycle phase, jitter, or random kick is supplied. A fixed 600 Pa
breath develops over 12 ms; its release begins at 0.8 s. The source runs for
one second, including onset and recoil.

This is a deliberately limited continuum approximation: lateral displacement
only, symmetric folds, simplified material law and boundary conditions. It is
not a full incompressible three-dimensional tissue model or a reproduction of
any published engine. The constants in the first round are provisional body
parameters, not measurements of the user.

## Airflow approximation and current limits

The first round uses quasi-steady viscous slit flow. Open paths across the
fold length contribute parallel conductance; cross sections along the air path
contribute series resistance. The kinetic exit loss and viscous resistance
determine the flow from the pressure difference. Surface pressure follows that
same pressure budget. Separation at 1.2 times the minimum area is a closure
assumption taken from reduced airflow research, not a discovered universal law.
Air inertia, volume displaced by the moving walls, reverse flow and turbulent
fluctuations are not yet represented. These limitations prohibit claiming an
exact physiological source or noise prediction.

Round 001 develops oscillation without a periodic driver, but it has irregular
cycles, overlarge excursions and abrupt steps in airflow. Its surface sequence
shows deformation and contact moving across the tissue. The airflow steps
coincide with a discretization defect: separation is snapped to an entire
cross section, changing the chosen exit area and integrated resistance. The
first correction to investigate is continuous separation within the surface
geometry, before adjusting material constants.

The full one-second enhanced depiction uses all 684 frequency rows and the
128-sample time lattice, with the original magma colors. No ordinary STFT
image is used. Supplementary displacement and force figures explain the
mechanism; they do not replace the full sound depiction. Each round retains
raw and playback WAVs, source snapshots, trace, complete surface history,
hash receipt, full enhanced field and native-cell PNG.

## Prior-source diagnostic

`tools/study_source_forces.py` observes the actual substep forces in the frozen
source-stage engine. Its newly rendered raw and playback audio are byte-equal
to that source-stage baseline, proving that observation did not change the
signal. The trace exposes abrupt contact damping and surface exposure, and
near-zero aerodynamic loading on the upper mass for the unloaded boundary.
The complete enhanced field for this signal was inspected. This identifies
assumptions; it does not prove which one caused the reported tinny percept.

## Research used to check the physical assumptions

- [Farahani and Zhang, 2016, experimental validation of reduced continuum
  phonation](https://pmc.ncbi.nlm.nih.gov/articles/PMC5384606/): geometry-based
  flow separation and the limitations of a reduced aerodynamic model.
- [Serry, Zanartu and Peterson, 2026 preprint](https://arxiv.org/html/2606.13480v1):
  a distributed membrane can retain tissue geometry at modest computational
  cost. Our experiment does not implement their beam/posturing model or use
  their source code.
- [Kazemirad and colleagues, 2014, in-vivo tissue shear
  measurements](https://pmc.ncbi.nlm.nih.gov/articles/PMC3961822/): material
  elasticity and viscosity require empirical grounding. The initial low
  viscosity has not been accepted as a human tissue configuration.
- [Titze, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9533394/): surface-wave
  motion and airway loading affect the pressure work that sustains vibration.

The first source-only timing was 0.108 s per second of generated signal,
excluding trace writes, while other analyses were running. This is about
9.3x real time and is not a qualified full-engine performance result.
Performance tuning remains parked as requested. If the distributed field is
retained, a smaller spatial basis derived from its elastic operator is the
representation to investigate; temporal waveforms must remain emergent.

## Surface-only rounds 002–005

Round 002 moved separation continuously between sections and integrated the
exposed fraction of each surface patch. The widespread flow steps disappeared,
although the full enhanced depiction retained an onset episode and intermittent
vertical streaks. Round 003 used the paired in-vivo cover estimate G=878 Pa,
eta=0.52 Pa s. It sustained a nearly sinusoidal open-gap cycle, not useful voice.
Round 004 initialized an adducted posture (central equilibrium gap -0.6 mm,
0.3 s of zero-pressure equilibration). It settled to static airflow. Round 005
allowed independent local cover rotation, separating rigid rotation from shear
strain. That altered the transient and equilibrium but still failed to sustain
voicing. All five complete enhanced fields were inspected; none is accepted.

The adduction continuity check caught a stale geometry accessor: positions had
advanced while the cached gap still represented the previous substep. Moving
the geometry refresh to the end of the step fixed the observation without
changing the generated raw audio (byte equality checked against round 004).
Checks for zero-pressure rest, unforced mechanical dissipation, and position
and momentum continuity under adduction now pass.

The scalar cover/backing support law remains an unvalidated reduction. In
particular, a measured shear material coefficient does not by itself validate
the assumed relation between normal surface motion and interior deformation.
The next representation will derive coupled lateral/vertical/longitudinal
motion from an elastic volume, retaining a small spatial deformation basis in
the running source. This supersedes further tuning of the scalar support law.
