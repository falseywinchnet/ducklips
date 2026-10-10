# Pressure-driven deforming tissue source

The user rejected the isolated two-mass source as nonhuman, then explicitly
directed development toward deformation and recoil governed by pressure and
elasticity. This experiment is source work, not another pronunciation or speed
tuning round. No naturalness result is accepted.

## What changed in the representation

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
