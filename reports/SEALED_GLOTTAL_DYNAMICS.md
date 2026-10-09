# Retained glottal air during elastic contact

Historical phase-equation report at commit `c7e3f0e`. The later
[contact-transition report](CONTOUR_CONTACT_TRANSITIONS.md) adds automatic cap
sealing/reopening controls and identifies the remaining pressure-driven failure.
The fixed-topology limitations below describe this earlier validation.

`ContourGlottalStorage` can now evolve either sealed-cap phase with two
positive-volume chambers. It retains their air masses, the active passage
momenta and the tissue state. A sealed passage has exactly zero flow and
requires zero momentum. The source also includes an elastic contact energy
and the corresponding tissue and rest-actuator work.

The dynamic solver still requires the set of open passages to remain fixed
during each advance. Automatic sealing and reopening are unfinished. These
are verified phase equations, not a completed oscillation cycle or a new
speech recording. The latest listening feedback remains attached to the
spatial-pressure recording: November improved; Foxtrot is jammed around /tr/.

## A sealed cap is not an empty chamber

The existing signed medial profile is integrated over its positive part.
A negative control opening represents tissue compression; the neighboring
sloping surface can still enclose air. `ContourKinetics` now computes volume,
air centroid, inertia and their derivatives for these clipped contours.
The two chambers must both retain positive volume and mass. An empty chamber
continues to be rejected, since removing its state requires separate transport.

Clipping makes volume nonlinear in the two control areas. Gas force now uses
a symmetric finite volume gradient from the same coordinate paths as the
kinetic-energy transfer. Thus

```
Delta V_i = sum_k DV_i,k Delta area_k
Delta H_i = mu_bar_i Delta M_i - P_bar_i Delta V_i.
```

The implementation differentiates the clipped interval integrals with paired
endpoint arithmetic. For an affine segment, the positive-part square and
cube primitives give its volume and first moment. Constant end segments use
the positive part directly. Zero-area plateaus use their closed-side
derivative, matching the existing geometry convention. This convention is
explicit: a contact boundary need not have a unique two-sided derivative.

A closed cap has zero inverse inertia and zero represented flow. Its
resistance field is an inactive zero sentinel, not a finite physical
resistance used in an open-passage law. The open-limit resistance diverges.
A nonzero stored momentum on a closed cap is rejected. This prevents a
hidden state with zero present energy from producing motion after reopening.
The blocked port supports an effort difference with zero power transfer.

During a sealed phase the solver imposes the zero-flow constraint, advances
only the active passage momenta, and solves the tissue and retained gas
states jointly. It checks the topology of every proposed endpoint. There
is no added area floor or automatic momentum reset.

## Contact from an elastic energy

The new `FoldContact` component uses the existing body stiffnesses, tension,
asymmetry, width and `collision_ratio`. For each lower/upper pair, let `kL`
and `kR` be its left/right tissue stiffnesses after tension scaling. Define

```
k_pair = collision_ratio / (1/kL + 1/kR)
gap    = signed_area / width
U      = (1/2) k_pair min(gap, 0)^2.
```

The series stiffness represents two compressed tissues carrying a common
contact force. Each side receives `-dU/dgap`, in its own outward coordinate.
For a symmetric body this reproduces the existing linear collision-force
magnitude. For an asymmetric body it supplies one shared reaction derived
from one energy. It is a declared lumped elastic assumption, not measured
contact stiffness or an anatomical fit. Base tissue damping is retained;
this component does not add a conditional collision damper.

The finite contact force is the exact secant of this potential, including
steps across zero in the constitutive test. The rest-area actuator also does
its share of contact work. Contact energy is included in source energy and
in the same reservoir/tract/external-force/dissipation ledger.

Lumped compression represented by overlap and collision springs also appears
in [Mora et al., Energy-based fluid–structure model of the vocal folds](https://pmc.ncbi.nlm.nih.gov/articles/PMC8210679/),
section 3.1 and Appendix A. Their contact potential is continuous through
collision. The current paired linear law is our explicit simpler assumption.
Their section 3.3 introduces a threshold that disables fluid dynamics near
closure; that treatment has not been adopted here. This reference provides
context for energy-based contact, not validation of our parameter values.

## What vanishes at sealing

For either constant-area cap segment, let its axial length be `c`, its small
positive opening `a`, and its chamber relative density `r=M/V`. When chamber
volume and density remain finite and positive:

```
L ~ rho0 c / (r a)
R ~ 12 viscosity width^2 c / (r^2 a^3)
tau = L/R ~ rho0 r a^2 / (12 viscosity width^2).
```

With positive viscosity and a bounded mass-effort drop, the corresponding
quasisteady cap flow, canonical momentum and energy obey

```
Q ~ a^3
lambda = L Q ~ a^2
K = (1/2) L Q^2 ~ a^5.
```

These are limiting constitutive relations, not a completed event simulation.
[The native-metric limit audit](contour_sealing_limit.json) checks both caps
from `1e-8` through `1e-13 m²`, with fixed positive retained masses and a
700 Pa mass-effort drop. Fitted small-area exponents agree with `(3,2,5)`
within `1.5e-5`. At exact contact the native geometry retains positive chamber
volume and returns an inactive cap. The audit does not assume this bounded
pressure limit holds when a whole chamber disappears, nor with zero viscosity.

The limit explains what a sealing transition must approach: zero cap flow
and energy, while retaining adjacent gas storage and other active momenta.
The rapidly shrinking relaxation time also explains why blindly extending
the open implicit step towards contact is difficult. It does not justify
resetting an arbitrary finite momentum or discarding a chamber at a cutoff.

## Verification

[The signed-contour reference](contour_sealed_reference.json) records 160
instantaneous checks and 92 finite transfers against independent 80-digit
calculations. The latter include 12 open-to-contact endpoint pairs. Maximum
scaled discrepancies are `1.08e-14` for instantaneous forces, `2.37e-13` for
finite kinetic forces, and `5.29e-14` for finite volume gradients. These are
constitutive checks; the dynamic solver still rejects an unowned topology
change.

Four fresh independent continuous controls evolve an asymmetric body with
one cap sealed: lower contact, upper contact, and each with a changing rest
opening. Unequal tissue stiffnesses, external closing forces, gas storage,
active passage inertia and contact work all participate. DOP853 with maximum
steps 0.2 and 0.1 microsecond gives normalized endpoint differences below
`1.02e-11`. Native trajectories preserve exact sealed-port zero flow and
zero momentum, with accumulated mass errors below `1e-19 m³` and work errors
below `1e-13 J`.

The initial seven native grids (100 through 6,400 steps over 0.2 ms) meet
endpoint accuracy for all four cases. The lower-contact cases still fail one
of their final convergence-ratio checks; that initial failure is retained.
[One further native refinement](contour_sealed_endpoint_reference.json),
using the same verified continuous endpoints, gives ratios **4.0453** and
**4.0549** at 12,800 steps. Together with the preceding accepted ratios,
all four cases pass the original second-order and endpoint criteria.
This extra check reuses its continuous oracle rather than integrating it again.

[The open-state regression](contour_open_regression.json) checks the changed
source against the independent open equations preserved at commit
`5bc0d1ace9273b158346921fd1a652edbb8c508f`. The original oracle report and its
equations/gas/body definitions are verified against committed bytes. Fresh
native forward, reverse, relaxation and rest-motion refinements all pass.

Twelve focused CTests pass. Contact-potential work, time reversal, absence
of contact force across an open gap, stale-momentum rejection and existing
state-restoration/no-allocation checks pass. Address/undefined-behavior
sanitizers pass the open, sealed and asymmetric controls. Sanitized signed
kernel and finite-transfer fixture outputs are byte-identical to native.
[Validation receipts](contour_sealed_validation.json) record the logs and hashes.

## Remaining transition work

The event owner still needs to join these phase equations. It must locate a
sealing event while integrating the fast viscous relaxation and its work,
arrive at an inactive cap with zero energy, and preserve both retained gas
states. Reopening must reintroduce the passage from that zero-energy limit
while retaining the rest of the source state. If a chamber loses all volume,
its air and energy need an explicit transport rule before the state can be
removed. Those operations are not implemented by accepting a sign change
in this fixed-topology step.

The current reference remains a reduced two-chamber geometry and an explicit
linear contact law. Physiological calibration, a full glottal cycle, connection
to the speech renderer, new pronunciation evidence and production efficiency
remain outstanding. Performance optimization is still parked.

To reproduce with the project's scientific Python environment:

```
cmake -S . -B build
cmake --build build --target ducklips_contour_sealed_probe ducklips_contour_glottal_probe
python tools/contour_sealed_reference.py
# The initial-grid ratio failure remains in its report; refine it next.
python tools/contour_sealed_endpoint_reference.py
python tools/contour_open_regression.py
python tools/contour_sealing_limit.py
build/ducklips_contour_sealed_probe
```
