# Viscous drag along the nonuniform lip profile

The developing continuous-source airway now has an explicit local slit-drag
model on the three open internal lip paths. The preceding closing-control
refinement discrepancy was concentrated at the last lip storage section.
Adding drag determined by the local aperture, while retaining the existing
integrated air inertia, reduces its finest-step pressure difference from
39.539 Pa to 0.03491 Pa. The largest closing difference anywhere is now
0.11336 Pa in the throat. This is numerical evidence for the modeled airway,
not a new intelligibility result or validation against human lip measurements.

## Geometry and physical scope

For each path between adjacent volume centroids, the existing axial profile
provides area A(x). A declared slit width W = 0.03 m gives half-gap
`a(x) = A(x)/(2 W)`. The path impedance is approximated by integrating local
stationary-wall slit impedances with shared volume flow. The inherited graph
already owns the plug-flow inertia `L = rho integral(dx/A)`; the added sites
supply transverse velocity-profile memory and drag only.

The new `lip_viscous_segments` helper places positive Gauss-Legendre weights
along the actual path. Each site has local half-gap and weighted inertia
`L_i = rho w_i/A_i`. The weights are normalized against the profile's existing
integrated inertia so that `sum(L_i) = L` exactly to rounding. This correction
is recorded, and is below 1.14e-7 relative at order eight in these controls.
It does not alter the graph's air-column inertia to obtain damping.

The existing 24-mode finite-slit kernel acts on the common branch flow. Its
positive relaxation terms retain profile memory, stored kinetic energy and
nonnegative viscous dissipation. Changing geometry retains each site's
canonical memory and accounts for its changed energy metric. Quadrature
identities follow fixed reference abscissae between the moving centroids.
The adapter applies sites in forward order before the graph/source step and
reverse order afterward; time refinement therefore remains necessary.

This is a local stationary-wall slit approximation applied to a moving profile.
It does not solve advective transport of the transverse velocity field or the
full moving-wall Stokes problem. The width is a declared body approximation,
not a fitted damping gain or a measured anatomical value. Exact seals and
invalid paths are rejected without an invented aperture floor. The exterior
lip path, lip-to-mouth attachment, mouth interior, and source/throat end paths
still need their appropriate loss and boundary treatment. In particular,
retiring or transporting viscous memory through disappearance is unfinished.

## Independent spatial check

`lip_viscous_reference.py` evaluates 18 profiles sampled through the two
closing and two opening intervals, with three internal paths per profile and
four quadrature orders: 216 native fixtures. Adaptive SciPy integration of the
exact local slit transfer supplies the independent comparison at 257 frequencies
from 20 Hz to 200 kHz. This range includes the previously identified 81.3 kHz
closing discrepancy. It is an accuracy test of this approximation, not a claim
that a local slit model remains a complete acoustic description over that band.

| Points per path | Maximum relative error in viscous correction from spatial quadrature |
| --- | ---: |
| 4 | 6.342% |
| 8 | 0.01918% |
| 16 | 0.000004468% |
| 32 | 1.27e-12% |

The last path changes aperture substantially. Replacing it with one uniform
slit that matches both integrated inertia and steady Poiseuille resistance
still produces up to 17.59% error in the frequency-dependent correction.
Thus the actual profile matters beyond its total inertia and steady drag.

The finite 24-mode approximation is a separate error source: at orders eight
and above the combined correction error reaches 1.952%, while the error relative
to the complete complex impedance is 0.01496%. Reporting only the latter would
hide error in the smaller viscous component. Native quadrature geometry agrees
with independently integrated centroids and inertia within 1.89e-15 relative.

## Actual moving-airway controls

The same continuous gas/tissue source, thermal throat storage, throat drag,
aperture radiation and actual Bravo geometry are retained. Only internal lip
drag is added. Each case consists of two continued one-millisecond intervals;
these are separately initialized closing and opening controls, not an utterance.

| Case | Before lip drag, 512 to 1024 subdivisions/ms | With eight-point lip drag, 256 to 512 | With eight-point lip drag, 512 to 1024 |
| --- | ---: | ---: | ---: |
| Closing | 39.5388 Pa | 0.41930 Pa | 0.11336 Pa |
| Opening | 2.75778 Pa | 7.46850 Pa | 0.64207 Pa |

Values are maximum endpoint pressure differences across all 14 storage nodes.
They measure numerical refinement, not physical pressure error against a human.
Closing's last lip node now differs by 0.03491 Pa. Opening's largest difference
remains at that node. Nonuniform ratios do not establish uniform second-order
convergence or a converged whole trajectory.

The frozen inviscid tract diagnostic now assigns the largest closing error
energy to a 9.268 kHz mode rather than the former 81.3 kHz lip-end mode. Total
projected error energy falls from 5.48e-11 to 5.99e-14 J. Opening falls from
5.97e-12 to 1.93e-12 J; its dominant projected energy is now at 196.8 Hz.
This diagnostic excludes source and loss-memory feedback and divergence-free
branch circulation. It is neither the coupled damped spectrum nor proof of
inaudibility.

The eight-point runs contain 7,171 accepted steps. The independent audit checks
geometry and replays graph, radiation, thermal and viscous updates using dense
matrix exponentials, driven by the recorded coupled-source flux. It is not an
independent whole-source oracle. Maximum discrepancies are 7.30e-9 Pa in graph
pressure, 3.01e-16 Pa s in viscous memory, 9.15e-20 J in step viscous dissipation,
and 1.09e-21 J in viscous geometry work. The largest interval air-balance residual
is 6.63e-20 m³. Geometry work is nonzero in the moving lip sites; it is retained
separately from nonnegative dissipation.

The sixteen-point controls also pass the independent replay with the same
acceptance limits. At 1024 subdivisions/ms, changing from eight to sixteen
points changes endpoint pressure by at most 0.00006695 Pa in closing and
0.00077834 Pa in opening. These differences are much smaller than the time-step
differences, although the sequential site's splitting error is also part of
this dynamic comparison. The separate frequency-domain integral isolates the
pure spatial quadrature error.

Eighteen focused CTests pass in Release and eighteen under ASAN/UBSAN. All four
moving intervals pass under sanitizers at 256 subdivisions/ms. Twenty-three
historical adapter outputs reproduce their complete JSON exactly, including
the preceding throat-viscous path. The unchanged spatial-pressure recording's
SHA-256 is checked by the [validation collector](continuous_lip_viscous_validation.json).

## Independent exponential arithmetic

The initial sixteen-point opening audit failed its existing exterior-momentum
limit: its full replay differed by 1.543e-12 Pa s against a 1e-12 threshold.
Investigation found roundoff in the independent full-matrix exponential when
many stiff profile modes required repeated scaling and squaring. A selected
matrix comparison with 50-digit mpmath gave 2.46e-12 maximum entry error for
SciPy's ordinary exponential. Switching to a conventional symmetric eigensolve
was worse for that matrix and was not adopted.

The reference now evaluates the small matrix increment `D = exp(A)-I` by a
scaled Taylor series, doubles it using `D <- 2D + D D`, and adds the identity
only at the end. This preserves small slow-mode increments through the scaling
stages. It remains a full-matrix reference independent of the native modal
propagator. Six selected matrices spanning the least, median and greatest
stiffness in closing and opening agree with 50-digit exponentials within
2.79e-16 per entry. This sample is not a global bound. The previously failing
opening replay now differs by 2.09e-13 Pa s in exterior momentum and 7.06e-9 Pa
in graph pressure, passing the unchanged limits. No native physics, native
propagator or acceptance threshold was changed for this correction.

See [the high-precision matrix receipt](lip_viscous_exponential_reference.json).

## Reproduction and remaining work

```sh
cmake --build build --target ducklips_lip_viscous_probe \
    ducklips_planar_inertive_motion_probe
python tools/lip_viscous_reference.py
python tools/planar_continuous_reference.py --radiation \
    --thermal-slit-width 0.03 --viscous-slit-width 0.03 \
    --lip-viscous-width 0.03
python tools/planar_continuous_reference.py --radiation \
    --thermal-slit-width 0.03 --viscous-slit-width 0.03 \
    --lip-viscous-width 0.03 --lip-viscous-order 16 --subdivisions 512 1024 \
    --output out/planar-lip-viscous-16 \
    --report reports/planar_lip_viscous_16_reference.json
python tools/lip_viscous_exponential_reference.py
python tools/viscous_refinement_location.py --input out/planar-lip-viscous \
    --report reports/lip_viscous_refinement_location.json
python tools/radiation_regressions.py --output out/lip-viscous-regressions \
    --include-radiation --include-thermal --include-viscous
```

Use the existing scientific Python environment. The source and approximation
checks are recorded in [the spatial receipt](lip_viscous_reference.json) and
[the moving-airway receipt](planar_lip_viscous_reference.json).
Whole-utterance integration, appropriate treatment of the remaining air paths,
and new human listening remain necessary. Performance optimization remains
parked by user instruction. No speech recording has changed in this step.
