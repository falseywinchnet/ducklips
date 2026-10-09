# Air stored between the moving vocal folds

`GasStorage` and `GlottalStorage` add two finite glottal air chambers attached
to the existing four tissue masses. Their pressure forces and displacement
exchange equal work. The lower inlet, inter-chamber passage and upper outlet
have independent flow ports, so incoming air can fill the opening as well as
pass through it.

This addresses a missing component exposed by the
[previous fold-work audit](GRAPH_PHONATION.md). It does not yet replace
`GraphPhonation` or render speech. Prescribed port flows are test inputs;
passage inertia, jet losses and coupling to the tract remain to be solved.
No pronunciation result or new recording is attributed to these controls.

## Constitutive model

The general use of compressible glottal storage and power-conserving coupling
has precedent in [Mora et al., Energy-based fluid–structure model of the vocal
folds](https://pmc.ncbi.nlm.nih.gov/articles/PMC8210679/). The component here is
our own restricted two-chamber implementation, independently checked below;
it is not a reproduction or validation of that paper's complete model.

Let `M` be air mass divided by ambient density (units m³), `V` its current
volume, `r=M/V`, `K=rho0*c0²`, and `gamma` the heat-capacity ratio. Use the
isentropic ideal-gas available energy relative to ambient:

```
H(M,V) = V f(r)
f(r) = K/[gamma(gamma-1)] * [r^gamma - gamma*r + gamma-1]
mu = dH/dM = K/(gamma-1) * [r^(gamma-1)-1]
P = -dH/dV = K/gamma * [r^gamma-1]
```

The mass-flow effort `mu` and wall pressure `P` are distinct derivatives of
one energy. They agree to first order near ambient density; treating them
as interchangeable at finite compression would lose the work identity.
Configured air gives `K=139650 Pa` and `gamma=1.4`, implying an ambient
pressure of `99750 Pa`. These controls use the existing body parameters,
not measured anatomical calibration.

For a finite step, use the density secant and its associated pressure:

```
mu_bar = [f(r1)-f(r0)]/(r1-r0)
P_bar = (r0+r1)*mu_bar/2 - [f(r0)+f(r1)]/2
H1-H0 = mu_bar*(M1-M0) - P_bar*(V1-V0)
```

At equal densities the secant becomes the derivative. Proportional changes
of mass and volume at ambient density produce zero effort, including a
large change of volume. The native evaluation uses stable logarithms and
series near ambient/equal densities; direct subtraction would discard the
small stored energy. The final-volume derivative of the discrete pressure
is also analytic and is checked independently.

## Geometry, mechanics and ports

For each lower/upper tissue pair:

```
V = thickness * [rest_area + length*(x_left+x_right)]
force_on_each_mass = P_bar * length * thickness
```

The chamber has the actual volume implied by this rectangular geometry.
The derivative of that same volume supplies the force area. Changing the
rest opening is a separate actuator and its work is recorded. The tissue
retains the configured masses, linear springs, inter-mass coupling and
damping, including side-specific tension. There is no inherited Bernoulli
force added on top of the chamber pressure.

Port signs are lower inlet, lower-to-upper passage and upper outlet:

```
Delta M_lower = dt*(Q_in-Q_between)
Delta M_upper = dt*(Q_between-Q_out)
port efforts = [mu_lower, mu_upper-mu_lower, -mu_upper]

Delta (gas + mechanical energy) + tissue_damping_loss
    = port_work + rest_opening_actuator_work + external_force_work
```

A four-coordinate implicit solve advances tissue and gas together using
midpoint mechanics and the finite-change efforts. It uses a fixed-size
Jacobian and a positive-volume line search. A feasible Newton starting
guess offsets the changed rest opening; it does not reset the accepted
volume. Failed steps preserve accepted physical state and diagnostics.
Stepping allocates no memory in the native controls.

## Independent checks

The [machine-readable receipt](glottal_storage_reference.json) records exact
source, binary, fixture and output hashes.

- **Thermodynamics:** 210 valid cases against 120-digit arithmetic, covering
  random compressions, nearly equal densities, nearly ambient states, a
  density ratio of `1e-100`, and proportional ambient filling. Three invalid
  cases are rejected. The maximum bulk/density-scaled work-chain residual
  is `5.53e-16`; the maximum scaled pressure-derivative error is `2.19e-14`.
  These are scaled absolute errors, not relative accuracy at zero pressure.
- **Dynamics:** five independent continuous ODE controls: sealed chambers,
  lower-chamber filling, changing rest opening, asymmetric tissue, and all
  three ports active with unequal external tissue forces. Each uses seven
  native refinements over 5 ms, from 120 to 7,680 steps. The independent
  DOP853 oracle is itself refined from 1 µs to 0.5 µs maximum step.
- **Narrow chamber:** an additional native control starts with an upper
  aperture of `1e-8 m²` and volume `5e-12 m³`. It checks positive volume,
  mass and work balance. It has no independent trajectory-convergence claim.
- **Implementation:** six native controls pass under AddressSanitizer and
  UndefinedBehaviorSanitizer; four focused CTests pass, including the existing
  graph, orifice and fold-coupling controls. Invalid net outflow and a NaN
  rest opening leave the accepted state unchanged.

At the finest step, halving the step reduces the normalized endpoint error
by 3.997–4.003 across all five ODE controls, consistent with second order.
The worst endpoint position error is `4.03e-10 m`; the worst velocity error
is `7.75e-6 m/s`. Across all 35 trajectories the largest accumulated work
residual is `6.61e-15 J`. The finest moving-rest case accumulates more solver
roundoff/residual than its coarser runs; a small energy residual does not
establish temporal accuracy. Some coarse runs have substantial phase error
despite preserving the work balance.

The 960-step native controls need at most four Newton iterations per step.
This is a solver observation, not an end-to-end rendering-speed measurement.
Performance optimization remains parked.

## Remaining physical connection

This component accepts positive masses and volumes only. It has no contact
spring, vanishing-cell transport, jet separation, fluid kinetic-energy state
or solved through-flow. The simple rectangular chambers are a constitutive
control, not a calibrated larynx. Their three flow ports must next connect
to momentum/loss elements and the tract with matching power variables.
The existing tract's linear acoustic pressure cannot simply be substituted
for finite-compression enthalpy without checking the energy connection.

The earlier instantaneous fold/orifice work gap remains evidence about the
earlier model. The new component demonstrates a consistent storage and wall
exchange; it has not yet reproduced or resolved that complete coupled
trajectory. Full closure transport, tract losses/radiation and full-utterance
integration remain required before a new Foxtrot listening comparison.

## Reproduction

```sh
cmake --build build --target ducklips_glottal_storage_probe
python tools/glottal_storage_reference.py
ctest --test-dir build --output-on-failure \
  -R 'ducklips_(glottal_storage|graph_phonation|graph_orifice|moving_graph)_reference'
cmake --build .build/asan --target ducklips_glottal_storage_probe
.build/asan/ducklips_glottal_storage_probe
```

The existing `.build/asan` configuration uses
`-fsanitize=address,undefined -fno-omit-frame-pointer`. Outputs are under
`out/glottal-storage/`, including `asan-checks.log` and `ctest.log`.
