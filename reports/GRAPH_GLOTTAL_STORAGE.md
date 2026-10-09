# Pressure-driven glottal chambers connected to the tract

`GraphGlottalStorage` connects the two gas/tissue chambers to a reservoir and
an acoustic graph through three inertive, dissipative passages. All three
flows are solved from pressure/enthalpy differences and retained momentum;
they are no longer supplied test inputs. The tract's average pressure feeds
back into the same implicit step.

This is a coupled engine component, not a new speech renderer. Passage
coefficients are fixed, independently configured controls. Their dependence
on moving fold geometry, jet separation and full contact remains unfinished.
The [previous chamber component](GLOTTAL_STORAGE.md) supplies gas storage and
matching tissue work. The older algebraic `GraphPhonation` is not replaced.

## Connection and work

The three flows `Q0,Q1,Q2` are ambient-reference mass volume per second:
reservoir to lower chamber, lower to upper chamber, and upper chamber to
tract. They need not match while chamber contents change. Each passage has
fixed positive inertance `L`, nonnegative linear resistance `R`, and
nonnegative quadratic loss coefficient `B`:

```
passage energy = lambda²/(2 L)
Q_bar = (lambda_new+lambda_old)/(2 L)
(lambda_new-lambda_old)/dt = drive - R Q_bar - B Q_bar |Q_bar|

drive0 = reservoir_mu - lower_mu_bar
drive1 = lower_mu_bar - upper_mu_bar
drive2 = upper_mu_bar - tract_mu_bar
```

Reservoir pressure is converted to enthalpy using the same gas equation of
state as the chambers. Supplying 700 Pa therefore does not mean assigning
700 Pa to both distinct gas efforts. The graph uses the linearized available
energy `q²/(2 C)` for charge `q=M-V`; its conjugate mass-flow effort is `q/C`.
This equals its acoustic pressure in the linear approximation. Connecting
those efforts preserves port work without asserting that the downstream
linear graph is a finite-compression ideal-gas model.

The graph is evaluated for free motion and unit source flow, yielding its
average load `mu_free + Z Q2`. The native solver finds all three passage
flows together, evaluating the gas/tissue solve for each candidate and
including that load. The 3×3 flow Jacobian is obtained by centered finite
differences; the inner tissue solve retains its analytic 4×4 Jacobian.
Both use bounded iteration and line searches. Accepted states and
receipts survive rejection; native stepping allocates no memory.

The source receipt checks:

```
Delta (chamber gas + tissue + passage kinetic energy)
  + tissue damping + passage dissipation + work delivered to tract
  = reservoir work + rest-opening actuator work + external tissue work
```

For the frozen graph controls, its energy change equals the received work.
The ground outlet has zero effort and therefore carries mass but no available
energy relative to that reference. Graph geometry work is not included in
this frozen-network theorem; the existing moving graph has its own metric
work requirements.

Energy-based compressible-flow vocal-fold models also distinguish relative
enthalpy, pressure and kinetic storage; see [Mora et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8210679/).
The present fixed-passage circuit is a restricted coupling implementation,
not a reproduction of their spatial fluid equations or phonation results.

## Checks and findings

The [independent receipt](graph_glottal_reference.json) hashes code, native
binary and outputs. The continuous oracle uses instantaneous gas pressure,
mechanical forces, passage-flow acceleration, and graph pressure/flow
coordinates. It does not use the native discrete effort or implicit solver.

Four 5 ms controls cover forward drive, reverse drive, undriven relaxation
with stored pressure/momentum, and a changing rest opening. Each has eight
native refinements (120–15,360 steps), for **32 coupled trajectories**.
The DOP853 oracle is separately refined from a 1 µs to a 0.5 µs maximum step.

At the finest step:

| Control | Error reduction on halving step | Maximum tissue-position error | Maximum tract-pressure error |
| --- | ---: | ---: | ---: |
| Forward | 3.970 | 5.74e-12 m | 4.27e-5 Pa |
| Reverse | 3.992 | 1.03e-12 m | 5.36e-6 Pa |
| Relaxation | 4.008 | 4.26e-11 m | 3.21e-5 Pa |
| Rest-opening motion | 3.970 | 5.75e-12 m | 4.29e-5 Pa |

Across all refinements, maximum accumulated mass-balance error is
`2.65e-21 m³`, maximum accumulated work-balance error is `1.86e-14 J`, and
maximum passage-law residual is `1.50e-7 Pa`. These are numerical checks of
these equations, not anatomical or perceptual validation.

A separate 100 ms forward control (19,200 steps) finishes with maximum mass
residual `6.44e-20 m³` and maximum work residual `5.90e-17 J`. During that
control the inlet and outlet flows differ by up to `1.003e-5 m³/s`, showing
why one algebraic through-flow cannot account for chamber filling. Its
long-duration control is not an independent continuous-trajectory audit.

Four native coupled controls and six underlying storage controls pass
AddressSanitizer/UndefinedBehaviorSanitizer. Five focused CTests pass. Native
checks include reverse flow, undriven energy decrease, state continuation,
invalid-input rejection and allocation-free stepping. The underlying
210 thermodynamic cases and 35 chamber/tissue continuous comparisons were
rerun after the solver correction below.

### Refinement exposed a tissue-solver error

The initial fine-step coupled comparison failed its convergence criterion.
The inner tissue solver scaled its tolerance with the magnitudes of two
large opposing inertial terms. Their individual magnitudes grow as the time
step shrinks even when the remaining physical force does not. Moreover,
recovering velocity from the difference of two nearly equal absolute
positions amplified position subtraction error.

The solver now advances displacement directly and scales convergence with
the remaining physical forces, separately allowing for floating-point
velocity subtraction. This changes the numerical representation, not tissue
parameters or the constitutive law. The 15,360-step forward endpoint position
error fell from `9.68e-11 m` to `5.74e-12 m`; accumulated work residual fell
from `1.53e-12 J` to `1.87e-18 J`. Final convergence checks pass for all four
controls. Coarse phase error remains distinct from work conservation.

## Remaining source physics

Fixed coefficients verify the connection and momentum exchange but do not
make the glottal passages narrow with the folds. The next physical requirement
is to derive passage kinetic energy and loss from the same moving geometry
that sets chamber volume. Updating `L` without its associated work/forces
would reintroduce an unaccounted energy exchange. Full disappearance/contact
also needs mass and momentum transport rather than an area or volume floor.
The tests here establish neither self-sustained phonation nor a physiological
jet law. The configured coefficients are controls, not fitted anatomy.

The coupling API accepts the moving graph's stage metrics, but this audit
uses frozen graph metrics; replaying the native mouth/lip geometry remains
an integration check. No new audio, enhanced field, NATO score or speed claim
is attached to these results. Performance optimization stays parked.

## Reproduction

```sh
cmake --build build --target ducklips_graph_glottal_probe ducklips_glottal_storage_probe
python tools/graph_glottal_reference.py
python tools/glottal_storage_reference.py
build/ducklips_graph_glottal_probe --trajectory 19200 0.1 forward
ctest --test-dir build --output-on-failure \
  -R 'ducklips_(graph_glottal|glottal_storage|graph_phonation|graph_orifice|moving_graph)_reference'
```

Outputs are under `out/graph-glottal/`. Supplemental native/sanitizer output
hashes are in [the validation receipt](graph_glottal_validation.json).
