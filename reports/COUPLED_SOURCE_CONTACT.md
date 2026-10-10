# Contact in the coupled continuous source

The continuous source can now close and reopen either end while connected to
the joint pressure/flow, thermal, viscous and radiation system. The new path is
opt-in. It has completed free contact cycles and exact-zero opening controls;
a pressure-driven continued trajectory has also crossed a closing and reopening.
This is still mechanics evidence, not a speech recording or a pronunciation score.

## The update at a contact

`GraphContinuousGlottis::advance_joint` now accepts an optional source update
choice. Its default remains the resistance-fitted update. Both source update
choices use the same passive affine pressure response and the same accepted
source/tract port-work calculation.

`JointAccuracyOptions::maximum_contact_step` defaults to zero, disabling the
new transition policy. When enabled, the controller first attempts the normal
update. After failure, and only within the requested time bound, it tries the
existing contact-capable dissipative source update. It accepts that fallback
only if a source opening actually changes between open and closed at the
endpoints. A same-phase fallback is rejected. Ordinary motion therefore does
not silently switch to the more damped discretization.

The complete full-step/two-half-step comparison still applies. The accepted
path consists of the two half steps without extrapolation, and the accepted
receipts include their physical loss and separately identified numerical loss.
Closing and reopening counters count accepted transitions, not rejected trials.
No positive gap floor is inserted. The controller remains transactional over
the entire requested interval and its prepared advancement allocates no memory
in the fixed-geometry contact controls.

The actual moving-mouth adapter exposes the same policy through
`--source-contact-step`, used with `--joint-tolerance`. Its new counters are
printed only when this option is selected. A continued real-geometry closing
interval with the policy enabled gives exactly the previous JSON values after
removing those two added counters; no source contact occurs in that particular
geometry control. This checks integration of the option, not pressure-driven
phonation in a moving mouth.

## A converged final correction was being rejected

The tighter upper zero-area opening control initially failed. Tracing the
source solver showed a different issue from a physical inability to reopen:
the final permitted Newton correction could satisfy the residual test inside
the line search, but the outer loop then rejected it because its correction
budget had ended. At a 9.766 ns trial, the final scaled residual was about
1.845e-9, already below the unchanged absolute threshold.

The budget check now rejects only if that final correction has not converged.
It does not relax the force or trajectory tolerances. The formerly failing
control completes at local tolerance 2.5e-5 and contact-step bound 25 ns. Both
lower and upper zero-area opening cases are registered regression checks.
A trial change to fold finite-difference perturbations was not needed and was
discarded; the final implementation retains the existing perturbation rule.

## Completed controls

`tools/joint_contact_probe.cpp` connects the source to an eight-section,
17.5 cm fixed tube with explicit thermal, viscous and radiation states. This
is a declared tract control, not the complete articulated body. The free
controls start one opening at 1e-7 m² with a closing tissue velocity and zero
reservoir drive, then follow 2 ms of motion with local tolerance 1e-4 and a
100 ns contact-step bound.

| Initially closing end | Accepted substeps | Closures / reopenings | Largest work residual | Largest air residual |
|---|---:|---:|---:|---:|
| Lower | 2828 | 1 / 1 | 1.19e-17 J | 9.16e-22 m³ |
| Upper | 2302 | 1 / 1 | 9.57e-19 J | 4.46e-22 m³ |

Both retain positive source gas and exactly zero exterior flow at a closed
source port. Physical and numerical losses remain nonnegative. The final
solver correction leaves all previously recorded fields of both coarse free
controls unchanged. Five focused checks pass in Release and ASAN/UBSAN builds,
including the tighter exact-zero regressions. Finer whole-cycle runs with the
same current binary also complete one closing and reopening each, using local
tolerance 2.5e-5 and a 25 ns contact-step bound. They take 4714 lower and 3654
upper substeps. Maximum common-time source-pressure changes are 4.28 and
4.12 Pa, while the maximum tract-pressure changes are 0.00145 and 0.0000418 Pa.
Opening-area differences stay below 2.06e-11 and 7.96e-12 m² respectively.
This distinguishes internal fast pressure sensitivity from its effect on the
tract and tissue in these controls. One pair of jointly refined tolerances
and contact bounds is not a proof of convergence or an isolation of either
error source.

## Continuing phonation without resetting memories

The fixed-tube control can now export its complete physical continuation
state and accept it through `--initial-state`. It retains all four tissue
positions and velocities, gas efforts and flows, rest opening, graph charge
and momentum, radiation state, and every thermal and viscous memory.
This format is specific to the declared fixed tube and body parameters.

A 1400 Pa run with a small initial lower-fold displacement completes the first
5 ms. Restoring its endpoint changes the scaled source state by less than
2e-16 and energy by 5.43e-20 J. Source gas masses are identical; graph and all
loss-memory values are identical. The subsequent run has already completed a
pressure-driven closure and reopening. An immutable observed prefix records
that evidence separately from the still-running request to extend the motion
by another 15 ms. Its work residual remains below 1e-12 J and its air residual
below 1e-18 m³ in that prefix. No prescribed glottal waveform drives this run.

The probe also records outlet flow and the radiation branch's relaxation flow
and resistance, so subsequent audio construction can use the physical
radiation port rather than substituting a source signal. No audio has been
published from these values yet.

The [contact validation receipt](joint_contact_validation.json) binds the
completed controls, exact source/binary identities, regression logs, actual
moving-adapter check, and bounded phonation prefix. Longer sustained motion,
articulatory streaming, realistic speech and production efficiency remain
unfinished. The original 10× realtime requirement remains in force; this
accuracy reference is not a performance result.
