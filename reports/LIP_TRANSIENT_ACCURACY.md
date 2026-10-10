# Separate physical damping from timestep damping in lip opening

The remaining opening discrepancy is concentrated in a fast lip-storage mode.
An independently derived damped mode of the actual endpoint geometry explains
why the two-stage passive update suppresses that motion at the tested steps.
A selectable three-stage update now preserves the same physical response much
more accurately. This changes the time approximation, not the passage geometry,
viscosity, thermal parameters or radiation load. No new speech recording follows
from this change.

## Identify the mode with its physical losses present

The [mode diagnostic](lip_relaxation_mode.json) freezes the final geometry of
the actual opening control and retains every thermal, viscous and radiation
memory. The source port is unforced. It is not the full nonlinear source's
coupled spectrum, and the earlier inviscid wave basis is not used as the answer.

Eliminating the memories at a complex frequency s gives

```
Y_node(s) = s [ C + sum(k lambda / (s+lambda)) ]
Z_branch(s) = s/G + R0 + sum(k lambda s / (s+lambda))
D(s) = diag(Y_node) + B diag(1/Z_branch) B^T
```

The local lip pole is a zero of the pressure operator's Schur complement.
This representation avoids forming a giant fast-memory matrix and subtracting
its large nearly balanced entries to recover a slow response. It retains the
full rational response of the existing finite memory spectrum.

The pole is approximately **-2433.12 + i 392026.14 per second**: a 62,392.90 Hz
oscillation with a physical decay rate of 2433.12 per second. The reconstructed
complex state satisfies both node-air and branch-pressure equations. Its real
part is then used as initial state in fresh native frozen-geometry runs.

For a fixed linear mode, the discrete response is independently known from the
Runge-Kutta stability function. Across both updates and all four tested steps,
the native state agrees with that discrete mode response within 1.71e-20 in
energy coordinates. Thus the isolated mode test distinguishes timestep error
from physical damping without fitting either to an endpoint pressure.

| Steps per millisecond | Extra numerical decay, two stages | Extra numerical decay, three stages |
|---:|---:|---:|
| 256 | 14,445.60 /s | 391.75 /s |
| 512 | 2,269.08 /s | 13.74 /s |
| 1024 | 299.50 /s | 0.4423 /s |
| 2048 | 37.96 /s | 0.01393 /s |

At 1024 steps per millisecond over the tested 0.25 ms duration, excess amplitude
loss falls from 7.214% to 0.01106%. Frequency error falls from -6.109 Hz to
-0.007268 Hz. The frequency being ultrasonic does not prove that the mode is
irrelevant to moving geometry or source feedback; this is not an audibility test.

## Update the approximation, retain the same states

`RelaxationGraphReference` now supports `RelaxationOrder::third` and `::fifth`.
Both use the same equations and eliminate local memory and flow stages before
solving for node pressures. The new three-stage Radau IIA option has fifth order
for smooth fixed-metric equations. It uses the standard
[Radau coefficients and stiff stability properties described by Hairer and Wanner](https://www.unige.ch/~hairer/preprints/coimbra.pdf).

The algebraic-stability matrix is rank one, `M = v v^T`, with
`v=((sqrt(6)-1)/18, -(sqrt(6)+1)/18, 1/9)`. The numerical energy loss is therefore
`dt^2/2 * ||sum(v_i f_i)||_H^2`, nonnegative and reported separately from physical
loss. Derivative differences implement the exact zero sum of v, avoiding a
spurious loss for constant stage derivatives. Memory elimination evaluates
`lambda (I+dt lambda A)^-1` in a bounded form when the relaxation is very stiff.

The actual 14-node tract now needs 42 pressure unknowns for the three-stage
option, while retaining all 1321 passive memory coordinates. It does not solve
a full system over those memories. No performance benchmark or speed claim is
made; performance optimization remains parked by instruction.

`GraphContinuousGlottis::prepare_joint` accepts the order, and the moving driver
selects it with `--joint-relaxation --joint-fifth`. The two-stage path remains
available as the explicit comparison. The geometry is still frozen at the
midpoint and the source still uses its existing implicit work-conjugate update.
The entire moving source/tract method is not fifth order.

## Independent accuracy and actual moving controls

The [generic reference](radau5_reference.json) checks 80 random full-state
systems against independent, unreduced three-stage equations. Maximum scaled
state difference is 1.88e-12. A separate matrix exponential and continuous ODE
agree within 6.81e-17 in energy coordinates. Native transient errors at 16, 32,
64, 128 and 256 steps are 2.827e-8, 8.833e-10, 2.753e-11, 8.585e-13 and
2.676e-14, approximately 32-fold reduction per halving. The stiff constant-
pressure slit tests retain the correct steady flow within 8e-14 relative error.

The [actual moving controls](planar_radau5_reference.json) use the same Bravo
geometry, losses and continuous source as the preceding work. At 64 subdivisions
per millisecond, independently reconstructed geometry and a full sparse stage
system check the native elimination, including memory retirement and geometry
work. This reference uses recorded source flux. Higher resolutions check
work, air balance and refinement, without rerunning the complete stage oracle.

The full-stage comparison is within 3.78e-8 Pa across 518 independently checked
geometry snapshots. Maximum energy-coordinate difference is 6.37e-13. Across
7428 native steps, the combined source/tract air residual is at most 7.24e-20 m3.

| Successive subdivisions per millisecond | Closing, two stages | Closing, three stages | Opening, two stages | Opening, three stages |
|---|---:|---:|---:|---:|
| 256 to 512 | 0.918697 Pa | 0.201510 Pa | 1.51040 Pa | 1.85001 Pa |
| 512 to 1024 | 0.447036 Pa | 0.0403273 Pa | 2.44321 Pa | 0.113044 Pa |

These are successive endpoint differences, not errors against an exact
continuous moving solution. The finest difference is about 11 times smaller
for closing and 22 times smaller for opening. The coarser opening difference
is larger; not every resolution benefits equally. Both complete sampled
trajectories are also compared at common endpoint times in the collected
validation receipt.

The common-time comparison shows the change is not confined to the endpoint.
For the 512-to-1024 pair, opening's maximum sampled pressure difference falls
from 12.0122 to 0.306602 Pa. Closing's maximum lip-section difference falls
from 12.0008 to 0.322919 Pa. However, closing's throat difference remains
9.97345 Pa, nearly unchanged from 9.97698 Pa. The three-stage update therefore
resolves much of the tested lip timing error while leaving a larger transient
near the source. These are native-to-native differences, not absolute errors.

The [continuous source comparison](continuous_radau5_reference.json) uses fresh
native runs and the previously saved independent 164-state ODE, whose artifact
and generator hashes are verified before reuse. Its finest sampled scaled
error remains 0.00399405, almost identical to the two-stage passive result.
Its endpoint scaled error is 0.000182935. This upgrade does not resolve the
remaining source-pressure evolution error.

## Verification and remaining scope

The [collected validation](radau5_validation.json) records four focused Release
and four AddressSanitizer tests covering both passive orders and their source
connection, four actual moving sanitizer controls, and 31 exact historical JSON
regressions. Eighteen preserved third-order fixture files are also replayed;
13 are bit-identical, and the largest scaled numerical difference in the other
five is 8.28e-25. Both orders retain the no-allocation and transactional-rejection
controls. Sanitizer outputs are compared with stated tolerances, not called
bit-identical.

The remaining work includes source gas-pressure accuracy, moving/contact
integration through complete utterances, remaining loss geometry and anatomical
validation. Local numerical agreement does not score human speech. Listening
feedback still belongs to the unchanged spatial-pressure recording: November
improved, Foxtrot still jammed around /tr/.
