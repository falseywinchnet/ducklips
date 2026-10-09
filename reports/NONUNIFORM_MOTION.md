# Continuous nonuniform geometry: reference and speed rejection

The production speech renderer still has the remaining geometry-clock lines
identified in [GEOMETRY_TRANSPORT.md](GEOMETRY_TRANSPORT.md). This work adds
independent moving-network references and a native numerical kernel. It does
not change production PCM, the human reference, or the enhanced display.

## What moves, and what is conserved

For the explicitly assumed lossless canonical network,

```
q = C p                  lambda = M u                 G = M^-1
q_dot = -B G lambda      lambda_dot = B^T C^-1 q
H = (q^T C^-1 q + lambda^T G lambda)/2
H_dot = (-p^T C_dot p + lambda^T G_dot lambda)/2
```

`B` is fixed incidence. Positive diagonal `C` and `G` vary smoothly and
nonuniformly. These equations account for the specified acoustic metric work;
they do not supply mean wall-displacement flow, convection, viscothermal
memory, branch losses, contact/closure mechanics, or empirical human realism.

The five-node, ten-branch fixture has four internal connections, one outlet,
and five grounded inertances. Its nominal frequencies are approximately
651, 1206, 2001, 2460 and 3264 Hz. Independent canonical and energy-coordinate
integrations agree to 7.04e-10 relative state error. Integrated work closes the
energy balance to 2.68e-11 in the low-mode case. The uniform-area limit agrees
with the prior exact propagator to 2.34e-10.

## A frozen energy subspace can mislead

With five pressure and five flow modes at the nominal shape, the fixed
orthogonal energy reduction has 17.0% RMS state error during the nonuniform
motion. A union of two poses needs 14 states and reaches 5.60e-5 RMS error;
three poses fill the full 15-state space. Each still satisfies its own
projected power identity. Energy balance alone does not establish trajectory
accuracy.

This does **not** mean that the lossless moving network needs five new
physical circulation variables. Its momentum equation gives a simpler exact
invariant when initial momentum is a gradient:

```
lambda = B^T eta
eta_dot = C^-1 q
q_dot = -B G B^T eta
```

There are only twice as many dynamic states as nodes, regardless of branch
count. The 10-state potential representation agrees with the full reference
to 8.60e-10. Canonical circulation stays below 9.70e-16. The apparent energy
in a nominal frozen circulation subspace arises because multiplication by
`sqrt(G(t))` moves the energy-coordinate subspace. It is not new canonical
circulation. Arbitrary series forcing and branch-specific damping generally
break the gradient invariant and require additional treatment.

See [nonuniform_metric_reference.json](nonuniform_metric_reference.json).

## Native direct reference: accurate convergence, insufficient speed margin

[MovingTreeReference](../include/ducklips/moving_tree.hpp) implements two-stage
Gauss integration of the full charge/momentum system. Eliminating branch
stage variables leaves a tree of 2x2 node blocks. Leaf elimination costs
O(nodes + branches), including fresh metric coefficients and factorization
at every step. The prepared step allocates no storage. Ground and parallel
branches are supported; cyclic or disconnected internal topology is rejected.
There is no per-block reprojection and no output smoothing or notch filter.

The native tests check stationary quadratic energy conservation, a parallel
branch network, fourth-order convergence against the exact moving uniform
mode, and rejection of cyclic topology. Independent nonuniform reference
checks use a low pressure mode, that mode with added canonical circulation,
and the highest nominal mode.

| Substeps per 48 kHz output sample | Low-mode maximum relative state error | Highest-mode error | 96-node kernel speed, three runs |
|---|---:|---:|---:|
| 1 | 9.11e-6 | 0.0630 | 15.99–16.33x |
| 2 | 5.69e-7 | 0.00398 | 8.09–8.11x |
| 4 | 3.55e-8 | 0.000249 | 4.03–4.05x |

Errors cover 40 ms of nonuniform motion on the five-node reference; speed
uses a separate 96-node, 191-branch, 287-state network for five simulated
seconds. Timings include stage-metric filling and every factorization, but
exclude the voice source, detailed losses, articulation and microphone.
They are not whole-speech rendering rates. Added-circulation errors are
similar to the low-mode values. Refinement ratios are 15.8–16.1, consistent
with fourth order.

At 48 kHz a stationary 8 kHz mode advances at approximately 7987.51 Hz under
this Gauss step. At 12 kHz it advances at 11913.10 Hz. This phase bias is the
reason the high-mode trajectory error cannot be inferred from the low-mode
success. Stationary energy conservation does not remove numerical dispersion.

**The accurate substepped configurations fail the 10x floor.** The verification
script deliberately returns nonzero after writing its receipt. The one-step
16x timing is not acceptance of the kernel as a speech replacement. Receipts,
source/binary hashes and phase-bias values are in
[moving_tree_reference.json](moving_tree_reference.json).

## Better representation: exact oscillation, slow geometric coupling

The speed rejection motivates separating the rapid oscillation from geometry
transport instead of numerically subdividing every acoustic cycle. The
continuous modal equations can be derived directly from the exact canonical
potential system, without an arbitrary state-transfer invariant.

Let `S=C^-1/2 B G B^T C^-1/2`, `S U=U Omega^2`, `R=C^1/2 U`, and
`F=R^-1 R_dot`. Transform `q=R a`, `eta=R^-T h`, then set `b=Omega h`:

```
a_dot = -Omega b - F a
b_dot =  Omega a + (Omega_dot Omega^-1 + Omega F^T Omega^-1) b
H = (a^T a + b^T b)/2
```

For separated eigenvalues,

```
D = U^T S_dot U
Q_ij = D_ij / (omega_j^2 - omega_i^2), i != j; Q_ii = 0
F = U^T diag(C_dot/C) U / 2 + Q
omega_dot_i = D_ii / (2 omega_i)
```

Uniform area motion reduces exactly to the previously verified opposite-sign
metric terms, `a_dot=-omega*b-r*a/2`, `b_dot=omega*a+r*b/2`. The fast diagonal
oscillation already has an exact inexpensive propagator. The remaining
geometry operator in this fixture has maximum norm about 59.2/s, whereas
acoustic angular frequencies span approximately 3586–25052 rad/s. That scale
separation is a concrete reason to investigate this representation.

The new [continuous modal oracle](../tools/moving_modal_reference.py) verifies
the equations against the full nonuniform network: maximum relative state
errors are 1.02e-9 for the low-mode start and 1.33e-9 for the highest-mode
start. Its instantaneous wall-power residual is below 5.93e-13. Results are in [moving_modal_reference.json](moving_modal_reference.json).
This oracle is not a native fast implementation. Near-degenerate modes,
closures, branch-specific losses and sources remain unsolved here. It rejects
ambiguous eigenvector orientation or small eigenvalue gaps rather than
silently using an unstable connection formula.

## Reproduce

```sh
cmake --build build
ctest --test-dir build --output-on-failure
python tools/nonuniform_metric_reference.py --out out/geometry_clock/nonuniform
python tools/verify_moving_tree.py --out out/geometry_clock/moving-tree
# The preceding command intentionally fails its 10x speed gate.
python tools/moving_modal_reference.py --out out/geometry_clock/moving-modal.json
```

The next implementation task is to exploit exact oscillation with continuous
geometry coupling, then restore and verify physical losses/source coupling
and measure complete speech throughput. The remaining horizontal lines are
not claimed removed by these numerical reference experiments.

Validation on this revision: core Release 4/4, core ASan/UBSan 4/4, and surface
Release 6/6 CTest suites pass. A fresh production Alfa render remains byte
identical to `out/geometry_clock/alfa.wav`; it has zero nonfinite or clipped
samples and renders 1.31 seconds in 0.079 seconds in this single-run check.
This is a PCM regression check, not a replacement for the previous multi-case
performance receipt. The physical reference prototype does not change the
speech default.
