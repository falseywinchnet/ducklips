> **OWNER RULE (2026-10-20): NO CALIBRATION, EVER.** No fine-tuning, fitting, per-word patches, per-speaker profile fits, or tuned constants. Progress comes only from structure.

# The Jacobian we depend on: its cost, its mathematics, a shared and checkpointed representation (2026-10-24)

Labels as AGENTS.md: GIVEN, OBSERVED, MEASURED, HYPOTHESIS, CANDIDATE,
REJECTED, DECIDED.

Scope (GIVEN, owner direction 2026-10-24, NOTES_pyramid.md step 45): the
puppet's Levenberg-Marquardt Jacobian over the equation error is central
finite differences over full renders, 6-8 columns a word and about 1,000
residuals a Jacobian for the voice problem. Refine it; do the formal
mathematics toward a better representation (the waveguide's adjoint,
checkpointed renders, the chain to the puppet's parameters). Module:
`analysis/adjoint.py`. No change to `puppet.py`, `voice.py`,
`outside_in.py`, `fields_word.py`, `pyramid_word.py`, `NOTES_pyramid.md`
(the main session's). No engine change was needed: the engine's snapshot
and resume (`engine.Run`, `tract_run_copy`, `Settings.divergence`,
NOTES_engine.md) already exist behind the ctypes interface, as do the peel's
(`observer.implied_record` / `implied_resume`) and the geometry's
(`reverse_core.moving_geometry_record`, NOTES_kernels.md); the puppet's
`equation_error` used none of them (OBSERVED).

Environment for every measurement: the M4 Mac mini (10 cores),
`~/miniforge3/bin/python3` (numpy 2.3.5, numba 0.63.1, scipy 1.16.3), the
engine library `build_engine/libtract_engine.dylib` of 2026-10-22, shared
with the main session's `puppet.py` (8 workers) and an unrelated 8-worker
job: load average 19-28 throughout. Times are wall and CPU seconds of this
process; ratios within one run are the measurement, absolute times are
inflated by the load (step 45's 0.29-0.50 s residual read 0.55-0.72 s here).
Word: Papa, `P@lips!,AA,P@lips!,AH`, `JETBOUND=1 SEALBLEND=1`: 51,840 steps
(0.540 s), 540 ticks, 46 nodes, window [1920, 50880), seal 0.01, no jet;
true events (s) P-0.100 AA+0.150 AA-0.220 P+0.290 P-0.360 AH+0.420; the
start the hinted times (+-7.3 ms). Commands: `adjoint.py profile | check |
columns | lags | tangent | sweep SPECS`; logs under `/tmp/esf/adjoint_*.log`
on the Mini.

## 1. The cost structure (MEASURED)

One plain residual (`puppet.PuppetProblem.residual`, cProfile, Papa at the
hinted start, -12.7 dB; wall 0.688 s, CPU 0.619 s, median of 5, load 28):

| piece | s | share |
|---|---|---|
| three `engine.render` calls (`pyramid_word.FieldFlow`: the free folds of each of the two closure windows, then the forced folds) | 0.526 (0.175 each) | 75% |
| `observer._implied_run` own time (closures, window bookkeeping) | 0.040 | 6% |
| the peels (`StatesPeel.advance`, three calls) | 0.062 | 9% |
| `reverse_core.moving_geometry` | 0.044 | 6% |
| `puppet.command` (49 `node_progress`, `read_postures`) | 0.005 | 1% |

One column (two residuals, secant 10 ms): wall 1.356 s, CPU 1.247 s.

Reading: under load the C++ render is the cost (step 45 measured 0.20 of
0.29-0.50 s at a lower load; here 0.53 of 0.69). The three renders are of
one causal system over one word; two of them are read only inside a closure
window (`observer._implied_run` reads `folds[hidden_from:hidden_to]`), the
third to the end.

## 2. The step as an operator (OBSERVED: `observer._forward_steps`, the replica of `Tube::step_lining`, `step_air`, `scatter_links` for this configuration; `src/tube_step.cpp`)

Indices: steps t = 0..T-1, nodes i = 0..n-1, links i = 0..n-2 between
nodes i and i+1, control ticks of 96 steps. Per step the geometry G_t
(`moving_geometry`'s GEODUMP row) gives each node's area A_i, the frame's
pushed-air flow Phi_i, each link's viscous resistance R_i and each node's
lining constants (surface S_i, cover mass m_i, body mass M_i, hoop and
coupling and support stiffnesses, cover damping c_i, body damping C_i).
State x_t: per link the arriving waves a+_i (toward node i+1), a-_i (toward
node i) and the smoothed outgoing waves s+_i, s-_i; per node the lining (q_i,
v_i, Q_i, V_i: cover and body displacement and velocity) and the pressure
p_i; the radiation's inertance flow F_L and last pressure p_prev. Input u_t
= U_in(t), the flow entering node 0 (the folds'). Outputs U_out(t) (the
sound) and p_0(t).

One step t, in order:

(L) Lining, implicit midpoint on the cover with the body by backward Euler
(constants of G_{t-1}; per node, AXIAL = 0 so no coupling between nodes):

    q' = [m (2 q + h v) + (h^2/2) k_c Q] / (m + (h^2/4)(k_h + k_c)) - q
    v' = 2 (q' - q) / h - v
    V' = [M V + h (k_c q' + k_s o - (k_c + k_s) Q + c v')] / (M + h (C + c) + h^2 (k_c + k_s))
    Q' = Q + h V'

Linear in (q, v, Q, V) given G.

(J) Junction at node i (admittances y_i = (A_i + A_{i+1}) / (2 rho c)):

    K_i = 2 y_{i-1} a+_{i-1} + 2 y_i a-_i - Phi_i          (node 0: U_in + 2 y_0 a-_0 - Phi_0;
                                                             node n-1: 2 y_{n-2} a+_{n-2} - F_L - g p_prev - Phi_{n-1})
    Y_i = y_{i-1} + y_i                                      (node 0: y_0; node n-1: y_{n-2} + 1/R_rad + g)
    rho_i = S_i^2 / Y_i + c_i                                (the lining's resistance at the junction)
    v_i' = [(m_i - h rho_i / 2) v_i + h (S_i K_i / Y_i + c_i V_i)] / (m_i + h rho_i / 2)
    p_i = (K_i - S_i (v_i + v_i') / 2) / Y_i

with R_rad = 128 rho c / (9 pi^2 A_{n-1}), L_rad = 8 rho / (3 pi sqrt(pi
A_{n-1})) + rho e / A_{n-1} (e the tube past the last node), g = h / (2
L_rad). Affine in (a+, a-, U_in, v, V) given G.

(R) Radiation: F_L <- F_L + g (p_{n-1} + p_prev); p_prev <- p_{n-1}; U_out =
p_{n-1} / R_rad + F_L.

(H) Hop low-pass (k = HOP = 0.02) on the outgoing waves:

    s+_i <- (1 - k)(p_i - a-_i) + k s+_i
    s-_i <- (1 - k)(p_{i+1} - a+_i) + k s-_i

(S) Series element on link i (z_i = 1 / y_i, l_i = 2 z_i + R_i, kappa_i the
Borda-Carnot coefficient (rho/2)(1/A_narrow - 1/A_wide)^2 in the direction
the air expands, phi_i in [0, 1] the seal's open share of the step):

    d_i = 2 (s+_i - s-_i)                      (+ the jet's series source, none here)
    U_i = phi_i sgn(d_i) 2 |d_i| / (l_i + sqrt(l_i^2 + 4 kappa_i |d_i|))
    a+_i <- s-_i + z_i U_i,   a-_i <- s+_i - z_i U_i        (sealed, phi = 0: a+ <- s-, a- <- s+)

The seal share phi_i (`Tube::set_seal_blend`, step 9): 1 while the link's
narrower area is above the seal, 0 while below, and in the step it crosses,
the share of the step on the open side with the area linear over the step:
(least - seal) / (least - before) opening, (before - seal) / (before -
least) closing.

Linearity (OBSERVED, correcting the brief's premise): stages (L), (J), (R),
(H) are linear in the state given G; stage (S) is not, through kappa (the
kinetic loss, `kinetic=True` in every path), and the engine's jet when
present. With kappa = 0, U_i = d_i / l_i and the whole step is linear. The
tangent below linearizes (S) on the trajectory.

The geometry Gamma: G_t = Gamma_t(a_0, ..., a_t), a per-node recurrence
(`reverse_core._geometry_steps`, `Tube::refresh_geometry`): the target
radius r* = sqrt(max(a, 1e-4) / pi) ramped over 12 steps at each 12-step
control tick, followed by a critically damped 160 Hz second-order follower
(a sleep branch below 1e-9 cm and 1e-6 cm/s), the lining constants compiled
from the radius (surface ~ r, hoop ~ 1/r, masses ~ r, dampings ~ sqrt of
products), Phi_i the follower's velocity times the surface through two 250
Hz poles, R_i the viscous law of the link's mean area (max of 0.03 A^-2.6
and 0.0065 sqrt(2250) / A) glided over 16 steps unless either node moves
faster than 0.5 cm/s (then exact). Smooth except the law's max, the fast
switch, the sleep branch and the area floor.

## 3. Tangent (forward mode)

For a direction d theta of the puppet's unknowns:

    d a_r      = (d a / d theta) d theta                                       (section 6: analytic)
    d G_s      = sum_{r <= s} Q_{s,r} d a_r                                     (Gamma's tangent, per node)
    d x_{t+1}  = A_t d x_t + b_t d u_t + P_{t,t} d G_t + P_{t,t-1} d G_{t-1}     (the step's tangent)
    d U_out(t) = (1 / R_rad) d p_{n-1} + d F_L + p_{n-1} d(1 / R_rad)

with A_t = dF_t/dx on the trajectory: stages (L), (J), (R), (H) as they are
(they are linear), stage (S) replaced by its linearization

    d U_i = gamma_i(t) d d_i,   gamma_i = dU/dd = phi 2 / (l + r) [1 - 2 kappa |d| / (r (l + r))],   r = sqrt(l^2 + 4 kappa |d|)

(continuous across d = 0, where gamma = phi / l from both sides although
kappa switches direction there), and d a+_i = d s-_i + z_i d U_i + U_i d z_i,
d a-_i = d s+_i - z_i d U_i - U_i d z_i. The parametric part P collects the
explicit derivatives of (J), (R), (S) with respect to the geometry
quantities of the step: in (J), dK_i/dy_i = 2 a-_i and dK_{i+1}/dy_i = 2
a+_i, dY/dy = 1, dK_i/dPhi_i = -1, and through rho_i and v_i' the lining's;
in (S), dU/dl = -2 d / (r (l + r)) (so dU/dR_i = dU/dl, dU/dz_i = 2 dU/dl),
dU/dkappa = -4 d |d| / (r (l + r)^2), dU/dphi = U / phi; in (R) and node
n-1's (J) the radiation load's dR_rad/dA and dL_rad/dA. The seal share's
derivative with respect to the areas in the crossing step is of size 1 /
(least - before): the step's slope across one engine step of event motion,
the "ramp" of step 9 (4-20 per 2e-6 phones), finite but a hundred times the
slope elsewhere.

Cost: one tangent lane is the primal's arithmetic about twice over; N lanes
carried together (`ad/dual.hpp`'s Dual<N>) amortize the primal's memory
traffic and the geometry's work (only the nodes a lane moves). MEASURED on
the keyframe family without seals (NOTES_gradients.md, AA,IY,UW, 143
columns): 4.5 s CPU against 43.6 s for forward differences, 0.031 s a
column against two residuals; the puppet's family adds the windows
(section 5) and the fold targets' tangent (`ff_tangent`'s `targets_dot`,
untested there).

## 4. Adjoint (reverse mode)

For a scalar Lambda = sum_t v_t . (outputs at t) (v the residual's weights:
v = r for the gradient of |r|^2 / 2, or any one vector), the adjoint state
lambda_t obeys, backward from lambda_T = 0,

    lambda_t    = A_t^T lambda_{t+1} + (d outputs_t / d x_t)^T v_t
    dLambda/du_t = b_t^T lambda_{t+1}
    dLambda/dG_s = P_{s,s}^T lambda_{s+1} + P_{s+1,s}^T lambda_{s+2}  =: g_s
    dLambda/da_r = sum_{s >= r} Q_{s,r}^T g_s                        (Gamma's adjoint, backward in time, per node)
    dLambda/dtheta = sum_r (d a_r / d theta)^T dLambda/da_r

one backward pass for every area at every step and every unknown at once.

Reciprocity (what A_t^T is). In power-normalized waves (a scaled by sqrt y)
the junction's lossless scattering, from (a+_{i-1}, a-_i) to the outgoing
(p - a+_{i-1}, p - a-_i), has the matrix

    [ 2 y_{i-1}/Y - 1        2 sqrt(y_{i-1} y_i)/Y ]
    [ 2 sqrt(y_{i-1} y_i)/Y  2 y_i/Y - 1           ]

symmetric; the lining's scalar implicit exchange at the junction is a
reciprocal one-port (a scalar admittance in series with the junction) and
keeps the symmetry; the linearized series element maps (s+_i, s-_i) to
(a+_i, a-_i) by

    [ 2 z gamma      1 - 2 z gamma ]
    [ 1 - 2 z gamma  2 z gamma     ]

symmetric already in the plain waves (the same z on both sides); the hop
one-pole s_t = (1 - k) w_t + k s_{t-1}, as an operator on sequences,
transposes to mu_t = (1 - k) lambda_t + k mu_{t+1}, the same one-pole run
anticausally; the radiation's trapezoid and the lining's two implicit steps
transpose likewise (their coefficients are symmetric per step, their
recursions reverse). A delay (a wave crossing a link in a step) transposes
to a delay in reversed time. Therefore A_t^T is the step of the same
waveguide with time reversed and the wave directions exchanged (an adjoint
wave that in the forward run arrived at node i+1 leaves node i+1 toward
node i), the one-poles run anticausally, and the only trajectory-dependent
coefficients gamma_i(t) and phi_i(t) read from the forward run: Tellegen's
adjoint network for a reciprocal, time-varying, mildly nonlinear line. The
adjoint pass needs d_i(t) (for gamma) and, for the parametric contractions
P^T lambda, the forward waves and pressures of every step: either stored
(T x 6 n doubles, 114 MB for Papa) or recomputed block by block from the
engine's own snapshots every 960 steps (a second forward pass).

What the folds' linearization needs. The cord (`src/cord.cpp`, VTL 2.4's
triangular glottis, GPL mathematics only) is x^c_{t+1} = f(x^c_t, w_t,
A_0(t), pi_t; zeta_t): w_t the returning wave from the vocal tube (forced:
FORCEWAVE's value; free: `vocal_.outgoing_left()`), A_0 the first node's
area, pi_t the posture target (FOLDTICKS, glided over 1 ms by two one-pole
stages), zeta the hypothesis-independent drives (lungs' ramp and
compensation, the trachea's 33-node Tube, variation, tremor). OBSERVED
(`VoiceBody::step`, NOTES_gradients.md section 1): w and A_0 are the only
vocal-tube quantities the folds read. Its tangent needs

- df/dx^c: the two masses' backward-Euler step (a 2 x 2 implicit system a
  step: its Jacobian is that system's), the contact fraction and the
  quadratic first-touch force (C1), the contact damping growing with the
  fraction, the Bernoulli flow with the glottal inertance (a quadratic
  solve, smooth; a kink at `max(area, floor)`), the posture glide (linear),
  the trachea's Tube (section 2's structure, 33 nodes, driven by the cord's
  wave below) — about 60 state doubles a step plus the trachea's;
- df/dw: w enters the supraglottal pressure as 2 w - Z_0 U (the two-port
  with the first node, Z_0 = rho c / A_0): one linear entry into the flow
  solve and the driving pressures;
- df/dA_0: through Z_0 only (and `forced_flow` = (sent - w) A_0 / rho c);
- df/dpi: the glide's gains, linear;
- not differentiable: the ear (`hear_pitch`, PITCHFB 0.08 per closure, the
  heard period a whole number of steps: a closure instant crossing a step
  boundary moves the tension by a quantum; NOTES_gradients.md section 6
  measured jumps every ~4 us along a time parameter), the triangular
  glottis's contact fraction when a rest gap goes negative (a 0 -> 1 jump).

The cord's adjoint is lambda^c_t = (df/dx^c)^T lambda^c_{t+1} + (the
revealed residual's weight at t), with dLambda/dw_t = (df/dw)^T
lambda^c_{t+1} feeding the peel's adjoint (section 5) and dLambda/dA_0(t)
feeding the geometry's. It needs the cord's trajectory stored (T x ~60
doubles, 25 MB) or recomputed from checkpoints.

## 5. The peel and the hidden windows

The peel (`reverse_core.SidePeel._side_link`, `observer._states_link`) is
the inverse of section 2 layer by layer, j = n-1 down to 0. With alpha_j(t)
the wave arriving at node j+1 from link j and sigma_j(t) the smoothed wave
node j+1 sent into link j (its one-pole run causally from node j+1's known
pressures), layer j is, for every t,

    U_j(t)       = (alpha_j(t) - sigma_j(t)) / z_j                               (the series element inverted: linear)
    eps_j(t)     = sigma_j + (z_j + R_j / 2) U_j + (kappa_j / 2) U_j |U_j|       (the wave that entered from node j: s+_j; explicit, no root)
    a-_j(t)      = eps_j(t) - z_j U_j(t)
    p_j(t)       = (eps_j(t+1) - k eps_j(t)) / (1 - k) + a-_j(t)                 (the hop one-pole inverted: one step ahead)
    alpha_{j-1}(t) = [p_j Y_j + S_j vbar_j - 2 y_j a-_j + Phi_j] / (2 y_{j-1})    (node j's junction solved for the other wave; the lining stepped causally on p_j)

and at j = 0 the last line is U_in(t) itself. Each layer is a causal
recurrence (sigma, the lining) with a one-step-ahead read, so the whole peel
reads n steps ahead (the lookahead). Given G, the peel is linear in the
sound except for the kappa term (explicit, differentiable: d eps = d sigma +
(z + R/2 + kappa |U|) dU + U dz + (U/2) dR + (U|U|/2) dkappa; note z + R/2
+ kappa |U| = 1 / (2 gamma) at the same point, the chain of inverses).

Tangent of the peel: the same layers in d, from the lips to the glottis.
Adjoint of the peel: the transposed layers in the reverse dependency order,
j = 0 up to n-1 (from the glottis toward the lips, the forward tube's
direction of information), each layer's recurrences run anticausally and
its one-step-ahead read becoming a one-step-behind write. One pass, like
the tube's. The peel's adjoint turns the revealed residual's weights at
node 0 into weights on every (G, t): dLambda/dG from the peel's side; the
cord's adjoint gives dLambda/dw_t, and w = (p_0 - U rho c / A_0) / 2 passes
it back into the peel's outputs (dLambda/dp_0 = lambda^w / 2, dLambda/dU
-= lambda^w rho c / (2 A_0), dLambda/dA_0 += lambda^w rho c U / (2 A_0^2)).

Hidden windows (`observer._implied_run`): a closure's stretch [h_k, h'_k)
(from the first sealed step less the seal node's index to the last sealed
step plus two) is not peeled; the tube runs forward from the state the peel
recorded at h_k - 1 with the folds free (the engine's own run, NaN forced
waves), and the residual there is the modelled sound against the sound.
Its tangent is the coupled tangent of section 3 and section 4's cord over
the window (the returning wave the tube's own: d w = d outgoing_left),
started from the peel's tangent of the recorded state, and its final state
seeds the next peel segment's initial condition; its adjoint is the coupled
adjoint backward over the window, its terminal lambda at h_k entering the
preceding peel's adjoint at the record. The window's edges are whole steps
chosen by thresholds (`observer.closures`: seal, HIDE_AREA): held fixed by
the derivative; a parameter that moves an edge by a step moves a residual
row from one definition to the other: a jump of |(D - U) - (M - y)| at that
step.

## 6. The chain to the puppet's parameters (analytic; `command_tangent`, `physical_tangent`)

Column c's clock (`puppet.node_progress`): on transition k, s = k + (1 -
cos(pi f)) / 2 with f = (t - L'_k) / D_k, L'_k = L_k + lag_c the column's own
leaving, D_k = max(A_{k+1} - L_k + extra_c, one step) its move's duration;
a hold has s constant. So

    ds/dtheta = (pi / 2) sin(pi f) df/dtheta,
    df/dL_k = (f - 1) / D_k,   df/dA_{k+1} = -f / D_k,   df/dlag_c = -1 / D_k,   df/dextra_c = -f / D_k

(where D_k is clipped at one step, df/dA and df/dextra are 0 and df/dL = -1
/ D_k). Then `fields.read_postures` with the identity space field: command
= (1 - share) rows[lower] + share rows[upper], so d command[:, c] / d s_c =
rows[upper, c] - rows[lower, c], piecewise constant in s with kinks at the
targets (integer s), where sin(pi f) = 0: the product is continuous. The
side columns (port, folds) read the last column's clock. Last, the physical
clip (`pyramid_word.physical`, max with the sealed area 1e-4 cm^2, with 0
for the port and the rests): C0, the derivative zero where the clip holds
and `slope x ds/dtheta` elsewhere, with a jump at the crossing (step 43's
"meter": the seal crossing is where the smooth command meets the clip).
Where a stretched or lagged move overruns the next hold the later
transition wins (as in `node_progress`), and the derivative follows the
winning branch.

Verification (MEASURED, `adjoint.py tangent`, Papa with every column's lag
and extra free, 6 + 49 + 49 = 104 unknowns, central differences of the
physical command at h = 1e-7 s): 88 columns with cosine 1.000000000, norm
ratio 1.000000000 and no tick off by more than 1e-6 of the column's
largest entry; the worst relative gap over them 1.6e-9 (the differences'
own rounding). The other 16 columns move nothing: the lags and extras of
nodes 0-5, whose areas are equal in all four of Papa's targets (so
rows[upper] - rows[lower] = 0 exactly), and of the fold columns 46 and 47
(section 9.3); their analytic tangent is exactly zero and their central
difference is the rounding of (1 - s) a + s a against a, below 4 eps
max|command| / h = 6.8e-8. No difference straddled a clip crossing or a
hold's edge at this h (a crossing lies within +-1e-7 s of the start with
probability ~2e-4 a tick); where one does, the two sides' slopes differ
(section 7's reading).

## 7. Where the linearization is honest (MEASURED, `adjoint.py sweep`)

Along each event time at the hinted start (-12.7 dB), the forward
difference |r(theta + d) - r(theta)| / d and the central difference over
+-d, in the residual's units (scaled by the true flow's swing), and the
cosine of each central difference with the one at the scale before:

| d (s) | P- forward / central | cosine | AA- forward / central | cosine | AH+ forward / central | cosine |
|---|---|---|---|---|---|---|
| 1e-7 | 7.77e3 / 7.77e3 | | 4.94e3 / 4.94e3 | | 4.51e4 / 4.51e4 | |
| 3e-7 | 7.77e3 / 7.77e3 | 1.0000 | 4.94e3 / 4.94e3 | 1.0000 | 4.51e4 / 4.51e4 | 1.0000 |
| 1e-6 | 7.78e3 / 1.37e5 | 0.165 | 4.94e3 / 4.94e3 | 1.0000 | 4.51e4 / 4.51e4 | 1.0000 |
| 3e-6 | 7.79e3 / 1.83e5 | 0.634 | 1.07e5 / 5.49e4 | 0.563 | 4.51e4 / 3.88e4 | 0.979 |
| 1e-5 (one step) | 5.26e4 / 1.03e5 | 0.635 | 6.57e4 / 3.39e4 | 0.635 | 4.56e4 / 7.38e4 | 0.509 |
| 3e-5 | 2.05e4 / 5.59e4 | 0.833 | 3.69e4 / 2.03e4 | 0.817 | 5.89e4 / 5.94e4 | 0.815 |
| 1e-4 | 1.41e4 / 2.94e4 | 0.785 | 1.84e4 / 1.07e4 | 0.863 | 4.73e4 / 3.88e4 | 0.814 |
| 3e-4 | 9.44e3 / 1.61e4 | 0.711 | 9.76e3 / 6.30e3 | 0.878 | 3.74e4 / 2.32e4 | 0.753 |
| 1e-3 | 7.60e3 / 8.93e3 | 0.853 | 5.02e3 / 4.40e3 | 0.826 | 2.78e4 / 1.93e4 | 0.703 |
| 3e-3 | 9.53e3 / 5.70e3 | 0.513 | 2.79e3 / 4.02e3 | 0.864 | 1.81e4 / 1.42e4 | 0.837 |
| 1e-2 | 5.62e3 / 2.70e3 | 0.434 | 2.11e3 / 2.68e3 | 0.740 | 1.28e4 / 9.36e3 | 0.703 |

(AA+, P+ and the second P- follow P-, AA- and AH+ respectively; the full
table is in `/tmp/esf/adjoint_sweep.log`.)

Reading:

- Below about 1 us (3 us for AA-, P+, AH+) the forward and central
  differences agree to 3 decimals and do not move with d: the derivative
  exists there and is what the tangent of sections 3-6 would compute.
- Between 0.3 us and 3 us every column meets a jump: the central
  difference's norm leaps 15-30x (P-: 7.8e3 -> 1.37e5 at d = 1e-6: a jump
  of about 0.14 in the scaled residual's norm within +-1 us), and from
  there the difference norms fall with d roughly as d^-0.5, not d^-1: many
  jumps of either sign accumulate. The smooth slope times d equals the
  jumps' contribution only at d ~ 1 ms.
- The cosine between consecutive secant scales is 0.5-0.88 from one step
  to 10 ms: each rung of the puppet's descent (SECANTS 10 ms, 1 ms, 0.1 ms,
  one step) measures a direction that differs by 30-60 degrees from the
  next rung's, and all of them from the pointwise derivative.
- HYPOTHESIS for the jump's source (not isolated here; NOTES_gradients.md
  section 6 isolated the same kinds on the keyframe family): the ear's
  whole-step closure timing (PITCHFB 0.08 per closure, jumps ~4 us apart
  along a time parameter) and a hidden window's edge moving by a step (a
  residual row changing definition). The seal crossing itself is
  continuous under SEALBLEND (a ramp across one engine step of event
  motion) and would show as a slope change, not as d^-1 growth.

Reading (with steps 8, 9, 25, 42 and NOTES_gradients.md section 6): the
equation error is C0 in an event time with its slope changing wherever the
seal crossing moves through an engine step (every 10.4 us of event time:
the SEALBLEND ramp), with jumps wherever a fold closure instant crosses a
step under the ear's whole-step timing (~4 us apart near closures), and
with the voicing onset's non-monotone forced response at ~3 us (physics).
A derivative at a point predicts the residual within about a microsecond of
event time. The puppet's finest secant is one engine step (10.4 us), its
coarser ones 0.1-10 ms: every rung is a secant over this structure, and
step 44's descent (-90.8 dB on Papa) is a descent on those secants. An
exact tangent would differentiate a different, rougher object than the one
the puppet descends on; it earns its place for the spectrum and precision
telemetry at a converged point (Phase C), or once the engine's own
discreteness is timed within the step (the ear's closure, the window edges:
NOTES_gradients.md section 6's proposals, the owner's call).

## 8. Representations compared

Per column of the Jacobian, with Gamma the geometry (0.044 s), P the peel
(0.10 s with `_implied_run`'s bookkeeping), E one engine render (0.175 s
under this load), and s_k the share of the word before column k's first
changed tick:

| representation | cost of column k | exactness | status |
|---|---|---|---|
| central differences, each residual from nothing (today) | 2 (Gamma + P + 3 E) | the function itself | MEASURED, section 1 |
| the same with the three renders as one run (`SharedFieldFlow`) | 2 (Gamma + P + 1.78 E) exact; 2 (Gamma + P + E) canonical | bit-identical (exact mode) / rounding (canonical) | MEASURED, section 9 |
| ... and each residual resumed from the last (`ResidualCache`) | 2 (1 - s_k)(Gamma + P + E) + the record's overheads | bit-identical | MEASURED, section 9 |
| forward tangent, blocks of N lanes (`ad/`, extended to seals, windows, fold targets) | ~0.12 of a central column (NOTES_gradients.md's ratio), + the windows' templated tube | the derivative of the piece the trajectory is on | CANDIDATE, not built: ~1,500 lines of templated C++ (a `Tube` with seal, blend, local linings, static admittance; `SustainRun::control`; `VoiceBody`'s forcing path) and its bit-identity gate |
| adjoint, the full Jacobian | one backward pass per residual row: m = 48,960 passes | - | REJECTED (as NOTES_gradients.md section 2) |
| adjoint, the gradient J^T r | ~3-5 residuals for all p at once | the derivative of the piece | CANDIDATE for gradient methods; LM wants J |
| Gauss-Newton by CG, matrix-free (J v tangent, J^T u adjoint) | ~2 passes a CG iteration, 20-60 iterations | the derivative of the piece | CANDIDATE only above p ~ 500 (the voice problem has ~130) |

## 9. What was built: the shared, checkpointed residual (MEASURED)

`analysis/adjoint.py`:

- `RenderCache`: one `engine.Run` a process, snapshots every 960 steps; a
  render whose settings agree with the last render's before step k
  (`Settings.divergence`, bit for bit over construction, ticks and forced
  waves) resumes from the last snapshot at or before k, or continues the
  live run when everything computed so far agrees; a render may stop at the
  step its caller reads to.
- `SharedFieldFlow`: `pyramid_word.FieldFlow`'s interface over the cache; a
  closure window's render stops at the window's end (`observer._implied_run`
  reads no further), the final render runs to the end.
- `ChainedGeometry`: `moving_geometry` into one buffer, the recurrence's
  state kept every 960 steps, the tail from the first changed tick recomputed.
- `ResidualCache`: `pyramid_word.equation_error` with the last evaluation
  as the base of the next: `observer.implied_resume` for the peel and the
  windows, `ChainedGeometry`, `RenderCache`; the first changed step read
  from the physical command (every column: areas reach the geometry and the
  peel, folds and port reach the engine).
- `jacobian_cached` / `worker_columns_cached`: the columns ordered by
  decreasing first changed tick (each then resumes at its own change), dealt
  to the pool's workers, each with its own caches.
- `command_tangent`, `physical_tangent`: section 6.

### 9.1 Residuals, bit for bit (`adjoint.py check`)

51 commands in the order a Jacobian's worker sees them: the hinted start,
then every event moved by +- each secant (10 ms, 1 ms, 0.1 ms, one step),
then all events moved by +-0.4 ms; each evaluated by `PuppetProblem.residual`
(plain, from nothing), by a `ResidualCache` in exact mode and by one in
canonical mode, the two caches chained from their previous point.

- Exact mode: **51 of 51 identical to the bit** (`tobytes`).
- Canonical mode: relative difference 7e-13 to 2.0e-12 in every residual
  (worst 2.01e-12; worst dB difference 2.2e-12 dB): the rounding of
  forcing the folds with their own recorded wave against letting them run
  free, carried by the oscillator over the rest of the word.
- Engine steps run per residual (T = 51,840): plain 3 T; exact 89,100 to
  93,200 = 1.72-1.80 T; canonical 17,280 to 43,200 = 0.33-0.83 T (the
  live run continued from the last snapshot before the first changed step:
  a column of the last transition (tick 368) runs the last third only).
- CPU totals over the 51: plain 30.4 s, exact 18.9 s (1.61x), canonical
  10.8 s (2.82x). Per residual: plain 0.52-0.72 s, exact 0.32-0.48 s,
  canonical 0.13-0.30 s (load 24-31).
- First changed steps (from the command): the first transition's columns
  9,408-10,368 (tick 98-108 of 540: 19%), the second's 19,488-20,448 (38%),
  the third's 34,368-35,328 (66%); the geometry resumes from the 960-step
  boundary below (9,600 / 19,200-20,160 / 33,600-34,560); the peel from its
  512-step chunk boundary (`implied_resume`).

The seal-less path (`peel_moving_record` / `peel_moving_resume`, one
render a residual, no windows), Alfa `AE,L@alveolar=0.25,F@lips~0.15!,AH`
(44 nodes, jet 0.4, 51,840 steps), the same 51 points: **51 of 51
identical** in exact mode; the canonical form is the exact one there (no
windows: relative difference 0); CPU totals plain 19.8 s, cached 11.9 s
(1.66x): the prefix resume alone (geometry, peel with its jet, the one
render), 0.29 s against 0.39-0.44 s a residual.

Why exact mode keeps 1.78 renders instead of one (OBSERVED, the negative
result of this piece): `observer._implied_run` calls `free_flow` for
closure k with waves that are finite inside the earlier windows, 0.5 (p -
U rho c / A_0) of the recorded forward run there (the folds' own returning
wave: forcing them with it reproduces the free run to rounding, not to the
bit), and NaN from window k's start; `pyramid_word.equation_error`'s final
call has NaN inside every window. The final render therefore diverges from
the last window render at the first window's start, which for Papa (a word
beginning with the sealed /p/) is step 0: the final render always runs
whole, and only the window renders are shared. The canonical form (every
earlier window NaN in every call: the folds free inside a hidden stretch,
which is how the stretch itself was computed) makes the three calls one
run; its cost in exactness is the 2e-12 above. Not adopted into
`pyramid_word.py` here (not this worker's file): CANDIDATE for the main
session, with the owner's rule in mind that this is a representation of the
same equation error, not a constant.

### 9.2 The Jacobian (`adjoint.py columns`, POOL = 6 workers, one column a worker)

At the hinted start, each secant, three Jacobians a secant after a warming
one: plain (`fields_word.worker_columns`'s way: two whole residuals a
column), cached with empty caches (each worker's first residual whole, its
second resumed at the column's change), cached again (the caches warm from
the previous Jacobian: the steady state of an LM whose base did not move
before the column's change). Wall seconds, 6 workers under load 24-31:

| secant | first changed tick (P-, AA+, AA-, P+, P-, AH+) | plain | exact, empty | exact, warm | canonical, empty | canonical, warm |
|---|---|---|---|---|---|---|
| 10 ms | 98 108 203 213 358 368 | 2.37 / 2.56 | 1.71 | 1.55 | 1.31 | 1.10 |
| 1 ms | 107 108 212 213 367 368 | 2.25 / 2.59 | 1.57 | 1.34 | 1.26 | 1.07 |
| 0.1 ms | 108 108 213 213 368 368 | 2.62 / 2.49 | 1.85 | 1.75 | 1.25 | 0.95 |
| one step | 108 108 213 213 368 368 | 2.41 / 2.44 | 1.89 | 1.77 | 1.12 | 0.96 |

(plain measured in both runs: exact's / canonical's.) Every exact column
is identical to the plain one (cosine 1.000000, norm ratio 1.000000, max
|diff| 0, all 4 secants, both cached Jacobians). Every canonical column has
cosine 1.000000 and norm ratio 1.000000 with max |diff| 3.5e-10 to 7.3e-10
(10 ms), 4.1e-9 to 9.4e-9 (1 ms), 4.0e-8 to 7.3e-8 (0.1 ms), 3.7e-7 to
8.5e-7 (one step), against column norms of 2.7e3 to 6.2e4: the residual's
2e-12 divided by twice the secant, as it must be.

The same six columns in one process at the one-step secant (the plain
Jacobian's wall is bounded by one column with 6 workers; sequentially the
sharing across columns shows): plain 7.25 s CPU, exact 4.67 s (1.55x,
identical), canonical 2.61 s (2.50x).

Reading: with one column a worker, the cache's gain is inside the column
(the window renders shared, the behind point resumed from the ahead point's
run at the column's change); the steady-state wall time of a Papa Jacobian
falls 2.4 -> 1.6 s exact, -> 1.0 s canonical. In a real LM the base moves
between Jacobians by the accepted step, which changes every event's
command: the warm columns then resume at min(the column's change, the
step's earliest change) rather than at their own; the "empty" column is the
bound from below, the "warm" one from above.

### 9.3 The voice problem's pattern (`adjoint.py lags`: every column's lag free, 12 columns in one process, the finest secant)

The voice Jacobian hands each of POOL workers a run of ~12 of its ~130
columns, every residual in one process. Papa with all 49 lags free (6 + 49
unknowns), the 12 columns with the latest first changed ticks (lag46 and
lag47 at 540: they move nothing, see below; P-, AH+ at 368; AA-, P+ at 213;
lag44 at 125; lag14, lag29, lag42 at 113; lag7, lag8 at 112), one process,
CPU seconds:

| | plain | cached, caches empty | cached, warm |
|---|---|---|---|
| exact | 14.55 | 8.35 (1.74x) | 8.35 |
| canonical | 14.34 | 4.87 (2.94x) | 4.75 (3.02x) |

Exact: all 12 columns identical. Canonical: cosine 1.000000, ratio
1.000000, max |diff| 3.7e-7 to 8.9e-7 (as 9.2). A lag column's first
changed tick is the first transition's (tick 112-125: 21-23% of the word),
so its prefix saving is small; the gain is the renders'.

OBSERVED on the way (the family, not this module): `fields.read_postures`
reads every column past the areas (the port, the three fold columns)
through the last column's clock, so in the puppet the lags and extras of
the fold columns 46 and 47 move nothing and only column 48's do. Their
columns are exactly zero in every Jacobian (the plain one too: `spectrum`
names them unheard); `voice.py` freezes columns no transition moves, which
covers this. A detail for the main session, not changed here.

Estimate for the voice Jacobian (HYPOTHESIS from these ratios, not
measured on `voice.py`): its ~1,000 residuals are dominated by traverse
columns whose change is at each word's first transition, so the gain is
the renders' share: about 1.7x in exact mode, about 3x canonical, at the
same wall-time pattern (12 columns a worker, each worker's caches warm
from its own run).

## 10. Decisions, negative results and open edges

DECIDED (this worker, within its scope): the shared, checkpointed residual
in exact mode is the representation to use for the puppet's and the voice's
finite-difference Jacobians: the same numbers to the bit, 1.5-1.75x less
CPU (MEASURED, sections 9.1-9.3), no engine change, no constant. Reversal:
delete the module; nothing else depends on it.

CANDIDATE (the main session's call): the canonical form (one engine run a
residual, 2.5-3x; differs by 2e-12 relative). It changes
`pyramid_word.equation_error`'s numbers by rounding; the owner's known
answers (step 44's dB figures) would move in their last digits only, but
they would move.

REJECTED:

- The adjoint for the full Jacobian: m = 48,960 backward passes (section 8).
- "The waveguide is linear in the wave variables": not with the kinetic
  series element on (section 2); the tangent and adjoint are those of the
  linearized element, which is still reciprocal.
- The exact tangent as the puppet's Jacobian at its present rungs: the
  derivative is honest below ~1 us of event time (section 7); the puppet's
  rungs are secants over jumps at 1-10 us, and its descent (step 44) is a
  descent on those secants. Building the puppet-family tangent (seals,
  windows, fold targets: ~1,500 lines of templated C++ and a bit-identity
  gate) was therefore not started here; it is the right tool once the
  engine's whole-step events are timed within the step (the owner's call,
  NOTES_gradients.md section 6) and for the Phase C spectrum at a converged
  point.
- Fewer workers with more columns each, to share prefixes across columns of
  one word: with 6 columns the wall time is bounded by one column either
  way; the sharing across columns pays only where a worker holds many
  (the voice problem).

Negative results kept:

- Exact mode shares 1.22 of the 3 renders, not 2, because of the
  finite-wave semantics inside earlier windows (section 9.1).
- A traverse (lag, extra) column's prefix is the first transition's: 20% of
  Papa; the checkpointing buys little there and the renders' sharing is
  the gain.
- Timing under load: wall and CPU swing with the core (the Mini at load
  19-31); every speedup above is a ratio within one interleaved run.

OPEN:

- Verify on a nasal word (the branched replicas have no cache path here:
  `ResidualCache` handles `word.seal > 0` with `implied_record` and the
  seal-less words with `peel_moving_record`; a velum goes through
  `pyramid_word.equation_error`'s branched path only if it exists there,
  which it does not today).
- The jump's source in section 7 on this family, stage by stage (schedule,
  geometry, peel, windows, folds), as NOTES_gradients.md section 6 did.
- `implied_resume` with a changed window count (a column that opens or
  closes a seal crossing): the record's segments are matched by index and
  start; the gate in NOTES_kernels.md passed its time variants; here no
  point changed the number of windows. Exact mode's bit identity would
  catch a failure; `check` should be run on any new word before trusting
  the cache on it.
- The canonical form's difference at a converged point (-90 dB), where the
  residual is 1e-4 of the start's: the 2e-12 relative is of the residual
  itself, so the dB figure would move by about the same 2e-12 dB, but this
  was measured at the start only.

## 11. House-style review (`planning/PROGRAMMING_HOUSE_STYLE.md`)

Scope reviewed: `analysis/adjoint.py` whole (the only first-party source
written here; no C++ was written or changed; no file of the main session's
was edited). Reviewed for: explicit types (every function's parameters and
return typed; `Callable` for the timed function, `multiprocessing.pool.Pool`
for the pools), named executable behavior (no lambdas, comprehensions or
generator expressions: the loops are written out), retained state named
(`RenderCache`'s run, snapshots, buffers and counters; `ChainedGeometry`'s
buffer and state records; `ResidualCache`'s records; the pool's per-process
`_problem` and `_cache`, `fields_word`'s pattern), ownership and lifetime
(`engine.Run` handles freed by their `__del__` when a snapshot list is
replaced; the geometry buffer owned by `ChainedGeometry` and returned as a
borrow valid until its next call, stated in its docstring), operation order
(the cache's fields published after the run advanced; the first changed
step read before any stage resumes), initialization (every field set in
`__init__`), conversions (`int(...)` at the step boundaries, `float` for
viscous), failure states (`RenderCache` raises if the engine does not
advance; a non-finite render is returned as such and `residual` walls it at
1e3 as `PuppetProblem.residual` does), repeated-loop storage (the geometry
buffer and the engine run reused across evaluations; the RAW buffer is
allocated a render, 830 KB, so the returned prefix and the kept array do not
alias).

Remaining deviations, recorded: `physical_tangent` uses `numpy.where`,
which evaluates both alternatives (the values are finite, no NaN semantics
are at stake; a numpy select is the idiom and a loop over 540 x 49 would
not be clearer); `SharedFieldFlow.__call__` builds an `engine.Settings` (a
C-side copy of the field and waves) per call, the engine's own interface;
`main_*` measurement commands print with `%`-formatting in long statements.
Legacy code (`observer.py`, `reverse_core.py`, `engine.py`, `pyramid_word.py`,
`puppet.py`) was read, not reviewed for style, and is not claimed compliant.
