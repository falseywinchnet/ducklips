> **OWNER RULE (2026-10-20): NO CALIBRATION, EVER.** No fine-tuning, fitting, per-word patches, per-speaker profile fits, or tuned constants. Progress comes only from structure.

# Pyramid control fields (2026-10-20)

Direction (GIVEN, owner 2026-10-20): a lightweight version of the Fourier
feature pyramid (Zhao et al., arXiv 2605.24278) without its MLP decoder,
solved by exact least squares.

## Step 1: field mode in the engine. MEASURED.

- AREATICKS / VELUMTICKS / FOLDTICKS: one row per 1 ms control tick (the
  engine reads its controls only at ticks), the last row held. They replace
  AREASEQ / VELUMSEQ / FOLDKEYS when set. In-process: `engine.render(...,
  field={"areas", "velum", "folds"})`, and `posture_recover.render(field=)`.
- Defaults are byte-identical: Papa and India keyframed truth renders,
  old library against new, compared with `cmp`.
- Keyframes sampled at the ticks and fed as a field: −234.0 dB (Papa) and
  −238.4 dB (India) against the keyframed render. Not bit-identical: the
  C++ fuses multiply-adds that the Python sampler does not.
- `analysis/pyramid.py`: time levels coarse to fine, each twice the last,
  cosine (DCT-II) interpolation (no seam at the word's ends); the finest
  level is the tick grid, so the band a word needs is measured, not chosen.
  Space is not a pyramid: a seal is one node wide.

## Step 2a: representation of the keyframed truth. MEASURED.

Coarse-to-fine exact least squares of the true field, first L levels:

| L (finest cells) | Papa sound | Papa seal crossings | India sound | India seal crossings |
|---|---|---|---|---|
| 5 (9 / 11) | 0.1 dB | 4 of 3 | −2.5 dB | 4 of 4, off by 7 to 64 ms |
| 8 (68 / 87) | −6.8 dB | 11 of 3 | −6.7 dB | 6 of 4 |
| 9 (135 / 173) | −10.9 dB | 3 of 3, off by about 0.7 ms | −11.2 dB | 4 of 4, off by 0.7 to 1.3 ms |
| 10 (270 / 345) | −19.0 dB | 3 of 3, off by 0.1 ms | −23.8 dB | 4 of 4, off by 0.1 to 0.5 ms |
| 11 (all ticks) | −230.4 dB | exact | −240.3 dB | exact |

- The keyframed truth is not bandlimited. Its piecewise-linear corners and
  fast release ramps need the full tick grid.
- Truncated levels ring around the seals (Gibbs): spurious crossings, and
  negative areas down to −1.06 cm².
- HYPOTHESIS: this is a fault of the truth, not of the pyramid. A physical
  command is smooth and aims past contact (a lip closure presses beyond
  touching; the contact clips the area at zero). With area = max(command,
  0), the overshoot is compression, not error, and the event is the moment
  the smooth command crosses the seal. Test next: a bandlimited truth with
  targets past contact, then the inverse.

## Step 2b: smooth physical truth. MEASURED.

- REJECTED, held targets read through one pyramid level as the truth. The
  ringing sealed /aa iy/'s node 15 mid-vowel, and the truth's own equation
  error was −18.8 dB (−32.4 dB for Papa).
- The truth is now targets held from arrival to leaving, with half-cosine
  moves (`pyramid_word.command_of`). A stop's node aims past contact (the
  opening vowel's area, negated). Contact clips areas at STOP_AREA and
  ports and fold gaps at 0 (`physical`).
- Equation error at the truth: /aa ah/ −245.7 dB, /aa iy/ −247.4 dB, Papa
  −247.7 dB. Papa's seal crossings: 125.05, 254.93 and 384.06 ms.
- Pyramid levels for Papa's truth (equation error): L7 (34 cells) −28.1 dB,
  L8 −40.4, L9 −59.9, L10 −83.8, L11 exact. Seal crossings are right to
  0.01 ms from L7.

## Step 3: inverse attempts. MEASURED.

- REJECTED, the radiated flow against the true sound as the residual: flat
  at 0.0 dB at every level; the voice's cycle phase swamps the tube.
- REJECTED, plain Gauss-Newton from the start projected onto level 1. The
  projection loses the stops (no seal), and steps along the null space sent
  areas to 4e5 cm².
- Pyramid corrections on top of the full hint start, with
  Levenberg-Marquardt (Nielsen damping), areas unknown and folds true:
  −3.6 → −5.6 dB through level 7, then stalled, with the second release
  still 7 ms late. 88.5% of the start's error is that release (370–420 ms).
- Scan of the second release at the true tubes: inside the closure window,
  flow comparison is a plateau on the early side (−21 dB from −2 to −10 ms)
  and a cliff on the late side. Radiated volume (the running sum over the
  window) funnels from both sides (−60 at −10 ms, −37 at +10 ms, exact at
  0). Its weight against the open part is not principled yet.
- The seal's flow read from the sound alone (`seal_flow`: the sound peeled
  through the front of the tube only, to the seal node). At the true
  geometry it is about 2e-4 through Papa's closure and 1487 at the release:
  about 140 dB of contrast.
- Observability: the engine reads controls at 1 ms ticks, so an event is
  observable only to the tick.
- `read_events` (largest tick-to-tick rise or fall, no threshold), with the
  hint geometry:
  - Lip releases come out right: Papa 126 and 385 (true first open ticks
    126, 385); /aa b aa/ 256 (true 256).
  - Closures come out wrong: Papa 257 (true 255), /aa b aa/ 134 (true 125),
    Echo 127 (true 125). A closure fades as the area shrinks; it does not
    jump.
  - Echo's velar release comes out at 214 against 249. The front cavity is
    17 nodes long, and the hint's areas there are wrong.
- HYPOTHESIS, next: the front cavity's areas are fixed by the closure
  itself. Nothing passes the seal, so requiring zero seal flow over the
  closure is an equation error for the front of the tube alone. That gives
  exact least squares on the front areas, and then the release and the
  closure read through the right geometry.

## Step 3b: front of the tube from the closure. MEASURED.

- REJECTED as built: LM on Echo's front areas (nodes 27–45, in front of the
  velar seal at node 26), minimising the seal flow over the start's sealed
  ticks. The seal flow went from 2.2 to −0.7 dB and the front area error
  from 1.59 to 9.8 cm² rms. With the true sealed ticks (diagnostic only):
  1.5 to 0.1 dB, 9.7 cm². Not identifiable this way.
- Cause: even at the true geometry the peeled seal flow is not zero just
  after a closure.
  - Echo: 156 over the first 18 ticks (the equation error's scale is 116),
    then 1.4e-2 through ticks 144–244, 9e-2 just before the release, 1341
    after it.
  - Papa: 185 over the first 18 ticks, then 2e-4, then 1487 after the
    release.
- The truncated peel does not model the seal: the engine's SEAL rule makes
  the node a closed end, while the peel passes a 1e-4 cm² node. The early
  transient is that mismatch.
- HYPOTHESIS, next: a seal-aware front peel (the seal node as a closed
  end, as the engine steps it) makes the seal flow exactly zero through the
  whole closure at the truth. Then zero seal flow fixes the front areas,
  and the reader gets closures and velar releases through the right
  geometry.

## Step 4: events read off the sound. MEASURED.

- What must be zero at a seal is the flow through the link in front of the
  seal node, as the engine steps that link. It is not the junction inflow:
  a sealed link reflects one step late and hop-smoothed, so the junction's
  inflow y(2a − p) is not zero.
- `link_flow`: peel the sound through the front of the tube to node + 1's
  pressure p and inflow; a = (inflow / y + p) / 2; s = the leaving wave
  p − a hop-filtered; flow(t) = y (a(t+1) − s(t)).
  - For a lip seal the front is one node, and the one-node peel paired the
    end node with itself. That case uses the open-end stage
    (`observer._states_open_end`) on nodes 44–45.
  - On Echo's 19-node front, the junction inflow so peeled matched the
    forward replica's to 1e-12.
- At the truth, the flow through the link is zero (≤ 2e-12) exactly over
  the engine's sealed steps: Papa 0–126.896 and 258.083–385.938 ms;
  /aa b aa/ 128.302–257.000; Echo 127.542–251.094.
- The engine seals later than the command crosses. The area glides over 16
  steps, so Papa's link seals at tick 260 against a command crossing at
  254.93.
- Reading (`read_events`, then `read_steps` to the step) holds the seal
  node sealed over the whole word, so the candidate's own opening does not
  enter the reading. From the hint start:
  - Releases are exact to the step: Papa 12182 and 37050, /aa b aa/ 24672,
    /aa t aa/ 24650.
  - Closures read 31–58 steps (0.3–0.6 ms) late, because a closure fades
    instead of jumping.
  - Echo's velar release still reads at 214 (true 251): the 17-node front
    is the hint's, and it starts moving at 213.
- Echo's front solved from zero link flow over the closure: area error 1.59
  → 0.30 cm² rms at level 1 (17 constants), link flow −122 dB. Finer
  levels add nothing, and the release then reads late (258): the front's
  motion during the transition (213–251) is not constrained.
- Nasals: the sound includes the nostrils, so reading through the mouth's
  front alone is wrong. OPEN.
- Oscar's smooth truth blew up at 0.257 s (the jet). OPEN.
- Release-shift scan at the true tubes: the forward-sound closure residual
  funnels on both sides (Papa −248 dB at 0, −22 at ±1 ms, −9 / −4 at
  ±10 ms). A link-flow closure residual is −17.7 dB even at the truth: the
  hidden window opens at 0.1 cm², but the link seals only at 0.01. Not
  used.

## Step 5: inverse from read events. MEASURED.

- Start = hint tubes with each transition shifted so the start's sealed
  steps begin and end on the steps read (`align_events`). Start −7.9 dB
  (was −3.6). Releases exact.
- Event times are quantised to the step, so the residual is a staircase in
  them. With the edge on a step boundary, a 1e-9 move jumped Papa's
  residual by 0.9 of its scale at the release steps; the Jacobian there was
  meaningless (|dr| ≈ 1.09 for every h from 1e-8 to 1e-3). Fixes:
  - align to the middle of the band of shifts that gives the read step;
  - hold the seal nodes out of the unknowns. A seal's step depends only on
    its own node, so the events stay where they were read.
- LM from there (areas unknown, folds true): −8.0 → −10.2 dB through level
  5, then stalled. The area error rose 0.97 → 1.10 cm² rms while the cost
  fell. The equation error can be lowered by wrong areas from this start,
  so the basin around the truth is narrow, or some area combinations trade
  off against each other. OPEN.

## Step 6: null space, basin and representation. MEASURED.

- Jacobian of the equation error at Papa's truth, area pyramid corrections
  at level 5 (900 unknowns, seal node held):
  - Singular values from 163 down to 3.8e-6 of the largest; the median is
    1.6e-4. Ill-conditioned (about 4e7), with no exact null space.
  - The weakest directions are area wiggles at nodes 1–5 and 35–40, spread
    over the word.
  - Node 0 (the glottis) is seen 4× more strongly than any other node.
  - Area rms is therefore not the measure: many tubes sound nearly alike.
- Basin, with offsets projected onto level 5's span and the truth's fine
  structure kept in the base:
  - From 30% of the hint's offset: −20.0 → −54.8 dB in 6 LM steps.
  - From the full offset: −11.3 → −28.7 dB, still descending. The area rms
    stays at 0.87 cm² (sound-alike tubes).
  - Wide basin.
- REJECTED: the field as the pure sum of its levels from level 1 (the hint
  projected, the seal nodes held at full resolution). It started at
  +0.9 dB (a tube nearly still over the whole word) and reached only
  +0.2 dB through level 5.
- Reading: the basin test worked because its base kept the truth's fine
  time structure (the transitions' timing and shape). The hint start gets
  exactly that wrong (vowel times off by 7.3 ms), and coarse levels cannot
  remove or rebuild it. Coarse to fine over the whole word is the wrong
  axis for a causal problem.
- HYPOTHESIS, next: solve causally in time, as layer peeling does. The
  equation error up to t depends only on the field up to t, so march short
  windows. Each window's field is solved at full resolution, starting from
  the solved past and the hint ahead. These are small, well-started exact
  least squares problems.

## Step 7: two fields, A = P(ξ, s). MEASURED.

- GIVEN (owner, 2026-10-21): "we need a time field", then "we need two
  fields. space and time." Area against fixed clock time made timing into
  fine area detail; a 7 ms shift needed every fine coefficient.
- `fields.py`: the postures P(u, s) (each phone's tube along the node
  coordinate u, linear in s between targets), the space field ξ(x, t) and
  the time field s(x, t).
  - Known answer: the truth read through ξ = identity and the true s
    reproduces the command to 7e-15. Equation error −250.3 dB (Papa),
    −249.5 (/aa iy/), −247.5 (Echo).
- Representation, time field only, against the area field: for Papa, L7 is
  71 unknowns at −31.6 dB (areas: 3266 at −28); L9 is 274 at −60.3 dB
  (areas: 12604 at −60). About 46× fewer unknowns. Events right to 0.03 ms
  from L7.
- BUG, fixed: the Levenberg-Marquardt damping was λ × (largest squared
  column norm) × each column's norm, which over-damped every step by about
  the largest squared column norm. That was the source of all the
  "accepted, moved nothing" steps, in both inverses. λ is now
  dimensionless (Marquardt).
- REJECTED: the time field carried directly and solved level by level. It
  crossed the lip seal five times where the word has three (progress ran
  backwards), and the coarse fits walked away from the truth.
- The time field is now the running sum of a clipped rate: s = Σ max(r, 0)
  per tick, with r carried by every pyramid level at once (as the paper
  trains). Monotone by construction. The truth's rate is exact at the
  finest level: time error 1e-16, −248.4 dB, 1084 unknowns.
- Inverse, time field only (true postures, ξ = identity), from the hinted
  times (transitions 7.3 ms off), difference step one engine step's worth
  of the field:
  - The equation error went −12.7 → −32.4 dB in 25 steps, then −38.6 dB,
    and stalled there.
  - Lip events within 0.004 ms of the truth.
  - Timing error (field clamped to the word's phones) 0.020 phones rms. By
    stretch: up to 0.071 inside the sealed voiceless closures, which is
    unobservable (correctly unrecovered); ≤ 0.009 in open stretches;
    ≤ 0.002 in vowel holds; 0.0000 after 420 ms.
  - Round trip (the recovered sound against the true sound): 33.7 dB,
    from 6.8 dB at the hinted start. By stretch: 46 / 36 / 31 / 26 dB
    across the first /aa/, 37–42 dB across /ah/, and 9.8 dB at 340–380 ms
    (inside the closure, where the true sound is the front cavity's faint
    ringing).
  - The stall: 69% of the residual is just after the second release
    (395 ms), where /ah/'s voicing restarts after the voiceless /p/. The
    folds' posture is read through the same shared time field, and voicing
    onset after a release is a bifurcation (measured earlier at ns–µs
    sensitivity). A finer difference step (one engine step instead of one
    tick) changed nothing.
- OPEN, next:
  - The folds on their own time field (the folds are their own
    articulator, not the tongue's).
  - Then the postures and the space field.
  - Then the time field per node.

## Step 8: the staircase is in the engine. MEASURED.

- The folds on their own clock (a second clipped-rate time field for the
  fold and port columns): −38.6 → −38.7 dB. Not the stall's cause.
- Barrier: along the straight line from the recovered field to the truth
  the cost rises −38.6 → −22.5 dB before falling to −248. A blend of two
  fields with the same releases crosses the seal somewhere else.
- Events as linear constraints (the time field held at each seal crossing
  to first order): a damped step of norm 16 moved the second release by 92
  engine steps. REJECTED: the clipped rate's kink in the holds is not
  linear.
- Events as knots in the mouth clock (exact by construction: between two
  read crossings, progress runs from one crossing's value to the next in
  proportion to the clipped rate's running sum). With the true postures
  all three Papa events read exactly to the step (closure 24776, releases
  12182 and 37050), and the knots land on the true crossings (125.052,
  254.933, 384.060 ticks). Yet no LM step lowered the cost, with or without
  the exact event check.
- Jacobian against the function at the knotted start: cosines 0.05–0.46
  between J d and the finite change along random directions; the finite
  change is 2–4× larger. Down the gradient: −13.42 → −13.44 dB.
- Uniform time shift at Papa's truth: the change grows exactly in
  proportion to the shift (1110 per phone) up to 1e-5 phones (about
  0.6 µs). Between 1e-5 and 3e-5 it jumps (error −86 → −51 dB). Beyond
  that the jumps dominate: at 1e-3 the change is 3.35, against 1.1 for the
  smooth part.
- So the equation error is a smooth function plus a staircase whose stairs
  are about 2e-5 phones apart and far larger than the smooth part. The
  stairs are events quantised to the 96 kHz engine step: a link sealing
  when its glided area crosses the seal, the folds' contact starting and
  ending. Any parameter that moves an event moves it by whole steps.
- HYPOTHESIS, engine refinement (the owner's category: the geometry or
  engine made to fit speech): contact and sealing located within the step
  (sub-step event time, the step split at the event, or a contact law
  continuous in the gap). The equation error then becomes continuous in
  every timing and area parameter, and the Jacobian means something. Before
  any further inverse work, find which events make the stairs.

## Step 9: the seal blended within the step; the inverse converges. MEASURED.

- Stair source (uniform-shift test at Papa's truth, 1e-5 against 3e-5
  phones): the first release's seal edge moved one step (12182 → 12181);
  the stair was 0.49 of the scale over the 2 ms after it (the smooth part
  0.033), with smaller after-stairs over 35 ms of voicing. The hidden
  windows did not move.
- ENGINE, `SEALBLEND` (`Tube::set_seal_blend`; replica
  `observer.seal_share`, passed to `_forward_steps`):
  - In the step a link's narrower area crosses the seal, the link is open
    for the share φ of the step on the open side (areas linear from the
    step before). Its flow is φ × the open flow, and the waves carry it.
  - No new constant. Off by default: the Papa and India keyframed renders
    are byte-identical against the pre-blend library.
  - With it, the truth's equation error is −248.4 dB (the replica matches).
  - Engine output across the release: unblended, it jumps 59.5 in one
    2e-6-phone increment; blended, it ramps (4–20 per increment) across
    about one engine step. Continuous, and steep.
- The branched (nasal) forward replicas do not have the blend yet. OPEN.
- BUG, fixed: Marquardt's column-norm damping left clipped coefficients
  (column norms near 7e-7) undamped, giving steps of norm 1e6 that raised
  the cost. The step is now measured in each coefficient's quantum
  (`problem.steps`) and damped uniformly, with μ0 = τ·max diag (Madsen,
  Nielsen and Tingleff). The pool's Jacobian columns matched in-process
  columns exactly (cosine 1.0000): the columns were right, the step was
  not.
- The knot values are now the exact level crossings of the postures
  (`fields.level_crossing`). Interpolated from the start between ticks,
  they had left the truth through the knots 1.3e-2 off.
- Papa, time fields only (mouth and folds clocks, true postures, ξ the
  identity), from the event-aligned hinted start, SEALBLEND:
  - −12.3 → −43.0 dB equation error in 33 LM steps (it stalled at −38.6
    before). Releases within 0.002 ms.
  - Mouth clock 0.0172 phones rms. By stretch: ≤ 0.0016 across the first
    /aa/; 0.05–0.09 around and inside the closures; 0.0000 from 440 ms.
  - Round trip 35.7 dB. By stretch: 51.2 / 42.9 dB across /aa/, 40 dB at
    the releases, and about 30–31 dB from 440 ms on, where both clocks are
    exact.
  - Reading: at the end of the word the timing is right but the waveform is
    30 dB off. That is the folds' limit-cycle phase: a self-oscillator has
    a neutral phase direction, so a tiny earlier difference shifts every
    later cycle. HYPOTHESIS: the waveform round trip there is bounded by
    phase, not by the recovered motion.
  - The folds' clock is 0.038 phones rms off, worst inside the closures
    and at the voicing onsets.

## Step 10: postures, space field and clocks together, from the hints. MEASURED.

- `fields_word.py full` (FullProblem):
  - Unknowns: the two clocks, the vowels' tubes (a stop's follows from its
    vowel), and the space field ξ = x + d(x, t) on the tube's 3 coarsest
    levels × 6 time levels.
  - The knots' values follow the postures.
  - Known answer exact: command 7e-15 off, −251.2 dB.
  - Hint start: vowels 1.01 cm² rms off, −6.9 dB; events read off the sound
    within one step.
- REJECTED: the space field unweighted. A still space field is a different
  posture, so the two traded: sound −8.1 → −16.1 dB while the vowel areas
  went 0.93 → 1.56 cm² off, and ξ moved the seal node and its events.
- Gauge: ξ = x + w(t) d m(x), with w = 4f(1 − f) (f the mouth clock's
  fraction: zero at every target, so a posture is the tube at its target)
  and m zero at the seal nodes (a stop's place is the transcript's). Events
  then stayed pinned, but the space field still took the error: |d| 10
  nodes on the first step, for a truth that has none. The vowels did not
  improve (1.01 → 1.14).
- Postures and clocks without the space field, each kind of unknown scaled
  to unit mean diagonal (one damping for both; before, the rates' larger
  diagonals over-damped the postures):
  - −6.9 → −21.3 dB in 10 steps, releases within 0.003 ms.
  - But the timing got worse (0.064 → 0.141 phones rms; the folds' clock
    0.061 → 0.112) and the vowels barely moved (1.01 → 0.94 cm²). The sound
    improves by compensation: wrong postures with compensating timing.
- Reading: this is acoustic-to-articulatory non-uniqueness. Many tubes sound
  alike (the Jacobian at the truth has a condition number of about 4e7),
  and with 46 free areas per vowel, timing and shape trade. The sound alone
  does not pin 46 free areas per posture.
- HYPOTHESIS: the postures belong in a body's coordinates (jaw, tongue
  body and tip, lips, a few muscles), not 46 free areas. A body has few
  degrees of freedom and cannot make most sound-alike tubes, so the
  non-uniqueness shrinks to what anatomy allows. This is the owner's
  puppet: the geometry learned across voices is the prior.

## Step 11: postures in body coordinates; sound-alike valleys. MEASURED.

- `BodyProblem` (`fields_word.py body`): each vowel's posture as the
  articulatory model's 9 parameters (posture_recover.FIELDS: jaw, tx, ty,
  tipx, tipy, ap, pro, lw, drop) through `posture_recover.nodes`.
  `posture_recover.profile` now prints parameters at %.17g (six decimals
  quantised finite differences; table values parse identically).
  - Known answer exact: tubes 0 off, −251.2 dB. 2186 unknowns, 18 of them
    posture.
- From the hints (space field off): −6.9 → −18.6 dB in 12 steps. The first
  step threw the mouth clock 0.064 → 0.221 phones rms and the folds' clock
  a whole phone off. The vowels stayed about 0.88 cm² off.
- Basin tests, true timing, articulatory parameters 20% / 50% of the way to
  the hints:
  - INVALID (my error): knots taken from the start's own crossings. The
    moved postures had moved the last release (+1.95 / +4.98 ms), so the
    knots pinned it wrong.
  - Redone at 20% with knots at the events the sound gives (the truth's
    crossings): −17.4 → −32.8 dB in 12 steps, releases within 0.02 ms.
    Yet the first step threw the mouth clock 0.020 → 0.130 phones rms (it
    recovered to 0.045), the folds' clock went 0.65–0.73 phones off, and
    the vowels went 0.20 → 0.26 cm² off.
- Reading: the sound is matched better and better while the telemetry
  wanders. The solution set is a valley of sound-alike motions, and LM
  moves along its floor. The folds' clock is the clearest case: 0.7 phones
  off at −32.8 dB, because the fold movements inside voiceless closures are
  not heard. The sound alone does not pin all of the articulation, even in
  body coordinates.
- OPEN, the owner's decision: how to choose within the valley.
  (a) A physical least-action principle (articulators move with least
  effort: minimum jerk or rate change per node; the per-node "friction" the
  owner wants learned).
  (b) Report only the observable components, each with its uncertainty,
  and mark the rest unrecovered.

## Step 12: least action, and the agreement gate. MEASURED.

- ACTION=1: after each accepted LM step, a step lowering the clocks'
  action (Σ (Δ rate)², least jerk of progress) within the directions the
  sound does not see (σ² below LM's own damping), halved until the sound is
  no worse.
  - Inside the LM step it was never accepted: the unseen set grew with
    the damping.
  - As a lexicographic second step it was never accepted either: in the
    unseen directions the sound rose by a hair.
  - The known-answer truth is half-cosine transitions with holds, not
    least jerk, so this truth cannot test the principle. A fair test needs
    a truth made by it. OPEN.
- AGREEMENT 0.25 (Nocedal and Wright's η): a step is kept only when its
  actual gain is at least a quarter of its model's prediction. Before, any
  gain kept a step that slid along the valley for 0.1 dB. 20% basin test,
  before → after the gate:
  - sound −32.4 → −38.2 dB;
  - mouth clock 0.055 → 0.036 phones rms;
  - folds' clock 0.16 → 0.032;
  - vowels 0.32 → 0.26 cm² (start 0.20).
  - The first step still pushed the timing to 0.107 before it came back.

## Step 13: which articulators the sound pins. MEASURED.

- Gated body inverse from the hints: stalled at −18.7 dB (timing 0.19
  phones rms off, folds' clock 0.33, vowels 1.0 cm²). Failed.
- Staged, stage 1 (postures from the holds' central halves, clocks held):
  hold sound −22.2 → −25.6 dB; vowels 1.01 → 0.95 cm². The holds do not
  pull the articulators to the truth from the hints. (The old pipeline sat
  near 1 cm² of tube error at −55 dB too; tubes converged only at
  −67 dB (0.026) and −84 dB (0.0006).)
- Per-articulator observability at Papa's truth: the hold sound's change
  for a 1% move of each articulatory parameter
  (`/tmp/esf/body_observability.py`):
  - /aa/: tongue body height +14.2 dB, jaw +4.7, tongue body front-back
    −4.4, tip −8 to −9, lips −15 to −19 (aperture −18.6, protrusion −17.0,
    width −14.9), drop −20.5.
  - /ah/ (its hold window shorter): tx +0.9, ty −1.9, jaw −13.8, tip about
    −27, lips −32 to −42.
  - Singular values from 1 down to 4.0e-8. The two weakest directions:
    - /ah/ lip aperture against lip width (+0.78 / −0.62, 4e-8): the sound
      hears the lip opening's area, not its height and width;
    - /ah/ protrusion (2.2e-7): the lips' length is nearly inaudible.
- Reading: the articulatory model has degrees of freedom the sound cannot
  separate. For telemetry and the motion analysis, the coordinates the
  sound pins are the opening areas and the tongue body and jaw. Lip shape
  is only its area; protrusion is barely heard. HYPOTHESIS: reparameterise
  the lips as opening area and length, and the inverse loses its exact
  nulls.

## Step 14: surveys of the postures. MEASURED.

- Knot fix: knot ticks are now where the unclipped command crosses the
  seal. The clipped area's interpolation between a clipped and a free tick
  put the knots off: the truth's hold sound was at −59.2 dB; now −256.0.
- Hold sound along the straight line from the hints' to the truth's
  articulatory parameters (clocks true): −27.2, −27.8, −27.3, −26.8 at 0
  to 30%, then −28.8, −29.9, −31.4, −33.9, −37.0, −42.0, −48.7, and −256 at
  100%. A plateau, not a barrier: no slope where the hints start.
- FREEZE_FIELDS=lw,pro holds the unheard lip parameters (step 13's nulls).
- Survey scored by the hold equation error (256 Sobol postures over the
  hint table's range per heard parameter): the best was no better than the
  start (/aa/ −28.5 vs −28.5, /ah/ −25.6 vs −24.0). Stage 1 then stalled
  at −25.7 dB with the vowels 0.91 cm² off. The equation error is flat away
  from the truth.
- Survey scored by the hold's harmonic levels up to 5 kHz (as
  tube_recover.SpectrumTarget; 1024 postures): /aa/ 0.576 → 0.284 nepers
  rms, /ah/ 0.579 → 0.218. But the postures went further from the truth
  (1.30 cm²), and LM from there reached hold sound −10.4 dB with the vowels
  2.14 cm² off. Spectrum-alike postures are not the true posture.
- Reading: every measure finds a sound-alike posture; none leads back to
  the truth from the hints. Articulatory telemetry from sound is pinned
  only in the combinations the sound hears. NEXT: map those combinations
  per vowel (the observability table), so the telemetry and the motion
  analysis use only what is pinned.

## Step 15: observability across vowels. MEASURED.

The hold sound's change for a 1% move of each articulatory parameter (dB
against the scale), two vowels a word, at the truth (`body_observability.py`):

| | jaw | tx | ty | tipx | tipy | ap | pro | lw | drop |
|---|---|---|---|---|---|---|---|---|---|
| AA (AA AH) | −0.7 | −5.3 | 9.4 | −13.3 | −13.9 | −19.9 | −22.1 | −16.5 | −24.7 |
| AH | −4.7 | 0.0 | 7.5 | −7.7 | −7.3 | −29.1 | −36.7 | −26.9 | −25.9 |
| IY (IY OW) | −2.3 | 27.2 | 22.7 | −18.7 | −3.3 | −0.4 | −14.5 | −10.3 | −10.0 |
| OW | 16.8 | 22.9 | 28.1 | 15.3 | 16.0 | −19.3 | −12.4 | −27.2 | 9.8 |
| EH (EH IH) | −30.5 | −16.9 | −10.5 | −31.1 | −28.9 | −35.0 | −54.5 | −34.4 | −38.2 |
| IH | −15.8 | −2.6 | −2.1 | −49.9 | −38.3 | −43.7 | −34.5 | −39.5 | none |
| AE (AE UW) | −26.8 | −18.5 | −1.4 | −42.1 | −54.6 | −37.8 | −82.4 | −38.6 | −27.6 |
| UW | −13.6 | −11.7 | −2.5 | −17.1 | −34.7 | −15.0 | −20.3 | −26.5 | −15.0 |

Weakest directions:
- AA AH: ap/lw (7.0e-10), then pro (1.9e-6).
- IY OW: lw/ap/pro (4.7e-8 and 1.0e-7).
- EH IH: IH drop (exactly 0: the parameter does nothing for IH), then
  IH lw/ap/pro (4.8e-9).
- AE UW: AE pro (1.2e-5), then UW lw (7.8e-5).
- Papa's AA AH (step 13): ap/lw 4.0e-8.

- The tongue body (tx, ty) is pinned in every vowel. The jaw and tip
  where the vowel uses them.
- The lips' shape is not pinned: aperture against width is the weakest or
  second-weakest direction in every word measured, and protrusion is nearly
  null. The sound hears the lip opening's area.
- Values differ by position in the word (the hold windows differ in length).
- ENGINE / GEOMETRY CANDIDATE (the owner's refinement list): the lips in
  the articulatory model have two degrees of freedom (aperture, width) where
  the sound hears one (area), and protrusion barely. A lip parameterised as
  opening area (and length) removes the exact nulls.

## Step 16: the heard subspace is small and nonlinear. MEASURED.

- Only jaw, tx, ty free (2 vowels), everything else at the truth, from the
  hints' values, LM on the holds: hold sound −27.5 → −34.1 dB, then
  stalled. Off the truth at the end: /aa/ jaw 5.93 (start 5.78), tx 0.14,
  ty 0.90; /ah/ jaw 4.31 (from 8.38), tx 0.11, ty 0.22. The jaw and the
  tongue body trade.
- Each vowel's own hold Jacobian (9 parameters, 1% units, at the truth):
  - /aa/ singular values 1, 2.0e-2, 4.9e-3, 2.5e-6, then 1e-17 (exact
    zeros).
  - /ah/ 1, 5.0e-2, then exact zeros.
  - The sound of a held vowel sees 2–4 combinations of its articulators,
    and several parameters change the tube not at all at that posture.
- Gauss-Newton on /aa/'s 3 heard coordinates only (the rest held at the
  hints), from the hints: (9.45, −30.72, −9.80) → (6.85, −24.50, −27.16),
  hold sound −39.6 dB, stalled. /ah/ did not move. The map from body to
  sound is strongly nonlinear over the hints' distance.
- CONCLUSION (MEASURED, this engine and articulatory model, known answer):
  - From a start as far off as the hints, the sound does not return the
    articulation.
  - Near the truth it pins about 2–4 combinations per held vowel, and the
    timing where events and transitions are heard.
  - Accurate articulatory telemetry from sound alone is therefore a few
    numbers per posture plus the heard timing. The body's other freedoms
    must come from the body: its mechanics, least action, and the
    cross-voice learning of the owner's plan.

## Step 17: Echo, the velar release. MEASURED.

- Time fields with the true postures and knots from the hinted read:
  Echo's velar release read 31.5 ms early (step 20431 against 24105), so
  the knot pinned it wrong. Start 0.58 phones off.
- Without knots (EVENTS=0), with SEALBLEND: −12.2 → −24.8 dB in 40 steps;
  the release stayed about 2 ms early (−1.99 ms; start −1.51); mouth clock
  0.034 rms. The sound alone does not pull a velar release into place.
- Release read through the hinted timing 20431; through the knot-free solve
  23882 (223 steps early); through the true timing 24105 (exact). The
  velar closure reads 12283 even through the truth (true 12244): its
  reading is biased by 39 steps.
- REJECTED: read-and-solve rounds (`fields_word.py iterate`). A knot at
  the read release (2.3 ms early) distorted the front's motion; the next
  reads went 21199, then 25983, and the field diverged (0.49 phones rms).
  A read through a solved field carries that field's bias back into the
  knot.
- OPEN: the velar release. Its reading needs the front cavity's motion
  (17 nodes) right before the release, and no stable way to get it yet.

## Step 18: the observability map of the timing. MEASURED.

`analysis/observability.py`: for each clock, a one-quantum raised-cosine
nudge (10 ms wide, every 5 ms) of the solved time field, and the change in
the equation error it makes (dB against the scale). Papa, the solved time
fields (step 9's solution, true postures):

- Mouth clock:
  - word-initial closure 5–120 ms: none at all (−6000: sealed lips, nothing
    radiates);
  - first release 125 ms: +17 dB;
  - /aa/ 130–250 ms: −1 to −14;
  - second closure 260–290: none; 295–375: −14 to −17;
  - second release and /ah/ 380–415: +9 to +21;
  - after the last target (425+): none (the field is clamped there).
- Folds' clock:
  - −77 to −85 dB through the word-initial closure;
  - −14 to −31 through /aa/;
  - −52 to −103 through the second closure;
  - −16 to −20 through /ah/;
  - none after the last target.
- Reading: the mouth's timing is pinned sharply at the releases and
  moderately through the vowels, and is unheard where the lips are sealed.
  The folds are heard only while voicing; inside the voiceless closures
  their timing is 80–100 dB down.
- This map is the telemetry's uncertainty: the motion analysis (the
  owner's puppet plan) should take the timing only where the map is
  heard, and treat the rest as chosen (by least action), not measured.

## Step 19: smooth truths of the other words. MEASURED.

- Truth equation error, smooth commands, SEALBLEND:
  - Golf (G AA L=0.25 F~0.15!): −254.3 dB;
  - Alfa (AE L=0.25 F~0.15! AH): −244.4 dB;
  - Echo: −246.8 dB;
  - Oscar (AA S~0.08! K! ER): −29.9 dB (the same without the blend).
- Oscar: the engine's own true render goes non-finite at 0.2573 s, as the
  velar seal (node 28) closes (area 0.0001 at tick 255) while the /s/
  jet's gate is still on. The forward replica matched the engine to 2e-13
  up to there, and the equation error was exact to 250 ms.
- ENGINE CANDIDATE (the owner's refinement list): the explicit jet is
  unstable when a seal closes with the jet gated on. The keyframed truth
  closed fast enough to pass; a smooth closure does not. (Already on the
  open list as the jet's stability margin.)
- Nasal words: the field inverse has no branched (nose) equation error,
  and the old pipeline's march carries tuned constants (WEAK_ORAL 0.7)
  that the no-calibration rule forbids. OPEN: a structural redesign, not
  a port.

## Step 20: the jet bounded by the lungs. MEASURED.

- Cause of Oscar's blow-up (`VoiceBody::reversible_jet`): the jet's speed
  at each link past the narrowest node is |U_link| / A_narrowest. As the
  velar link closes under the /s/ jet's gate, A_narrowest heads toward
  zero while the links past it still carry the /s/ region's flow, so the
  speed runs away (the strength goes as speed³). The closing fade
  (REVJETCLOSE) reaches zero only below 0.01 cm².
- ENGINE, `JETBOUND` (`VoiceBody::set_jet_bound`; replicas:
  `reverse_core._jet_source`, `kernels_jet.jet_source`, threaded through
  `observer._forward_steps`, `_states_link` / `_states_advance`
  (StatesPeel) and `reverse_core._side_link` / `_side_advance`
  (SidePeel); `reverse_core.JET_BOUND`):
  - the speed is bounded by the fastest the lungs can drive air,
    √(2 p_lung / ρ) = 3746 cm/s for 800 Pa;
  - a physical bound from known quantities (no new constant), and it keeps
    the jet reversible (the bound is the same at every link);
  - off by default: the Papa and India renders are byte-identical.
- With JETBOUND=1 and SEALBLEND=1, truth equation error over the whole
  word: Oscar −255.6 dB (before: non-finite at 0.2573 s), Papa −249.8,
  Golf −254.3.
- The nasal march (marching.py) does not yet pass the bound. OPEN.

## Step 21: Golf, the timing inverse. MEASURED.

- Golf (G@velar AA L=0.25 F@lips~0.15!), time fields with the true
  postures, knots from the read events, JETBOUND=1, SEALBLEND=1.
  - The word-initial velar release read off the leading silence: step
    12175, knot at tick 125.332 (exact against the truth's crossing).
  - −25.7 → −38.6 dB equation error in 40 steps; the release exact
    (+0.000 ms).
  - Mouth clock 0.0166 phones rms; folds' clock 0.0242.
  - Round trip 35.8 dB. By stretch (40 ms each):
    - 120–160 ms (the release into /aa/): 26.6 dB;
    - /aa/ and /l/: 32–35 dB;
    - 400–440: 55.1;
    - the voiceless /f/ 440–560: 90.4 / 102.2 / 96.1 dB (the hashed jet
      noise, reproduced exactly).
  - Timing error by stretch: up to 0.102 phones inside the initial velar
    closure 80–120 ms (unheard); ≤ 0.015 through /aa/; 0.0000 in /f/. The
    folds' clock is 0.04–0.06 off through /aa/ and /l/, where voicing is
    heard. NOTE: the per-stretch "sealed" column in roundtrip_fields.py
    reads node 44 (the lips) and does not see Golf's velar seal.
- As with Papa: where the motion is heard the timing comes back to
  0.002–0.015 phones, and the waveform round trip in voiced stretches sits
  near 32–36 dB. The unvoiced /f/ is exact to 90–100 dB, which supports
  the limit-cycle-phase reading of the voiced ceiling.

## Step 22: the voiced waveform round trip measures phase. MEASURED.

- Golf's true timing shifted uniformly, the voiced stretch 160–400 ms:
  | shift (phones; ms in a transition) | voiced round trip | equation error |
  |---|---|---|
  | 0.0005; 0.03 | 14.5 dB | −44.0 dB |
  | 0.001; 0.06 | 10.7 dB | −41.0 dB |
  | 0.002; 0.12 | 7.2 dB | −33.2 dB |
  | 0.005; 0.30 | 3.4 dB | −31.4 dB |
  | 0.01; 0.60 | 2.0 dB | −28.6 dB |
- A 30 µs shift collapses the voiced waveform match to 14.5 dB. The
  recovered timing gives 32–35 dB there: where it matters it is
  effectively finer than 30 µs.
- Conclusion: in voiced speech the waveform round trip measures the folds'
  oscillation phase (which carries every past timing difference forward),
  not the recovered motion. The equation error tracks the motion. The
  unvoiced /f/ (no oscillator) reproduces to 90–102 dB. "Both ways" for
  voiced speech is judged by the equation error and the recovered motion,
  not the raw waveform.

## Step 23: the velar release, more. MEASURED.

- Echo's knot-free solve (−24.8 dB): 90% of the residual at 240–280 ms,
  the velar release (true 249–251 ms); the mouth clock 0.075 phones off
  there.
- The release axis alone is a funnel (the true tubes, K's leaving moved):
  −4 ms −15.3 dB, −2 −19.0, −1 −22.8, −0.25 −26.8, 0 −246.8, +0.25 −26.2,
  +1 −22.8, +3 −15.0. The stall 2 ms early is a joint valley: the rest of
  the timing compensates the release.
- REJECTED: reading with the front held at the stop's posture. It leaks
  from 220 ms (5e6–1.6e7 per tick), when the true front starts toward
  /ow/, and reads 21199 (true 24105).
- The ACTION step's unseen set is now the numerical null space (σ ≤ √ε
  σ_max), where the sound is exactly silent (step 18's −6000 dB stretches).
  NOTE: whether bodies move by least action cannot be tested on a
  synthetic truth; a truth made by the principle tests only the mechanism.
  It needs real voices.

## Step 24: Alfa, roughness without seals. MEASURED.

- Alfa (AE L=0.25 F~0.15! AH; narrowings, no seal, so no events or knots),
  time fields with the true postures, JETBOUND, SEALBLEND: from −19.5 dB,
  no LM step was kept (step 0: "no damping lowered it").
- Uniform shift at Alfa's truth: linear up to 1e-5 phones with a slope of
  2.8e4 per phone (Papa's 1110: 25× steeper); then super-linear (change
  3.41 at 1e-4 against 1.7 expected) and non-monotone (−36.1 dB at 1e-4,
  −39.7 at 3e-4).
- Not the jet: with the jet's gain at 0 the scan is nearly identical
  (−36.2 / −39.8).
- HYPOTHESIS: the folds' contact. With no seals left, the remaining
  step-quantised events are the folds' closures each cycle, and their flow
  kinks at contact (the fold area floored at kMinimumArea). ENGINE
  CANDIDATE: the folds' contact located within the step, as SEALBLEND did
  for the seal. Not yet located; first find which steps change in the
  jump between 3e-5 and 1e-4.

## Step 25: voicing onset is oscillatory, not a stair. MEASURED.

- Alfa's non-proportional change lies entirely in 400–440 ms (56% / 44%):
  the voiceless /f/ into /ah/, voicing restarting. Peak at 427.458 ms,
  among the true flow's first near-zero minima (the folds' first contacts
  as they adduct).
- The folds' contact law is continuous: the wedge's contact length and area
  move smoothly, the soft contact has no corner, and the area floor never
  binds (the chink keeps 0.0002 cm² open). Not step-quantised like the
  seal.
- Fine scan of the onset window's error (rms over 400–440 ms against the
  scale) against a uniform shift: 3e-5 0.0082; 4.5e-5 0.0073; 6e-5 0.0036
  (a dip); 7e-5 0.0094; 8e-5 0.0205; 9e-5 0.0507; 1e-4 0.0557. Continuous
  and non-monotone, varying at the scale of about 5e-5 phones (~3 µs).
- Reading: the folds' forced response near the voicing threshold is
  extremely sensitive to the phase of their first contacts against the
  returning waves. Physics of onset, not an artefact. A start 0.03 phones
  off sees many such ripples. This explains Alfa's stuck start, Papa's
  stall at /ah/'s onset (step 7) and the after-stairs (step 9).
- OPEN: voicing onsets need a measure that does not ride the onset's phase.
  Candidates: treat each onset's first cycles as the folds' own decision
  (a window compared by a phase-free quantity, as volume did for a
  release), or read the onset off the sound as an event (the first
  contact's step), as releases are read.

## Step 26: volume across a voicing onset. MEASURED.

- Alfa's onset window (400–440 ms), uniform shifts of the true timing;
  flow rms against volume rms (the running sum of the difference over the
  window, divided by its length):
  | shift (phones) | flow | volume |
  |---|---|---|
  | −3e-4 | 0.0264 | 1.02e-2 |
  | −1e-4 | 0.0166 | 2.98e-3 |
  | −5e-5 | 0.0067 | 3.37e-4 |
  | −2e-5 | 0.0048 | 3.21e-4 |
  | 0 | 0 | 3.2e-13 |
  | +2e-5 | 0.0067 | 4.23e-4 |
  | +5e-5 | 0.0059 | 6.34e-4 |
  | +1e-4 | 0.0557 | 4.22e-3 |
  | +3e-4 | 0.0344 | 1.35e-2 |
  | +1e-3 | 0.1384 | 1.12e-2 |
  | +3e-3 | 0.0361 | 2.21e-3 |
- Volume is monotone within ±3e-4 phones; the flow measure is not.
  Beyond that, volume is non-monotone too (at 3e-3 a glottal cycle's flow
  nearly cancels in the sum). It widens an onset's funnel about 15×; it is
  not a global measure.
- OPEN: read a voicing onset off the sound as an event (the first glottal
  contact in the implied flow) and knot it, as releases are knotted. That
  needs a map from the onset to the folds' clock, which the folds'
  dynamics, not their posture, decide.

## Step 27: masking voicing transitions; the march carries the jet bound. MEASURED.

- ONSETMASK (`pyramid_word.voicing_mask`): the equation error is zero where
  the folds' rest is strictly between the word's fold targets (abducting or
  adducting) and for one glottal period after (the true flow's period). It
  is fixed during a solve, from the start field: a mask moving with each
  candidate would put its edges' jumps back into the cost. Alfa: 13942 of
  51840 steps masked.
- Alfa with the mask: from −30.3 dB, still no LM step kept.
  - Jacobian against the function along random directions: cosines 0.54,
    −0.43, −0.40.
  - Toward the truth: 0 −30.26, 0.1 −27.70, 0.3 −27.70, 0.5 −33.38,
    0.7 −37.65, 0.9 −34.58, 1.0 −254.51 (non-monotone).
  - 88% of the mismatch at 420–440 ms, just past the mask. After the folds
    reach the voiced posture their oscillation builds over several cycles,
    and the phase sensitivity lasts through the build-up, not one period.
- OPEN: an onset's end measured, not chosen: where the implied flow's cycle
  amplitude stops growing. That still needs a criterion for "stopped".
- `marching.jet_table`'s tuple carries the jet bound as its ninth element
  (reverse_core.JET_BOUND), passed at the march's three jet-source calls.

## Step 28: the onset mask until the first collision. MEASURED.

- ONSETCOLLISION: the mask runs on to the first cycle whose minimum falls
  to the median of the cycle minima over the stretch that follows (the
  true glottal flow's cycles as a known-answer proxy), plus one period.
  Alfa: masked 213.0–428.2 ms (the whole /f/).
- Alfa at its start: −35.6 dB. Toward the truth: −35.63, −34.24 (0.1),
  −34.74, −40.17, −44.03, −47.20, −255.27: nearly monotone. Yet the
  Jacobian's cosines along random directions are 0.83, −0.79, −0.81, and
  92% of the mismatch remains at 420–440 ms.
- True cycle minima: 70.4 at 402.8 ms, −1.6 at 410.5, then −3.0, −3.2,
  −2.7, −3.1 (at 440.9), and steady at −3.1 from there. The folds reach
  the closed phase at about 410 ms but settle only at about 441. The
  mismatch lies in those settling cycles, after the mask ends.
- Reading: the onset's phase sensitivity lasts until the limit cycle has
  settled, not until the first contact. Ending a mask at "settled" needs
  another threshold. Route 1 (masking the onset) stops here.
- OPEN: route 2, reading each onset as an event and knotting it; or a
  folds' model whose forced response near threshold does not turn on
  microseconds.

## Step 29: masks cannot carry onsets. MEASURED.

- Mask extended to the first collision at abductions too (a bug): it ran
  across the whole voiceless /f/ (20662 of 51840 steps). The masked cost
  went −35.6 → −58.4 dB in 10 steps while the timing got worse (mouth
  0.029 → 0.031 phones rms, folds 0.029 → 0.059): the mask hid the timing.
- Fixed (extended only at adductions, where the folds' rest falls toward
  voice): 14673 steps masked, start −34.0 dB, and no LM step kept.
- Conclusion: masking the onset either leaves its settling cycles, rough
  enough to block the solver, or hides the timing. Neither boundary is
  principled. Voicing onset is a structural problem of the folds' model or
  of how its onset is observed, not of the equation error's bookkeeping.
  OPEN, the owner's decision: (a) read each onset off the sound as an event
  (the first glottal pulse in the implied flow) and knot it; (b) a folds'
  model whose forced response near threshold does not depend on
  microseconds of phase.

## Step 30: onsets read as events (route a). MEASURED.

- `fields_word.py onset`: each voicing onset after a voiceless phone is
  read as the first collision in the glottal flow the sound implies through
  the start's tube (`pyramid_word.implied_flow`, `first_collision`); the
  folds' adduction is moved by bisection until the candidate's own rendered
  first collision lands there; then the folds' clock is held and only the
  mouth clock solved, no mask.
- Alfa: onset into /ah/ read at 412.979 ms through the start's tube,
  against 420.594 through the true flow (7.6 ms early). The bisection ran
  to its −20 ms bound (the candidate's first collision never crossed the
  read step inside it). The folds' clock ended 0.57 phones off, the start
  at −11.6 dB, and no step was kept.
- REJECTED as built: an onset read through a wrong tube is biased, as the
  velar release was. Reading events through the start's tube works for
  lip releases (a one-node front) and not for events that depend on much
  of the tube (velar releases, voicing onsets).

## Step 31: outside-in on a velar release. MEASURED.

- `/tmp/esf/outside_in.py`: the cavity in front of Echo's velar seal
  (nodes 27–45) gets its own clock (one pyramid clock; the rest of the tube
  at the hinted timing), fitted to zero flow through the link in front of
  the seal (`link_flow`, seal node held) over a window from the read
  closure.
- Window to 20335: link flow −146.7 dB; front clock 0.034 phones rms off
  the truth over the window. But a read just past the window jumps at the
  window's edge (the front beyond it is unconstrained).
- With the front's clock carried past the window at its last rate: the
  reads climbed 20431 → 21199 (220.8 ms, where the true front starts to
  move), then overshot to 28557 and 32007 (the moving front's leak
  outweighs the release's jump), and the windows that included open link
  could not be fitted (−9.1, −4.8 dB).
- Release as the latest window end at which zero flow can still be fitted
  (bisection; sealed fits reach −143 to −147 dB, open about −31; the test,
  within 100× of the established floor, is insensitive between 10× and
  1e4×), 20 Gauss-Newton steps a window:
  - sealed fits to 21303 (−146.8 dB);
  - degrading from 21319 (−108.3) to 21327 (−96.4);
  - release placed at 21319 (222.1 ms), true 24105 (251.1 ms).
- Reading: outside-in recovers the front's posture and its stillness
  exactly through the hold. It does not follow the front's motion once it
  begins (about 220 ms). The fit stalls there (a single clock can represent
  the truth, so it is the fit, not the representation), and the test
  takes "cannot fit" for "released". OPEN: fit the moving front (more
  steps, or the front's clock from the transition's own shape).

## Step 32: outside-in, the moving front. MEASURED.

- The same breakdown at about step 21320 (222 ms) with:
  - a direct, unclipped clock correction (DIRECT=1, no dead zone);
  - the seal node following the front's clock, capped at the seal (its
    area grows within the sealed range near the release).
- Decisive test: the true clock put into the front (direct correction)
  fits zero link flow to −148.9 / −149.4 / −149.4 / −149.9 / −150.6 dB for
  windows ending 20331 / 21303 / 21327 / 22339 / 24000. The front model is
  exact; the fit was failing.
- Direct-correction steps from each coefficient's own reading (the rate's
  steps were far too small): the window fit then holds to 21412–21420
  (−125.7 / −112.3 dB), the release is placed at 21420 (223.1 ms; true
  24105), and windows past it land at −73 to −102 dB from a warm start.
- Reading: once the front moves (from about 220 ms), zero link flow is a
  non-convex landscape in the front's clock. The warm-started fit cannot
  get from a held front to the moving one, though the true motion fits to
  −150 dB. Outside-in pins the front's posture and its hold exactly; its
  motion before a velar release is the open part. OPEN: a start for the
  moving stretch from the transition's own course (the hinted shape moved
  as a whole) rather than the held front carried on.

## Step 33: the front's motion onset is readable. MEASURED.

- ONSETSTART=1: after the hold is fitted, a read through the held front
  (carried at its last rate) gives the front's motion onset: 21199
  (220.82 ms; the true front starts at 220 ms). From there the moving
  stretch starts from the hinted transition's shape (its hinted duration).
- The bisection still breaks at about 21319 (−110.8 dB at 21319, −98.9 at
  21327): release placed at 21319, true 24105.
- Cause (in the set-up): the pyramid's coefficients are global in time,
  so fitting the earlier, held windows also rewrites the clock after the
  onset and erases the shape placed there. The front's motion onset reads
  right; the moving stretch needs a clock basis local in time (fitting the
  hold must not touch the motion), or its own coefficients fitted alone.

## Step 34: the moving front on a local clock. MEASURED.

- After the onset read (21199), the hold is fitted right up to it
  (−149.3 dB) and frozen, carried on at its last rate. The motion after the
  onset gets its own pyramid (381 coefficients, local in time), starting
  from the hinted transition's shape.
  - (Frozen straight after the first hold window instead, the hint's
    motion from 212.7 ms stood in the hold and every window sat at
    −24 dB.)
- Bisection: windows fitted to 21999 (229.2 ms; −109 to −134 dB), then
  failing (−73 to −87 dB); release placed at 21999 against the true 24105.
  The window to 23994 reached −83.5 dB.
- Continuation (the window grown a tick at a time, each fit warm-started
  from the last) failed at its first tick: the window to 21295 fitted only
  to −77.3 dB.
- Reading: the moving front's fit depends on how it begins (bisection's
  far-ahead first fits did better than tick-by-tick continuation). The
  front's posture, its hold and its motion onset are recovered exactly by
  outside-in; its motion over the last 22 ms before a velar release is not
  yet. OPEN.

## Step 35: the moving front's fit, bugs and the peel's reach. MEASURED.

- Check: the true motion put into the local pyramid fits zero link flow to
  −149.4 / −149.8 / −150.6 dB for windows ending 21295 / 22159 / 24000,
  with the frozen hold exact at the onset (1.0000 against 1.0000). The
  representation and the hold are right.
- In-process at the local start (no carry-on): −86.2 dB over the window to
  21295; nudging one local coefficient by 1e-3 moves it to −50.1. The
  coefficients do act.
- BUG, fixed: the carry-on point was the window's end rounded down to a
  tick, so the last partial tick inside the window, the one a continuation
  step adds, was overwritten. Now rounded up.
- Smooth-function difference steps (ε^(1/3)) for the local coefficients,
  since the link flow is a pure peel with no events: no change.
- Continuation still fails at the first tick (−71.7 dB). HYPOTHESIS: the
  front peel at a step reads the sound up to about 19 engine steps later,
  through the geometry at those later times, so the carried-on clock past
  the window still reaches the window's last steps. The window's fit and
  its carry-on are coupled through the peel's reach.
- Status: outside-in recovers a velar stop's front posture, its hold
  (−149 dB) and its motion onset (220.82 ms against the true 220) from the
  sound. The front's motion over the last 22 ms before the release is not
  yet recovered (best: bisection to 229.2 ms, true release 251.1). Paused
  here: the remaining issues are in the method's bookkeeping, not its
  physics.

## Step 36: the moving front, paused. MEASURED.

- Comparing the window only up to one front's length before its end (the
  front peel's reach, 19 steps): the first local window still fitted only
  to −77.9 dB. It also changed the hold. Its last steps were no longer
  scored, so the frozen hold came out at 0.9950 against the true 1.0000 at
  the onset, and that error, before the onset, dominated: a 1e-3 nudge to
  the local clock left −77.9 dB unchanged.
- The local pool's workers returned no columns (all zero, so the fit
  exited before its first step); computing the local columns in-process
  did not help, because the dominant error was then in the frozen hold.
- Status (paused at diminishing returns; the physics is verified, since
  the true motion fits to −149 dB):
  - outside-in recovers a velar stop's front posture, hold and motion onset
    from the sound;
  - fitting the front's motion over the last 22 ms before the release
    trades one piece of bookkeeping against another: the hold's freeze, the
    carry-on, the peel's reach, the local basis.
  - A clean formulation would fit the hold and the motion together on one
    clock local in time (a spline or B-spline basis instead of the global
    pyramid), with the peel's reach handled by fitting the window's clock
    one reach beyond its compared stretch.

## Step 37: a velar release recovered outside-in. MEASURED.

- `analysis/outside_in.py` (copied from scratch): the cavity in front of
  Echo's velar seal gets its own clock, the hinted clock plus a cubic
  B-spline correction, one knot a tick (local in time). The seal node
  follows it, capped at the seal.
  - Each window from the read closure is fitted to zero flow through the
    sealed link (`link_flow`, damped Gauss-Newton, ε^(1/3) central steps:
    the peel has no events).
  - Free coefficients: the last 10 ticks before the window's end plus one
    front's length past it (the front peel reads ahead; a local basis
    leaves the earlier fits alone).
  - The window grows a tick at a time, warm-started. The release is the
    first end where the fit leaves the sealed floor (within 100×;
    insensitive between 10× and 1e4×).
- Check: the true clock through the basis fits to −150.6 dB (window to
  24000).
- Echo, from the hinted timing:
  - hold −148.6 to −148.9 dB, front clock 1.0000 against the truth's;
  - slight slip as the hint's own early move is cancelled (0.9978 at
    20923, 0.9926 at 21307);
  - then the moving front tracked exactly: 1.0080 / 1.0319 / 1.0708 /
    1.1235 / 1.1883 / 1.2631 / 1.3455 / 1.3670 against the truth's
    identical values, link flow −114 to −125 dB;
  - at 24187 the fit broke to −32.6 dB.
  - Release found at step 24091 (true 24105): 14 steps, 0.15 ms, within
    the window's one-tick resolution.
- Reading: outside-in on a clock local in time recovers a velar stop's
  front posture, hold, motion onset, the motion itself and the release
  from the sound alone. The local basis was what step 36 lacked: on the
  global pyramid, fitting the hold rewrote the motion.
- NEXT: knot Echo's release at the recovered step in the full time-field
  inverse (step 17 had stalled 2 ms early without it).

## Step 38: Echo with the recovered release knotted. MEASURED.

- `fields_word.py` takes recovered events (EVENT_STEPS="phone:closure:release",
  −1 for none) and knots only those. FIXED: `align_events` and
  `event_steps` searched a release-only stop's sealed span at half the
  release step (written for word-initial stops, sealed from 0); now just
  inside the closure.
- Echo, release knotted at outside-in's 24091 (knot at tick 248.946; the
  true crossing is at 249.08):
  - start −12.9 dB, release −0.29 ms;
  - after 27 steps −28.1 dB (knot-free, step 17: −24.8), release −0.37 ms;
  - but the mouth clock drifted 0.028 → 0.051 phones rms and the folds'
    clock to 0.080.
- Reading: with the velar release readable, Echo is limited by the voicing
  onset after the voiceless /k/ into /ow/ (steps 25–30): the mouth clock
  bends to compensate the folds' phase-sensitive restart. Velar releases:
  solved by outside-in. Voiceless-to-voiced onsets: the remaining blocker.

## Step 39: both of Echo's events from outside-in; Echo's timing converges. MEASURED.

- Diagnostic: Echo with the folds' clock held at the truth and only the
  release knotted: −25.6 dB, mouth clock worse (0.028 → 0.037), the
  unknotted velar closure drifting +1.9 ms. Echo's other weak point was
  its closure.
- `outside_in.py`: after the release, the window's start is moved earlier
  a tick at a time while zero flow still fits (backwards, with the
  coefficients near the start free). Then each end is refined to the step
  by bisection.
  - Echo: closure 12246 (true 12244, 2 steps), release 24101 (true 24105,
    4 steps). Tick-level first: 12283 and 24091.
- Echo's timing inverse with both knotted (EVENT_STEPS="1:12246:24101";
  knots at ticks 124.465 and 249.024):
  - −13.2 → −45.8 dB in 28 steps (knot-free −24.8; release alone −28.1);
  - events within 0.05 ms; mouth clock 0.025 phones rms.
  - Mouth clock off by stretch: ≤ 0.0045 through /eh/; 0.156 at the start
    of the closure, 120–160 ms (unheard); 0.0006–0.0136 in and around the
    release; 0.0000 through /ow/.
  - Round trip 28.0 dB. By stretch: 56–73 dB through /eh/; 49.7–85.5
    inside the closure; 32.0 at the release; 22.6–24.7 through /ow/, where
    the timing is exact.
- Reading: with a velar stop's closure and release read outside-in, the
  inverse recovers the timing wherever it is heard. The voiced /ow/ after
  the voiceless /k/ round-trips at only ~24 dB with exact timing,
  confirming that the voiced waveform round trip measures the folds'
  phase, not the motion (step 22).

## Step 40: Oscar's velar stop, outside-in. MEASURED.

- `outside_in.py`: the true seal and its check window are now found from
  the read closure (they had been set from Echo's positions).
- Oscar (AA S~0.08! K! ER; JETBOUND, SEALBLEND): closure read 24793,
  against the true span (24792, 38904): exact to one step through the
  hint.
- But zero link flow with the true clock reaches only −64.6 dB (Echo
  −150). Inside the closure the error decays from 2.0e-3 just after the
  closure to 1.5e-5 by the release (time constant ~25 ms): a transient
  carried in from before the closure.
- `link_flow` now gives the front peel the word's jet: the narrowest node
  and area from the whole tube (mapped into the front's numbering, a
  narrowest node behind the front clamped to its first link) and the noise
  hashed by global link index. Result −60.8 dB, then −64.7: no change.
- Against the forward replica's true link flow (the replica matches the
  engine to 2.3e-13), the front peel is wrong before the closure: about
  1e-2 through /aa/ (25–100 ms, the jet's gate zero) and 2.7 against a
  true 1.5 during /s/ (150–225 ms). Echo's front peel had been checked to
  1e-12, but Echo has no fricative.
- OPEN, a replica-exactness bug: the front peel with a fricative in front
  of the seal does not reproduce the engine. Suspects: the jet's state in
  the side peel when its gate is zero, the side peel's narrowest table
  from the front-only geometry before the override, or the peel through
  the /s/ narrowing. Until fixed, outside-in on Oscar's /k/ cannot reach
  the sealed floor.
- Isolation (front peel against the replica's true link flow, first
  135 ms):
  - with and without the jet, identical through /aa/ (8.8e-3 to 1.5e-2 of
    the scale; the jet's gate is zero there). That is the reading's held
    seal-node area changing the link's admittance before the closure, as in
    Echo, whose closure still reached −150 dB;
  - from the /s/ onset the jet path is worse: 1.62e-1 against 1.54e-2 at
    100 ms, 2.64e-1 against 9.55e-4 at 125 ms.
  - So the whole-tube jet mapping in the front peel (narrowest clamped,
    global noise rows, the side peel's jet state) is inexact, and its error
    carries into Oscar's closure. OPEN: reproduce the engine's jet on a
    front-only peel exactly (the replica's jet source needs a narrowest node
    behind the front, not clamped, and its filter states per global link).

## Step 41: Oscar's front peel exact; closure exact; release blocked. MEASURED.

- The side peel's jet in global numbering (`kernels_jet.jet_source`'s
  `offset`, `SidePeel.link_offset`; the narrowest node and noise from the
  whole tube, no clamping) and with the engine's closing fade
  (`SidePeel.jet_closing` = JET_CLOSING when the word has a seal). The
  unheld front peel then matches the replica's true link flow to 1e-13
  throughout Oscar, /s/ and closure included.
- The remaining gap was outside-in's own command: the seal node capped
  over the whole word changed the link before the closure and, under the
  /s/ jet, became the narrowest node there. Capped only from the window's
  start: the true clock through the basis fits zero flow to −278.6 dB
  (before: −64.6).
- First window from the read closure: −34.8 dB (the /s/-to-/k/ transition).
  Anchored in the hold's middle with the past at the hint: −56.4 dB (the
  front peel's state carries the past). The first fit freeing the clock
  over 30 ticks before the window (HISTORY): −92.4 dB.
- Oscar: closure 24792 (true 24792, exact). Release 35649 against the true
  38904 (34 ms early). The front clock sat at 1.932 through /k/'s hold (true
  2.000; −86 dB, a sound-alike front posture, 7% toward /s/'s). When the
  true front starts toward /er/ (about 360 ms), the fitted clock runs
  backward, 1.93 → 0.97, until the fit breaks (−35.5 dB).
- OPEN: the front's clock must be monotone (progress never reverses, as
  the main time field's clipped rate guarantees) and anchored at its
  target through the hold; the direct B-spline correction is neither.
- REJECTED: the front clock's running maximum (monotone): its kinks broke
  the fit (−53 dB at 300 ms, clock 1.82–1.87), release at 28993.
- Anchored hold (coefficients wholly inside the hinted hold held at zero,
  the front on its target there): the fit holds to 28825 (−81 dB) and
  breaks at 28921 (301 ms), where the window enters the anchor. The clock
  just before is 1.94 against the truth's 2.0 (true arrival 290 ms). The
  slip toward the /s/ shape happens in the transition (260–300 ms), where a
  front 6% short of its target changes the link flow too little to be
  told apart: the posture non-uniqueness of step 16, inside the front
  cavity.
- Status, Oscar's /k/: closure exact (24792/24793); release not recovered.
  The front's transition from the fricative into the stop is sound-alike
  near its target. OPEN.
- Regression, Echo with the current outside_in.py (the past freed in the
  first fit, the cap from the window's start, the anchored hold): closure
  12245 (true 12244), release 24096 (true 24105). First window −152.0 dB;
  the true clock through the basis −288.1 dB.
- HOLDSTART=1 (the first window in the anchored hold, the past free):
  - 30 ticks free: −54.9 dB; the growth fails at once.
  - The past free back past the read closure: −95.2 dB; then the growth
    fails both ways early: backward at 311 ms, still inside the hold
    (29893, against the true closure 24792), and forward at 361 ms, where
    the true front starts to move (34711, against the true release 38904).
- Reading: the freed past is fitted to a sound-alike history that makes
  the first window consistent but is not the true one; moving the window
  into it exposes it. Oscar's release stays blocked by the front cavity's
  sound-alike posture as it leaves the fricative. The closure is exact from
  the read-closure start. Paused.
- From the read closure with a 40-tick free tail (RECENT=40): identical to
  10 ticks (release 28894; clock 1.9418 just before the anchor). Over the
  window the fit holds the 1.94 posture at −81 dB, where the true clock fits
  to −278 dB: a local minimum the closure-side start does not leave.
  Oscar's release paused.

## Step 42: Alfa, the mouth itself stuck at the onset. MEASURED.

- Alfa with the folds' clock held at the truth (FOLDS_TRUE=1): −24.5 dB,
  and no LM step kept. Not only the folds' timing.
- Mouth-only Jacobian at the start (true folds): cosines 0.54, −0.48,
  −0.39; 93% of the mismatch at 380–420 ms, the build-up before the true
  voicing restart (~410 ms). The folds' onset response depends on the tube's
  shape as voicing restarts, so the mouth's landscape is rough there too.
  Toward the truth the cost falls steadily (−24.5, −25.0, −26.2, −28.4,
  −29.1, −38.2, −243.4).
- QUANTUM=tick (a secant over one tick of clock, 5.6e-3 phones, to average
  the onset's ~5e-5-phone ripples): still no step kept from −19.5 dB.
- Reading: at every difference scale up to a tick the local gradient points
  the wrong way; the onset's roughness reaches at least that far, while the
  descending trend lives at the scale of the hint's own error (~0.1
  phones). Echo converged because its events were pinned, keeping the start
  near the truth; Alfa has no seal to pin. OPEN (the owner's call): a
  derivative scale tied to the start's uncertainty, or another anchor for
  seal-less words (the fricative's onset and offset).

## Step 43: Delta's /t/ by outside-in; the residual is a meter. MEASURED.

- Closure read exactly from the hinted start: 35913 (true 35912).
- Release: with the hint's own early move in the front's base clock the
  fit had to cancel it; the added ticks crept −348 → −108 dB and the
  release read 46659 (true 50047, 35 ms early). Fixed structurally: the
  base clock holds the stop's target from its hinted arrival on, anchors
  taken from the hint's own hold (not the held base, which froze every
  coefficient after the arrival). The fitted clock then tracked the truth's
  departure to the tick: 3.000386 at tick 491 (truth 3.000385).
- Release still read 47145 (30 ms early), because of the criterion, not
  the tracking. Probes of the zero-flow residual at the true clock (window
  to 47241):

  | clock | last tick |
  |---|---|
  | truth | −340 dB |
  | truth +1e-6 at one tick | −109 dB |
  | truth held at 3 from tick 492 (3.9e-4 short) | −113 dB |
  | truth held at 3 from tick 491 (1.5e-3 short) | −57 dB |

  The residual is a linear meter of front mismatch: a clock error d reads
  about 20 log10 d dB. The midpoint criterion (−111 dB, halfway from the
  −222 floor to the −0.3 open level) demanded a *moving* front fitted to
  3e-6 of a phone; the fit reached −115 to −135 dB. A release is a −0.3 dB
  event, 110 dB above anything the front can mimic within a tick. The
  detector confused numerical precision with physics.
- REJECTED today (knobs, removed): forward differences at the phone-index
  kink (same result), a floor at the stop's index, a smaller difference
  step (1e-8: same −135 dB). The fit's −135 dB floor against the truth's
  −340 is a solver matter, not physics, and is not worth another knob.
- Reading: what bounds a front's mimicry within a tick is how far a node
  can move in a tick — a per-node speed, the puppet's "friction". The
  detector, the sound-alike fronts (Oscar), the onset block (Alfa) and the
  null directions all want the same structure. Program turn below.

# Program turn (2026-10-24): the puppet as the inverse's family

The free time field (one knot a tick, hundreds of coefficients) is a fitting
machine: the sound barely hears most of its directions and the onset's
roughness sends it where it will. The owner's puppet plan — per-node
traverse speed, geometry that adapts and geometry that does not — is not a
later phase to run after telemetry; it is the prior that makes telemetry
well-posed. Phases:

- **A.** Close the tick era: knobs removed, step 43 logged (this entry).
- **B.** `analysis/puppet.py`: the inverse over event times (two a
  transition) with every column moving between targets by a half-cosine at
  its own lag and stretch (held at 0 and 1 here). Secant by decades from
  the hint's error (10 ms) to one engine step. HYPOTHESIS: Alfa descends
  from −19.5 dB; Echo, Golf, Papa do not regress. (Note: the pyramid truth
  moves every column, the folds included, on one clock — `command_of` — so
  the truth is in this family exactly.)
- **C.** Observability as telemetry: the Jacobian's singular spectrum over
  the few event parameters says what the sound determines and to what
  precision; the rest is filled by the puppet and flagged. Retires hints and
  anchors as priors.
- **D.** Friction across words: per-column lag and stretch fitted jointly
  over the word set, held-out words predicted from text and timing. No
  per-word tuning; cross-validated. This is "represent voices without
  exhaustive fine-tuning".
- **E.** Engine refinement from unphysical recovered motions only.

## Step 44: Phase B, the puppet family on the six words. MEASURED.

`analysis/puppet.py`, no knots, no events read, no outside-in: the hinted
event times (±7.3 ms alternating) as the start, every column on the shared
clock (lags 0), secant rungs 10 ms → 1 ms → 0.1 ms → one step, 8 LM
iterations a rung (ITERATIONS, the existing default). JETBOUND, SEALBLEND.

| word | before (best, how) | puppet | events off the truth |
|---|---|---|---|
| Alfa | −19.5 dB, no step kept (step 42) | −59.2 dB | 0.01–0.07 ms (1–7 steps) |
| Papa | −43 dB, true postures + knots | −90.8 dB | ≤ 4 µs |
| Echo | −45.8 dB, both events knotted from outside-in | −84.7 dB | ≤ 13 µs |
| Golf | −38.6 dB | −90.3 dB | ≤ 3 µs |
| Oscar | release blocked (sound-alike front, step 41) | −84.7 dB | ≤ 3 µs |
| Delta | closure exact, release detector failed (step 43) | −52.0 dB | 6 of 8 ≤ 8 µs; /l/ leaves +1.56 ms, /t/ arrives −7.67 ms |

- Every run was still descending or at its floor within about a minute of
  Mini time; the free time field took hours and stalled.
- HYPOTHESIS of the program turn: CONFIRMED on all six. The event family
  averages the onset's roughness over a move; the sound hears each event's
  midpoint best (the largest singular direction of every word is a same-sign
  pair, shifting a whole move) and each move's duration least.
- Delta's /l/→/t/ pair moved along the one direction that keeps the seal
  crossing where it is: the closure's completion past contact is unheard
  (the observability map's closures-unheard, now as a parameter direction).
  This is Phase D's case: with the /t/ node's traverse known from other
  words, the arrival follows from the leaving. Until then its precision is
  what the spectrum says.
- Alfa's floor (−59.2 dB, truth −238.9 in the family) and Delta's (−52.0,
  truth −247.6) are the two with voicing onsets inside the error. OPEN:
  whether more iterations lower Alfa (run queued) or the onset's 3 µs
  non-monotone response (physics) sets the floor. Not to be fixed by a
  knob; to be reported as precision.
- The precision line (spectrum): the move of each event, in µs, that
  raises the residual's sum of squares by its own remaining share
  (s sqrt(diag (JᵀJ)⁻¹), s² = |r|²/(N − p)). Oscar: 0.00–0.03 µs at
  −84.7 dB. A statement of what the sound fixes, not of the truth.

## Step 45: the cost we depend on. MEASURED (Mini, Papa, 51,840 steps).

- One equation-error residual: 0.29–0.50 s. Three full engine renders
  (`engine.render`, 67 ms each: FieldFlow for the hidden stretches and the
  forced folds) 0.20 s; the peel (`observer.implied` → `_implied_run`)
  0.19 s; `moving_geometry` 0.02 s; the puppet's command 2 ms.
- A Jacobian column is two residuals. Puppet words: 6–8 columns; the voice
  problem (`voice.py`, five words, ~100 shared traverse unknowns): ~1,000
  residuals a Jacobian.
- Owner direction (2026-10-24): a worker is to refine this algorithm and do
  the formal mathematics toward a better representation (the waveguide's
  adjoint for all sensitivities in one backward pass; checkpointed renders
  from the last unchanged tick; the chain to event times, lags and extra
  move times, which is analytic through the half-cosine). Its record:
  `NOTES_adjoint.md`, module `analysis/adjoint.py`; engine changes behind
  switches with byte-identical defaults.

## Step 46: the floors are local minima; precision is local. MEASURED.

- Alfa with 24 LM iterations a rung (ITERATIONS=24): −59.2 dB, the same
  events as with 8 (F leaves +0.023 ms, AH arrives −0.070 ms). Not the
  iteration cap: a local minimum. Between rungs the minimum hopped (at the
  0.1 ms rung F leaves −0.071, AH arrives +0.009): the landscape is rough
  at the tens-of-µs scale around the voicing onset — the folds' forced
  response near threshold (step 42's 3 µs non-monotonicity) sets the onset
  phase, and the onset phase sets the residual (a 30 µs fold-phase shift
  alone costs 14.5 dB in the voiced round trip).
- The spectrum's precision at Alfa's point: F− 0.07 µs, AH+ 0.10 µs — while
  both are 23 and 70 µs off the truth. At Delta's: T+ 3.07 µs while 7.67 ms
  off. The precision is the *local curvature*: a statement about the
  minimum found, not about the truth. Labeled so from here on. A rough
  landscape makes the local statement sharp and wrong at once.
- Delta's precision otherwise: D− 2.2, EH+ 2.2, EH− 7.4, L+ 6.8, L− 0.6,
  T− 0.5, AH+ 0.9 µs.
- OPEN, being measured: the residual along AH+ (Alfa) and T+ (Delta) through
  the found point and the truth, ±200 µs at 2–10 µs spacing (`landscape.py`
  in /tmp/esf): how wide is the truth's basin, how many minima.
- HYPOTHESIS (structural, for sealing transitions): parametrize each
  sealing move by its seal crossing and its duration instead of (leaving,
  arrival). The crossing is what the sound hears sharply; the 10 ms rung's
  linear model slid Delta 15 ms along the crossing-preserving pair and the
  finer rungs could not bring it back.

## Step 47: the landscape along one event. MEASURED (`landscape.py`).

Alfa along AH's arrival, through the truth (dB of the equation error):

| offset | −200 µs | −40 | −20 | −10 | −4 | −2 | 0 | +2 | +4 | +10 | +20 | +40 | +200 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| level | −37.3 | −51.7 | −46.2 | −50.3 | −57.6 | −63.4 | −238.9 | −63.1 | −57.1 | −49.3 | −43.5 | −36.2 | −39.8 |

- A needle: the basin is ~4 µs wide at −60 dB; beyond ±10 µs the surface is
  rough with false minima (−51.7 at −40 µs; through the found point −62.6
  at +50 µs, which is still 20 µs short of the truth and not on the way).
  LM's local model from 70 µs away sees only the roughness. Not a solver
  fault to fix with a knob: the measured shape of the problem.
- Shape: the level is linear in log|offset| (−63 → −57 → −50 over 2, 4,
  10 µs): the residual grows linearly with the offset, as a waveform
  shifted in time does — the voiced flow's phase (the fold phase that a
  30 µs shift alone put at 14.5 dB in the round trip).
- Delta along T's arrival: flat at the found point (−52.0 ± 0.2 over
  ±200 µs: the floor masks it), sharp at the truth (−101 dB at ±2 µs,
  −62 at ±200). The duration is heard; the found point is a wrong
  approach shape (/l/ leaves 1.56 ms late, the move 9.2 ms short) with the
  crossing kept.
- HYPOTHESIS (structural): the voiced stretch after an onset is, to first
  order, our folds' flow shifted in time against the implied flow. The
  shift is a smooth, linear meter of the onset's phase and can be measured
  directly (cross-correlation to sub-sample precision), where the equation
  error is a needle. Newton on the measured shift, with a coarse secant for
  d(shift)/d(event), should walk to the needle. Test: Alfa's found point.

## Step 48: the onset-shift meter REJECTED; the landscape has three scales. MEASURED.

- The shift of our driven flow against the implied flow over 5–60 ms
  after AH's arrival (least-squares d − i = τ di/dt), through the truth:
  +2 µs → −0.09, +10 → −0.86, +30 → −4.67, +70 → −2.58, +150 → +0.02,
  −70 → −0.33 µs. Not linear, not monotone, a few percent of the offset.
  Newton on it from the found point drifted −70 → −80 µs. REJECTED: the
  folds are forced by the implied returning wave, so the driven flow stays
  phase-locked to the implied flow; the residual is a difference in forced
  response (the onset transient), not a free-running phase.
- The landscape's three scales (step 47's scan with the rungs' record):
  a smooth bowl at ≥100 µs (the 0.1 ms rung brought AH's arrival to +9 µs,
  −41.8 dB); a rough layer at 10–50 µs with false minima deeper than the
  bowl (−50 to −63 dB against −49 at +10 µs: the one-step rung moved the
  event from +9 to −70 µs while lowering the residual to −59.2 dB); the
  needle within ±2 µs (−240 dB).
- Reading: below about −50 dB the equation error at a voicing onset is
  fitting the onset transient's roughness, which is the engine's forced
  folds' response near threshold (step 42: non-monotone at ~3 µs, physics
  of the model). On real speech no needle exists; the bowl is the
  information: an onset to ~10 µs. Telemetry should report the rung at
  which an event stopped moving coherently; a finer rung that moves it by
  more than the coarser rung's precision is fitting roughness.
- CANDIDATE for Phase E (engine refinement from strange behaviour): if the
  roughness comes from a hard contact or threshold in the fold model, a
  compliant contact (tissue) would make the onset's response smooth — the
  same C0-to-C1 move as the seal clip (step 43). Being checked in
  `src/voice_body.cpp`; the roughness grain measured at 0.5 µs.

## Step 49: the grain of the landscape. MEASURED (`grain.py`, Alfa along AH's arrival, 0.5 µs).

| offset (µs) | 0.5 | 1 | 2 | 4 | 7 | 12 |
|---|---|---|---|---|---|---|
| level (dB), either side | −75.3 | −69.2 | −63.3 | −57.4 | −52.7 | −48.9 |

- Within ±12 µs the surface is smooth, symmetric and exactly linear in the
  offset (6.0 dB a doubling): the meter of step 43, now at a voicing
  onset. No roughness there at all. Step 48's "roughness" is the waviness
  at 20–100 µs (−43.5 at +20, −36.2 at +40, −43 at +90 µs; −51.7 at −40),
  and the false minima of −59/−63 dB were found with *other* events 23–71 µs
  off: the coupling, not the single coordinate.
- So the picture is: a linear region of about ±12 µs around the truth in
  each event; a wavy band at 20–100 µs with a ~50 µs grain (source OPEN:
  not a glottal cycle, 8 ms; not one step, 10.4 µs); the smooth bowl
  beyond. The 0.1 ms rung lands inside the band (AH +9, F −71 µs); the
  one-step rung then resolves the band and can settle in a false minimum
  (F +23, AH −70 µs, −59.2 dB).
- The fold contact is already compliant (`voice_body.cpp:215`,
  contact_softness_2025, the Geometric glottis 2025 mathematics): step 48's
  hard-threshold CANDIDATE is weakened. The waviness's source remains OPEN.
- Recommendation (for the owner; reversal cost nil): report onset-coupled
  events at the band's scale, ±70 µs, as the telemetry's honest precision;
  seal crossings at the few-µs the runs give. 70 µs is 7% of a control
  tick and far inside anything the puppet plan (per-node traverse at the
  ms scale) needs. Chasing the needle below the band on a known-answer word
  would be fitting the model's own onset transient.

## Step 50: the worker's refinement of the cost (2026-10-24). MEASURED; record in `NOTES_adjoint.md`.

- Step 45's brief was wrong on one point: the engine already has snapshot
  and resume (`engine.Run`, `tract_run_copy`, `Settings.divergence`,
  `engine.Incremental`), as do the peel (`observer.implied_resume`) and the
  geometry; `equation_error` used none of them. No C++ was needed.
- `analysis/adjoint.py`: `ResidualCache` resumes each evaluation from the
  last at the first changed command tick. Exact mode bit-identical to
  `PuppetProblem.residual` on 51/51 Papa and 51/51 Alfa points and every
  Jacobian column tried; CPU 1.6× less per residual, 1.55–1.75× per
  Jacobian. A canonical form (the folds free inside every earlier hidden
  window) runs one engine pass a residual, differs by ≤ 2e-12 relative,
  2.5–3×. Not yet adopted by `puppet.py`/`voice.py`: to be verified on the
  other four words first (below).
- The derivative's honesty (its §7): pointwise forward and central
  differences agree only for secants ≤ 0.3–1 µs; between 0.3 and 3 µs every
  event column meets a whole-step jump (the ear's closures, the hidden
  windows' edges — HYPOTHESIS on the source), after which difference norms
  fall as d^−1/2. Consecutive puppet rungs measure directions 30–60° apart.
  This is the mechanism behind step 49's wavy band: the residual is
  piecewise smooth with step-timed jumps, and the rungs are secants across
  them. An exact adjoint would differentiate the rough object; REJECTED as
  the puppet's Jacobian at the present rungs; derived and costed for Phase
  C telemetry or for an engine whose whole-step events are timed within the
  step (the owner's call).
- OBSERVED by the worker: `fields.read_postures` reads every side column
  through the last column's clock, so lags of fold columns 46–47 move
  nothing (their Jacobian columns are exactly zero; `voice.py` names them
  unheard).
- Verified on the other four (`adjoint.py check`, Mini): Delta 67 points,
  Echo 35, Golf 51, Oscar 51 — 0 differ in exact mode; CPU plain → exact
  49.7 → 31.2 s, 7.7 → 3.9, 25.3 → 11.8, 30.2 → 14.7 (1.6–2.2× less). The
  canonical form differed by 2.6e-12 on Delta and not at all on the three
  others, and gained over exact only on Delta (17.3 s). DECIDED: the exact
  resumable residual is `puppet.py`'s and `voice.py`'s default
  (`PuppetProblem.cache`, opt-out `plain`); the canonical form stays a
  CANDIDATE. Reversal: one word on the command line.

## Step 51: Phase D's first attempt — the truth was not in the voice's family. MEASURED, fixed.

- `voice.py` over five words (Alfa, Papa, Echo, Delta, Golf), 47 shared
  node lags: "the truth in the family −23.9 dB". Cause: the words' node
  counts differ (Papa and Golf 46, Echo and Delta 44 — the count follows
  the first vowel's length), so a lag named by its node from the glottis is
  a different place in each word. Fixed: articulators named by distance
  from the lips (−231.9 dB; with extra move times −236.0).
- Everything the run then did was fitting an inconsistent model: lags of
  −57/+21 ms where the truth is 0, events worse than the hint (Delta
  29 ms). Discarded, kept here as the negative result.
- Single-word friction (Alfa, 47 free lags + 6 events): the lagged truth in
  the family at −224.9 dB, but the fit from the hint reached only −31.1 dB
  with events 1–5 ms off and lags ±20 ms wrong; the lags' local precisions
  0.5–86 ms. One word cannot determine 46 node lags: the lags absorb what
  the events should explain. Structure, not a knob: the shared unit is the
  articulator group (lips 2 nodes, blade 6, front 8, back 10, pharynx, velum,
  folds — `voice.GROUPS`, by distance from the lips), six unknowns for the
  five words. Caveat: the synthetic pattern (`friction_truth`) was drawn on
  the same anatomy, so a grouped recovery is a known-answer check of the
  machinery, not evidence about real articulators.
- Delta with (crossing, duration) for its three sealing transitions:
  −47.2 dB; the crossings D|EH, L|T, T|AH to 0.04, 0.23, 0.07 µs (local);
  the /l/→/t/ duration still wrong (+11.7 ms, after the 10 ms rung moved it
  +21 ms), at a different minimum from step 44's. The approach from a
  narrow lateral to closure is barely heard; its duration is multimodal and
  shallow. The traverse from other words is its proper source (Phase D).
- OBSERVED (worker): fold columns' lags move nothing (`read_postures` reads
  every side column through the last column's clock); the folds' lag is
  therefore carried by that last column only.

## Step 52: Phase D on grouped articulators — identifiability. MEASURED.

Five words, six shared lags (truth: blade −8, lips +6, folds +4 ms, the
rest 0), from the hint:

| run | level | lags found (ms) | events |
|---|---|---|---|
| joint, lags | −23.8 dB, stalled at the 10 ms rung | pharynx −1.9, back +0.8, front −3.3, blade −2.3, lips +4.5, folds +3.8 | 1.4–22 ms off |
| joint, lags + extras | −31.9 dB | lags 7–12 ms off; extras to 11 ms off | 7–11 ms off |
| hold-out Delta, true traverse | −53.2 dB | — | 6.9 ms (the /l/→/t/ duration) |
| hold-out Delta, learned traverse | −22.9 dB | from the joint run | 15.9 ms |
| hold-out Delta, no traverse | −21.9 dB | zero | 14.3 ms |

- The cross-validation machinery works and shows the stakes: 31 dB between
  the true traverse and none on a held-out word. The learned traverse was
  wrong because the joint solve stalled at the first rung.
- Why (OBSERVED in `sequence_recover.build_word`): a stop's or
  narrowing's tube is the vowel before's with one node changed, so in a
  stop transition one node and the folds move; a group's lag and the
  transition's event times are then collinear within the word. Only a
  transition where several groups move (a vowel change: P→AH, K→OW, T→AH,
  F→AH) separates them, about one per word. The joint landscape has long
  curved valleys mixing lags and events; LM from ±7.3 ms stalls in one.
- Structural answer (not a knob): lags are observable where articulators
  move together. A vowels-only word (IY, AA, UW: every node moves at every
  transition, each group with its own lag) should determine the group lags
  from one word; the stop words then contribute their events with the lags
  pinned by sharing. Running: `voice.py "IY,AA,UW" groups`, then six words.
  Real speech is rich in such transitions; the synthetic set was not.

## Step 53: group lags from a vowels-only word. MEASURED.

`voice.py "IY,AA,UW" groups` (every node moves at both transitions):
−36.2 dB; events 0.75 ms off; blade −7.69 ms (truth −8), lips −1.68 (+6),
pharynx +2.77 (0), front +3.51 (0), back +0.67 (0). The 0.1 ms rung ended
"no damping lowered it"; the one-step rung moved nothing.

- The blade's lag is observable and found; the lips' is not found, traded
  against spurious pharynx and front lags and 0.75 ms of event error. A
  valley at −36 dB (the truth −248.7 in the family).
- Reading: in every joint run (steps 52, 53) LM stalls at −24 to −36 dB
  where event-only fits reach −50 to −90 dB. The lag/event pairs are
  near-collinear; resolving a long valley's floor needs a Jacobian whose
  columns are consistent to within the collinearity, and the secant
  Jacobian is not (NOTES_adjoint.md §7: every column meets whole-step jumps
  at 0.3–3 µs; consecutive rungs measure directions 30–60° apart).
- Structural answer, not a knob: the problem is separable. Variable
  projection (Golub–Pereyra): for a candidate traverse, solve each word's
  events exactly (the event-only puppet fit, which converges to µs); the
  outer problem is then over the 5–6 shared lags alone, smooth at the ms
  scale because the collinear direction is solved away inside. Cost: an
  outer Jacobian column is one inner solve per word. Building:
  `voice.py ... varpro`.

## Step 54: the traverse from six words. MEASURED.

| run | level | blade (−8) | lips (+6) | folds (+4) | others (0) | events |
|---|---|---|---|---|---|---|
| vowel word alone, joint (step 53) | −36.2 | −7.7 | −1.7 | — | −3.3…+3.5 | 0.75 ms |
| vowel word alone, variable projection | −34.2 | +5.0 | −8.6 | — | −4.2…+3.9 | 1.6 ms |
| six words, joint (IY-AA-UW + the five) | −28.1 | −6.46 | +5.73 | +4.42 | ≤ 0.44 | 0.8–12 ms |
| hold-out Delta: five words' training | — | +2.4 | +5.05 | +4.74 | −4.7…+2.9 | — |

- Six words determine the group traverse to 1.5 ms from a zero start. One
  vowel word does not: its two solutions are a mirror pair (blade and lips
  exchanged in sign, −34 dB), and variable projection chose the mirror.
  Sharing across words breaks the pair. REJECTED for now: variable
  projection as the cure for the stall (it was not the cause).
- Gauge (OBSERVED): adding c to every lag and subtracting c from every event
  leaves the command unchanged. "Traverse off the truth" is read in the
  truth's gauge (three groups at 0); report lags relative to the pharynx
  from here on.
- Without Delta the blade's lag is not found (+2.4): the blade's evidence is
  Delta's T→AH, the one transition where the blade and the body move
  together. The hold-out therefore says what it should: a traverse is
  learned from the transitions that exercise it, and no others.
- The six-word fit's remaining error is on the events (Delta 9.5 ms, Golf
  12.3 ms) where the event-only fits reach µs: the joint LM stalled with
  the traverse nearly right. Next: alternate event-only and traverse-only
  solves from this point (`voice.py ... alternate`).

## Step 55: alternating solves from the six-word traverse — a fixed point. MEASURED.

Three rounds of (events with the traverse held; traverse with the events
held) from step 54's traverse: −35.5 dB at round 0 and unchanged after;
traverse 1.54 ms max off (relative to the pharynx: blade −6.22 (−8), lips
+5.98 (+6), folds +4.67 (+4), back +0.68, front −0.08); events: Papa 0.73,
Echo 0.72 ms, vowel word 1.6, Golf 1.5, Alfa 2.5, Delta 7.0 ms.

- A true local minimum of the joint problem: neither block moves. With the
  true traverse the same words fit to −50…−90 dB each (step 44), so the
  last 1.5 ms of blade lag costs 20–50 dB and is compensated by events —
  Delta's /l/→/t/ duration above all (7 ms), the direction the sound
  barely hears (steps 44, 51), in the word that carries the blade's
  evidence (step 54's hold-out).
- Reading: the traverse is learned from the transitions that exercise each
  articulator together with others; the synthetic set gives the blade one
  such transition (T→AH), and that word's own weak direction eats the
  error. Not a solver fault; a design fact about the data. Real speech
  offers many such transitions per articulator.
- Next (running): nine words — the six plus IY,T@alveolar!,AA (the blade
  releasing into a wide move), UW,D@alveolar,IY (the blade sealing and
  releasing across a vowel change), AA,P@lips!,IY (the lips across a vowel
  change) — joint, then alternating from its traverse. HYPOTHESIS: the
  blade's lag closes to < 0.5 ms when two more transitions exercise it.

## Step 56: nine words — the blade's lag is biased, not undersampled. MEASURED (Mini, after its reboot; POOL=6, niced).

Nine words (the six plus IY-T-AA, UW-D-IY, AA-P-IY), group lags relative to
the pharynx (truth in parentheses):

| run | level | back (0) | front (0) | blade (−8) | lips (+6) | folds (+4) |
|---|---|---|---|---|---|---|
| six words, alternating (step 55) | −35.5 | +0.68 | −0.08 | −6.22 | +5.98 | +4.67 |
| nine words, joint | −23.6 | +1.10 | −0.03 | −5.92 | +5.74 | +3.65 |
| nine words, alternating from zero | −31.8 | +1.09 | −0.03 | −5.93 | +5.75 | +3.65 |

- HYPOTHESIS of step 55 (more blade transitions close the blade's lag to
  < 0.5 ms): REJECTED. Two more blade transitions moved it −6.22 → −5.93.
  The bias is the same in every run: a structural cause, not too few data.
- The events carry the error: Delta 6.9 ms, Alfa 5.2, AA-P-IY 2.5 ms (the
  joint fit had UW-D-IY at 16.6 and Delta at 29 ms). These are the sealing
  moves' durations. In this family a stop's tube is the vowel before's with
  one node changed (`build_word`), so in a vowel→stop move only the seal
  node (and the folds) move; past the seal crossing the rest of the move is
  unheard, so the move's duration is heard only through the node's
  approach. A free duration per transition is a direction the sound barely
  fixes (steps 44, 51), and each such direction absorbs part of the
  shared lag.
- Reading (structural): durations should not be free per transition. In
  the owner's puppet plan the traverse speed belongs to the articulator
  ("friction": how fast each node moves). Then a move's duration follows
  from the articulator and the distance it travels, and each transition
  keeps one free time (its crossing or midpoint). That removes exactly the
  unobservable direction and gives Phase D what it was meant to learn.
  CANDIDATE, for the owner's decision: duration = the articulator group's
  traverse time (shared across words), times a distance law to be chosen
  by measurement, not fitted.

## Step 57: the traverse belongs to the articulators (owner, 2026-10-24: "yes, do the per-articulator traverse speed"). DECIDED; runs MEASURED below.

`puppet.TraverseProblem`: each transition keeps one free time, its
midpoint m_k; column c moves through every transition in its articulator
group's own move time D_g, starting at m_k + lag_g − D_g/2. A column's
clock is the sum of its half-cosine steps (`half_cosine_steps`): equal to
node_progress where moves do not overlap (8.9e-16), continuous and
monotone where they do (coarticulation; node_progress jumped there).

- Distance law: a constant move time per articulator, whatever the
  distance — the mass-spring relation of speech kinematics (peak velocity
  proportional to displacement). CANDIDATE law, to be tested on real
  recordings; not fitted.
- Known-answer truth in the same family (`retime_truth_traverse`): the
  plan's midpoints, `friction_truth`'s lags, `traverse_truth_durations`'
  move times (blade 50, folds 60, lips 70, front 80, back 100, velum 110,
  pharynx 120 ms). Truth in the family: Papa −250.8, Delta −248.5,
  IY-AA-UW −249.8 dB. The hint: midpoints ±7.3 ms alternating, every move
  time at the plan's mean move (65 ms), lags 0.
- Unknowns: nine words 23 midpoints + 6 lags + 6 move times = 35 (the
  event family had 46 events + 6 lags). The dropped direction is each
  transition's free duration, the one the sound did not fix (step 56).
- `voice.py ... traverse` (mode replaces the with_extras flag: lags,
  extras, traverse). Cache, alternate, varpro and hold-out work in every
  mode. `adjoint.ResidualCache.residual` now asks `problem.valid_vector`.
- Caveat (as before): the truth is drawn from the family it is solved in;
  this measures identifiability and the solver, not the law.

## Step 58: the traverse family — four groups found, a pharynx/back exchange. MEASURED.

- Delta alone, midpoints free, the traverse held at the truth: −109.0 dB,
  every midpoint exact to the µs from ±7.3 ms. In the event family the same
  word stalled at −52 dB with the /l/→/t/ duration 7.7 ms off (step 44): the
  free duration was the obstacle, as step 56 read.
- Nine words from the hint (23 midpoints + 12 shared), joint: −25.3 dB.
  Raw lags (truth's gauge): blade −8.24 (−8), lips +4.99 (+6), folds +6.01
  (+4), front −1.86 (0), back +0.86 (0), pharynx −7.78 (0). Move times:
  blade 50.1 (50), lips 72.4 (70), front 75.8 (80), folds 65.8 (60), back
  124.6 (100), pharynx 87.0 (120). Midpoints 0.25–1.5 ms off (Alfa 5.5).
  The alternating solve from there did not move (−23.8 dB).
- The blade, the step 56 problem, is now found (−8.24, 50.1 ms): the
  structural reading of step 56 is CONFIRMED.
- The pharynx and back look exchanged: back long and early, pharynx short
  and late, against the truth's pharynx 120 / back 100. The two groups move
  together in every vowel move and nowhere else; from one shared start (65
  ms each) the first steps decide which becomes the longer. A mirror pair
  again, like the vowel word's blade/lips (step 54).
- Pharynx-relative lags are the wrong report when the pharynx is the worst
  determined group; read raw lags (gauge fixed by the truth) or relative to
  the best-determined group.
- Test running: the alternating solve from the found traverse with the
  pharynx's and back's (lag, move time) exchanged. If its level falls well
  below −25 dB, the remedy is structural: after convergence, try each
  pairwise exchange of group traverses and keep the one the sound prefers
  (a discrete choice made by the residual, not a constant).

## Step 59: the cost of a voice solve. MEASURED (Mini, nine words, traverse family, from the hint).

| solver | wall | Jacobians | all-word residuals (trials) | level |
|---|---|---|---|---|
| step 58 (Jacobian every step, columns spread over 6 workers) | 715 s | 32 | 32+ | −25.3 dB |
| fast path (`voice.levenberg_broyden`, one task a word, 9 workers) | 355 s | 11 | 74 | −23.7 dB |

- Fast path: one pool for all rungs; each word's Jacobian rows in one task
  (its own midpoints and the shared columns its articulators use, every
  residual resumed in that worker's cache); trials evaluated a word a
  task; Broyden's rank-one update after an accepted step, a full Jacobian
  only when the gate rejects a step on an updated one (11 Jacobians for
  ~25 accepted steps); `report` reads word levels from the residual it has
  (it re-evaluated nine residuals a step); the canonical one-pass residual
  is the default (`puppet.attach_cache`; 0.21 s CPU against 0.32 exact and
  0.41 plain on Delta; ≤ 2.6e-12 relative from plain). 2.0x, not more.
- Where the 355 s go: about 3,600 word residuals (Jacobians ~270 each,
  82% of them the 12 shared columns × 9 words × 2) at ~0.9 s CPU each.
  The M4 has 4 performance and 6 efficiency cores; with 9 workers most
  residuals run on efficiency cores (0.21 s was measured on one, unloaded,
  likely a performance core).
- The same pharynx/back exchange (pharynx 86, back 114 ms; truth 120,
  100); blade 50.0 ms and −7.71 ms, lips 72.0 and +5.43, folds 63.0 and
  +5.68, front 75.5 and −1.89. The exchanged start (stopped at round 1,
  for the cost work) put the move times right (pharynx 123, back 94) but
  the back's lag at −7.1 ms: not yet resolved.
- Next levers, in order of measured headroom: the residual itself (the
  C++ render is 41% of a canonical residual at ~33 ns a node-step; the
  numba peel 29%, Python bookkeeping 10%, geometry 8%); forward
  differences at the coarse rungs (half a Jacobian); and amortization —
  a voice's traverse is learned once, a new utterance then solves only its
  midpoints (Delta: 4 midpoints, 39 s).

## Step 60: shapes from sound — the experiment built. MEASURED (smoke tests); runs pending.

`analysis/shapes.py` (`ShapesProblem`, a `voice.VoiceProblem` in traverse
mode): each vowel's posture is one set of the body's 9 parameters shared by
every word and transition it appears in; a stop's tube is its vowel's with
the seal node closed. A posture's difference step at a rung is the same
fraction of its field's spread over the vowel table as the rung's time
secant is of the plan's mean move (15% at 10 ms, 0.016% at one step).
`voice.VoiceProblem.columns_of` and `steps_for` are now methods a subclass
extends; `shared_extras` is bounded (it ran into appended unknowns).

- Truth in the family (IY-AA-UW and Papa): −250.1 dB. Starts: 20% of the
  way to the hints, tubes 0.17–0.21 cm² rms off; the hints, 0.83–1.05.
- Frozen lip fields cap what can be reached (everything else at the
  truth): protrusion and width at the start −38.7 dB (20%) / −22.7 (hints);
  width alone −36.5 / −21.9; protrusion alone −39.9 / −23.5; none frozen
  −250.1. Step 13's near-nulls are combinations, not fields: REJECTED,
  freezing a heard field. Nothing is frozen.
- Plan (runs wait for the speed worker, to keep its timings honest): nine
  words, 7 vowels × 9 = 63 posture unknowns, from 20%, 50% and 100% of the
  way to the hints; first with the traverse held at the truth, then free.
  Compare with steps 11 and 15 (one word: 20% partly converged, the hints
  failed).

## Step 61: the worker's residual speed-up, checked and adopted. MEASURED.

Record: `NOTES_residual_speed.md`. Native peel, geometry and hidden-window
kernels (`NATIVE_PEEL=1`), and the engine's `ENTERINGONLY` (`ENTERING_ONLY=1`):
1.6–2.4x less CPU a canonical residual over ten words (sum 1.229 → 0.657 s),
199 of 199 points bit-identical, default renders 22 of 22 byte-identical,
a nine-word voice Jacobian identical (warm wall 6.98 → 4.05 s). The engine
is now instruction-bound at ~16 ns a node-step over 77 nodes (the trachea
runs too); 5x would need the model changed, which bit-identity forbids.

- My check on the shapes residual (three words, 47 midpoint moves and a
  posture move, `/tmp/esf/spot2.py`): 0.138 → 0.068 s CPU a residual,
  bit-identical over 2,827,200 values. (A first check read no gain: the
  Mini's zsh does not word-split `env $sw`, so NATIVE_PEEL was "1
  ENTERING_ONLY=1" and off. Invalid, redone.)
- DECIDED: both switches on in every run from here (set at launch; the
  defaults stay off so the library's other callers are untouched).
- Cumulative today on the residual: plain 0.41 → exact 0.32 → canonical 0.21
  → native 0.09–0.10 s CPU (Delta); on a nine-word solve, 715 → 355 s before
  the native kernels.

## Step 62: why a run "took three hours", and the rework. MEASURED.

- The 100% shapes run sat idle 19.5 h. Cause (crash report
  python3.12-2026-10-25-013314.ips): a pool worker aborted inside `fork()`
  (`_malloc_fork_parent` → `_os_unfair_lock_unowned_abort`) while starting
  a `geometry_probe` process for a new posture — macOS's malloc fork
  handler in a process with other threads running. `multiprocessing.Pool`
  then waited forever for the lost task. Not the native kernels.
- Fixes: `posture_recover.ProbeServer` keeps one probe process a Python
  process, fed a posture a line (0.6 ms a posture against 4.4 ms a fork);
  the voice and shapes solvers use `concurrent.futures.ProcessPoolExecutor`,
  whose dead worker raises `BrokenProcessPool` instead of hanging.
- Precision: the probe printed areas at four decimals and arcs at three;
  full-precision tubes (`geometry_probe full`, the default output
  byte-identical) differ by up to 1.27e-3 cm², more than the finest
  posture steps move them. `shapes.py` now builds the truth and every
  tube at full precision.
- Full-precision shapes runs (traverse known, 9 workers): from 20% of the
  way to the hints 164 s, tubes 0.10–0.30 cm² off, −25.5 dB; from 50%
  181 s, 0.21–0.42, −15.6 dB; from the hints 176 s, 0.35–1.31, −16.7 dB.
  Every start moves far toward the truth; none converges. OPEN: stall or
  basin.
- Chunked Jacobian (`voice.jacobian_tasks`: each word's columns cut into
  chunks, about three tasks a worker, longest words first): the 20% run in
  135 s on 6 workers against 164 s on 9, results identical to the digit.
  A nine-word solve: 715 s (step 58) → 355 (step 59) → 164 → 135 s.
- Owner (2026-10-26): the Mini crashes when overloaded. A 10-worker pool
  drove its load to 26; stopped by PID. Cap: POOL ≤ 6 on the Mini, niced.
- The x86 compute pod (4 threads of a Threadripper 7960X, 8 GB): the
  worker's NEON code had no platform guard. `include/tract/lanes.hpp`
  is `<arm_neon.h>` on ARM (the Mini's libraries rebuilt byte-identical)
  and the same intrinsics on SSE2 elsewhere (FMA3 or std::fma for vfmaq);
  `CMakeLists.txt` takes -march=native on x86. The three vector files
  compile on the pod; the build and smoke test are running.

## Step 63: memory, the 8-vCPU pod, and everything free. MEASURED.

- A warm worker held 2.1 GB (nine words' caches, ~230 MB each: engine
  snapshots, peel records, moving geometry; measured on the pod with
  `resource`). Four workers and the main process overran the 4-vCPU pod's
  8 GB (a worker killed: `BrokenProcessPool`, the robust pool working as
  meant); six on the Mini's 16 GB went into swap (1.3 GB used at idle
  after). Fixes: the main process holds no caches; `voice.PinnedPool`
  (single-process slots, each word's columns pinned to slots by
  McNaughton's wrap-around rule on columns × residual length): 1–3 words
  a slot, load within 2% of even (6 and 8 slots).
- Owner (2026-10-26): an 8-vCPU, 16 GB pod (EPYC 9575F; CPU may be
  pushed to full). Jobs ship as one tarball (14 s); each job builds the
  engine, peel and probe for the host; NATIVE_PEEL gains nothing on x86
  (numba 0.193 against native 0.204 s a Delta residual on the Threadripper:
  the CPU, not the port).
- Everything free from 20% (23 midpoints, 12 traverse, 63 posture; 8
  workers, 8 iterations a rung): 216 s, peak 7.9 GB. −19.7 dB. Postures
  0.17–0.39 → 0.15–0.49 cm² (IY, AH better; AA, AE worse). Move times:
  blade 51.5 (50), lips 66.6 (70), folds 61.7 (60), front 71.9 (80),
  back 121.7 (100), pharynx 79.1 (120). Raw lags 3–7 ms high together;
  Alfa's events 38 ms off. The gauge (all lags +c, all midpoints −c) was
  free and the solve drifted along it.
- DECIDED: the gauge fixed in `shapes.py` (the pharynx's lag held at zero,
  `GAUGE`): a convention, not a fitted value. Rerun with 24 iterations a
  rung (the Mini's run from the hints, 24 iterations, traverse known: 304
  s, tubes 0.26–1.18 cm², −19.0 dB; the 1 ms rung still descending when it
  ran out, the 0.1 ms rung found no step).

## Step 64: a posture needs contexts. MEASURED (pod, gauge fixed, 24 iterations a rung).

Nine words, everything free from 20% of the way to the hints: 755 s, peak
8.2 GB, −20.4 dB. Tube error start → end, with the number of words that
contain the vowel:

| vowel | words | start | end |
|---|---|---|---|
| IY | 4 | 0.189 | 0.069 |
| AH | 3 | 0.197 | 0.130 |
| EH | 2 | 0.230 | 0.148 |
| OW | 1 | 0.389 | 0.160 |
| AA | 5 | 0.165 | 0.168 |
| UW | 2 | 0.184 | 0.239 |
| AE | 1 | 0.295 | 0.591 |

- AE, in Alfa alone, went the wrong way while Alfa's events went 38 ms off:
  inside one word the posture and the timing traded (steps 10–15's
  non-uniqueness). The vowel heard in most contexts bar AA came back best.
  Not monotone (AA, OW), but the two ends are clear.
- The gauge fixed, the common drift went; the lags still carry 2–3 ms of
  compensation (back +2.3, front +2.7, lips +8.0, folds +6.7; truth 0, 0,
  +6, +4). Move times: blade 50.9 (50), lips 67.7 (70), folds 61.7 (60),
  front 77.2 (80), back 118.4 (100), pharynx 91.5 (120).
- HYPOTHESIS: a posture is pinned by the number of transitions it is
  heard in, not by the solver. Test running: twelve words, the nine plus
  AE-P-IY, UW-T-AE, AE-K-AA (AE then in four words).

## Step 65: twelve words — contexts rescue a lone vowel, not convergence. MEASURED (pod).

Twelve words (the nine plus AE-P-IY, UW-T-AE, AE-K-AA), everything free
from 20%, gauge fixed, 24 iterations a rung: 867 s, peak 9.2 GB, −21.1 dB
(13 Jacobians, 55 all-word residuals). Tubes start → end (nine words' end):
IY 0.189 → 0.053 (0.069); AH 0.197 → 0.090 (0.130); UW 0.184 → 0.181
(0.239); AE 0.295 → 0.210 (0.591); OW 0.389 → 0.209 (0.160); EH 0.230 →
0.216 (0.148); AA 0.165 → 0.226 (0.168). Alfa's events 38 → 25 ms off.
Move times: pharynx 125.4 (120), back 87.6 (100), front 78.7 (80), blade
49.6 (50), lips 64.5 (70), folds 65.6 (60); lags back 0.02, front 1.20,
blade −6.27 (−8), lips 8.65 (6), folds 5.95 (4).

- The step 64 hypothesis holds for its sharp case: AE, heard once, failed;
  heard four times, it improved. It does not hold as a cure: every vowel
  stays 0.05–0.23 cm² off while the truth fits at −249 dB.
- In every shapes run the 0.1 ms and one-step rungs change nothing (the
  postures identical between them); they cost about a third of a run.
  The stall is at the coarse-to-fine handoff, ~−20 dB.
- NEXT: the landscape around the stalled point, along the posture and
  traverse directions (as step 47 did for timing), before any change.
