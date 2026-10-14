# DuckLips design

Labels follow the repository's epistemic vocabulary. Nothing here is
**DECIDED**; the mechanism is a **CANDIDATE** that runs.

## 1. The body is a branched acoustic network

The airway is a chain of oral sections from the glottis to the lips (64 by
default) with a nasal chain (24 sections) joined at the velum. Each section is
a compressible air volume, a node with compliance `C_i = A_i l_i / (rho c^2)`.
Each passage between sections is an air mass, a branch with inertance
`M_j = rho (l_a/2A_a + l_b/2A_b)`. The lips and nostrils are branches to still
air through a radiation end correction. Every node also has a yielding-wall
branch to ground: mass and resistance per unit wall area (Ishizaka-French-
Flanagan cheek values). The velar port is a branch whose conductance `1/M` is
proportional to its opening, so closing the velum simply removes the
coupling.

With pressures `p` on nodes and flows `u` on branches, and incidence `G`,

    C dp/dt = -G u + s_node        M du/dt = G^T p + s_series

and the stiffness operator is `K = G M^-1 G^T`. **OBSERVED:** the tract,
port and nasal branches form a tree; radiation and wall branches only add to
the diagonal. In the symmetric form `S = C^-1/2 K C^-1/2` the off-diagonal
pattern is that tree.

Closure is well-behaved in this form: a closing passage drives its
conductance `A/(rho l)` to zero; nothing becomes infinite.

## 2. Eigenstates by sparse tree solves

The architect's request was to "moose" the geometry: obtain eigenvalues so
the wave calculation is a sparse solve that can be deformed anywhere.

- **Tree elimination.** Eliminating nodes leaves-first gives a shifted
  `S - sigma I = L D L^T` with no fill, in O(n). This is the wrench transport
  solver's elimination of bodies over the contact graph, applied to air.
- **Inertia.** By Sylvester's law the number of negative pivots equals the
  number of eigenvalues below `sigma`. One factorization therefore says how
  many modes lie below 12 kHz and verifies the index of any candidate.
- **Deformation.** Each control step starts every mode from its previous
  shape. The first Rayleigh quotient of the old vector under the new matrix
  is exactly the first-order eigenvalue perturbation, i.e. the classic
  formant sensitivity function (kinetic minus potential energy density).
  Rayleigh quotient iteration refines it, the inertia count verifies its
  index, and an overlap test rejects a vector that duplicates an accepted
  mode. Failures fall back to inertia bisection plus deflated inverse
  iteration.
- **Sensitivities.** `A_i d(lambda_k)/dA_i` is computed analytically from the
  mode shape. It drives vowel inversion (section 5).

Pressure shapes are C-orthonormal, so with `z_k = a_k + i b_k` the acoustic
energy is `(1/2) sum |z_k|^2` and each mode obeys

    dz_k/dt = (-alpha_k + i w_k) z_k + phi_k^T s_node + i psi_k^T s_series,
    psi_k(j) = g_j (phi_k(from) - phi_k(to)) / w_k.

That is integrated exactly per sample. `alpha_k` collects viscous, thermal,
wall and radiation losses evaluated at the mode frequency, so losses are
frequency-dependent without extra filters.

**Moving from one geometry to the next.** The physical state, node pressures
and branch momenta, is reconstructed and projected onto the new modes.
"Passive" reprojection keeps `min(C_old, C_new) p` per node and
`min(g_old, g_new) M_old u` per branch. By Bessel's inequality in the new
C- and M-metrics this cannot increase energy. Modes may appear, vanish or
cross without being tracked. When the area function changes fast, the
control block is split into up to 8 interpolated geometry steps.

## 3. Sources: membranes, air and turbulence

- **Vocal folds.** A two-mass model (Ishizaka-Flanagan parameters in the
  Steinecke-Herzel form): two masses per fold, springs, a coupling spring,
  collision stiffness, and Bernoulli pressure on the lower mass. Tension Q
  scales stiffness up and mass down; pitch comes from a feedforward map
  (calibrated at prepare) plus a laryngeal feedback loop on the measured cycle
  period. The glottal flow is solved implicitly against the tract's own
  instantaneous response (a closed-form quadratic), so source-tract
  interaction is part of the physics.
- **Flow split.** The modes carry the fluctuation of the glottal flow about
  its two-pole 20 Hz mean. A lumped intraoral cavity (air plus
  static wall compliance) carries mass: it fills with the true glottal flow
  and vents through the oral constriction and the nasal path with Bernoulli
  orifices. The implicit step reduces to a quadratic in `sqrt|P|`. Its
  smoothed outflow derivative radiates as the low-frequency release pop.
- **Turbulence.** A jet beyond the critical Reynolds number makes a series
  pressure source `K * 0.5 rho v^2 * gate(Re) * sqrt(A/0.1 cm^2) * noise`. At
  the glottis this is aspiration, which also buffets the lower fold surface
  and is what starts the folds. At the oral constriction, or at the teeth for
  sibilants, it is frication. Gains are calibrated, not measured.
- **Perturbations.** Pitch jitter (low-passed noise on tension), tremor, and
  lung-pressure shimmer.

## 4. Air: the photonic transport discipline, made coherent

Radiators are the lips, nostrils and throat walls; their volume-flow
derivatives are emitters. `AirPath` compiles sparse links once for fixed
geometry, from each emitter (and its room images, ordered by order and cut
off by gain) to each of two microphones. A link carries a delay, a pressure
gain `rho / (4 pi d)` and a per-link one-pole wall loss. This follows the
photonic transport field: emitters, sparse links with per-link extinction,
receivers, unit responses compiled once and scaled by the live source, and a
march truncated below a cutoff.

The deliberate difference is that sound at a microphone is coherent, so
links carry delay and sign where light transport carries energy fractions.
A diffuse late tail, where an energy treatment would be valid, is not built
(**CANDIDATE**).

## 5. Behaviour: from text to cues

- **Dictionary.** A hand-entered core lexicon (about 250 words), an optional
  CMU-format file, and letter-to-sound rules as a last resort, which are
  reported.
- **Vowels.** A vowel is a posture: six tapered cosine modes of log area
  between the epilarynx and the lips, plus lip opening and protrusion. Each
  voice inverts formant targets (recalled Hillenbrand et al. 1995 male
  averages, scaled by tract length) into postures. The method is
  Levenberg-Marquardt with Tikhonov regularization toward neutral, using the
  analytic eigen-sensitivities. These are the "formant information gleaned
  from a language dictionary".
- **Consonants.** Gestures in three articulator slots (lips, tongue tip,
  tongue body), with place, width, target area and an optional teeth
  obstacle. A slot glides when re-targeted while active. Contact follows a
  smoothstep of activation, so closure meets zero area with zero rate.
- **Timing and coarticulation.** Duration rules cover stress, tenseness,
  voicing of the next consonant, phrase-final lengthening, cluster shortening
  and small deterministic variation. Gestures overlap, and vowel-to-vowel
  coarticulation comes from the tongue heading for the next vowel partway
  through the consonants. Critically damped articulator dynamics do the rest.
- **Glottis, velum and lungs.** Voiceless consonants abduct the glottis, with
  a voice-onset lag after stressed stops. Nasals open the velum with spill
  into neighbours. Each breath group inhales (audible when the voice breathes
  audibly), raises lung pressure, and lets it fall at the end.
- **Intonation.** Declination, accents on stressed vowels with downstep, and
  a final fall, or a rise on questions.
- **Streaming.** Cues are streamed through a lock-free SPSC queue a quarter
  second ahead of the render position.

## 6. Voices are bodies

| Voice | Tract | Folds | Habits |
|---|---|---|---|
| baritone | 17.5 cm | reference two-mass | 110 Hz, 800 Pa |
| alto | 14.8 cm | mass x0.45 | 205 Hz, breathier |
| child | 12.2 cm | mass x0.25 | 265 Hz, slower |
| gravel | 18.2 cm | pressed adduction | 88 Hz, high jitter, shimmer and tremor |
| whisper | 16.0 cm | folds held still | turbulence only |

## 7. Real-time contract

`Synthesizer::render(std::span<float>)` fills interleaved L,R 48 kHz float
frames. It does not allocate, lock or throw (all storage is sized in
`prepare`), which is the `gui_forms::AudioGenerator` contract.
`ducklips::GuiFormsVoice` adapts it.
