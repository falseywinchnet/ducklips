# VocalTractLab 2.4 — gestural control (study notes)

Source: VTL 2.4 API `Sources/` (GPL; read for the mathematics only), mainly
`GesturalScore.cpp`. Line numbers refer to that file unless named. OBSERVED =
read in code; DERIVED = arithmetic from its constants; CANDIDATE = our proposal.

## Tiers and gestures

Gesture `{duration_s, dVal, slope, sVal, tau_s, neutral}`; target
`{duration, value, slope, tau_s}` (`GesturalScore.h` 44–68). Tiers with clamped
ranges (`GesturalScore.h` 152–163, 707–734; `limitGestureParams` 612–661):

| Tier | Value | Slope | τ range |
|---|---|---|---|
| VOWEL | named shape | – | 8–20 ms |
| LIP, TONGUE_TIP, TONGUE_BODY | named consonant shape | – | 5–40 ms |
| VELIC | VO −0.1…1.5 cm² | – | 5–40 ms |
| GLOTTAL_SHAPE | named glottis shape | – | 10–40 ms |
| F0 | 48…108 st re 1 Hz | ±200 st/s (normally ±80) | 5–40 ms |
| PRESSURE | 0…16000 dPa | – | 5–40 ms |

Each tier partitions time (`putGesture` 263–407). Default τ = 10 ms "for 100 ms
transition durations".

## Target approximation (calcParamCurve 5192–5315)

Per target j (value b, slope m, τ), on local time t ∈ [0, d_j]:

    x(t) = (b + m t) + e^{a t}(c0 + c1 t + c2 t² + c3 t³ + c4 t⁴),  a = −1/τ

A 5th-order critically damped system (five poles at −1/τ) tracking a linear
target. At each boundary the value and derivatives f0…f4 are carried over
(slope added to f1) and the next coefficients solved (5299–5303):

    c0 = f0 − b'
    c1 = f1 − c0 a − m'
    c2 = (f2 − c0 a² − 2 c1 a)/2
    c3 = (f3 − c0 a³ − 3 c1 a² − 6 c2 a)/6
    c4 = (f4 − c0 a⁴ − 4 c1 a³ − 12 c2 a² − 24 c3 a)/24

C⁴-continuous. If τ increases, f_k *= (τ_old/τ_new)^k (5273–5285) to avoid
overshoot. First target starts on its line (c = 0); last extends.
Step response = Gamma(5) CDF, `1 − e^{−t/τ} Σ_{k=0..4}(t/τ)^k/k!` (5326–5379).
DERIVED quantiles: 10% 2.43τ, 50% 4.67τ, 75% 6.27τ, 90% 7.99τ, 97% 9.96τ;
ramp lag 5τ. Curves sampled at 400 Hz; tube areas interpolated between frames.
F0 filtered in semitones. Direction-dependent τ factors exist (5048–5081),
default 1.

## From tiers to per-parameter targets (calcTractParamTargets 4707–5083)

1. Cut the five tract tiers at the union of boundaries into slices.
2. A neutral consonant/velic gesture right after an active one is a **RET**
   (return) phase.
3. Start each slice from the vowel shape (fallback schwa, then neutral),
   τ = vowel τ.
4. Consonants override (below). In a RET slice values revert to the vowel but
   τ of all 19 parameters = the return gesture's τ.
5. Return gestures after every closing gesture: duration min(10·offsetTau,
   following neutral run), τ = offsetTau (`putClosingGesture` 3043–3101), so
   release speed is independent of approach speed.

## Consonant context interpolation

Map the current vowel x onto the a–i–u plane
`x ≈ x_a + α(x_i − x_a) + β(x_u − x_a)` (`mapToVowelSubspace` 4159–4281):
tongue by least squares over 17 parameters (pseudo-inverse (MᵀM)⁻¹Mᵀ,
4093–4146); lips by an exact 2×2 solve on (LP, LD). Clamp 0 ≤ α, β ≤ 1; if
α + β > 1 subtract half the excess from each (4290–4353).
Consonant `C = C(a) + α_T(C(i) − C(a)) + β_T(C(u) − C(a))`; LP, LD with α_L, β_L
(`getContextDependentConsonant` 4364–4408).

Dominance by place (4921–5009): with any consonant active, every target starts
as the mean of active consonant shapes; then the tongue-body consonant owns
TCX TCY TBX TBY TTX TTY JX JA LP LD TS2 TS3; tongue-tip overwrites TBX TBY TTX
TTY JX JA LP LD TS3; lips overwrite JX JA LP LD. The most anterior consonant's
τ governs all parameters. Anticipation comes from the vowel tier leading.

## Automatic score from phones, durations, F0 (createFromSegmentSequence 854–2480)

- Affricates split into stop + fricative (stop share pf 0.4, ts 0.33, tS 0.33,
  dZ 0.4). Pauses ≥ 100 ms at start, 100 ms at end. /d t n l/ become
  postalveolar near /S Z/.
- **Lung pressure**: background 0 (τ 20 ms); each pause-free stretch 8000 dPa
  (800 Pa), τ 10 ms. Onset lead 40 ms (default, /h/), 60 ms before a fricative,
  plosive midpoint − 100 ms before a voiceless plosive. Offset solved so the
  20 ms decay reaches 200 Pa (voiced), 150 (plosive) or 100 (fricative) at the
  last segment's end (leads ≈ 125/137/152 ms, DERIVED).
- **Vowels**: gesture [t1, segment end], τ = 0.1(t2 − t1) clamped 8–20 ms.
  After a pause t1 = begin − 100 ms, t2 = begin. Hiatus: begin ∓ 50 ms. Across
  consonants (Kleiner et al. 2015), C = consonant stretch, V1 = previous vowel:
  `t1 = V1end − 0.035 + 0.208 C − 0.441 V1dur`,
  `t2 = V2begin + 0.076 − 0.462 C`, then t1 ≥ max(V1end − 50 ms, V1 midpoint,
  previous movement end, 0); t2 ≤ V2begin + 10 ms; t2 − t1 ≤ 200 ms.
- **Diphthongs** (Fietkau et al. 2022): transition D = 0.080 + 0.317·dur,
  centred; V1 gesture τ ≤ 10 ms, V2 τ = 0.1 D.
- **Fricatives**: oral gesture [begin − 60 ms, end − 50 ms], τ 10 ms, then a
  100 ms return τ 10 ms (+40 ms before a plosive of another articulator).
  Glottal gesture voiceless-/voiced-fricative over [begin − 50, end − 50] ms
  (+40 ms if voiceless before p t k). /h/: glottal "h" [begin − 60, end − 60].
  /l/: [begin − 40, end − 40], τ 8 ms, 80 ms return.
- **Plosives/nasals** placed last. VOT (s), x = segment duration capped at
  150 ms (voiceless) / 100 ms (voiced): p 1.577x² + 0.215x; t 1.641x² + 0.343x;
  k 1.338x² + 0.387x; d −2.430x² + 0.483x; g −1.867x² + 0.507x; b 0 with modal
  glottis. Onset/offset τ (Kohler 1980): voiceless 8/14 ms, voiced 10/6,
  nasals 12/12. Velic opening m 0.3 cm² (±20 ms), n 0.2 (±30), N 0.3 (±10);
  opening τ 12 ms, closing 10 ms. Closure [begin, end − VOT] ≥ 10 ms. Plosive
  glottal gesture from closure − 25 ms to end − 65 ms (voiceless) or end − 80 ms
  (voiced), with stated adjustments before voiced fricatives, homorganic
  nasals and /l/.
- **F0**: nodes at 20/50/80% of each phone; one gesture per node pair with
  value = node (st), slope = Δst/Δt, τ = 5 ms (≈25 ms lag); first gesture
  shortened 10 ms. Glottal stop dips 1 st per 10 ms.
- **Closure timing by search** (`addClosingGesture` 2768–2933): move the
  gesture start/end in 2.5 ms steps (and shorten offset τ to 5 ms) against the
  forward geometry until the tract is closed exactly over the intended
  interval. Closed = section area < 2·1e-4 cm² formed by the right articulator.
- Aspiration = glottal opening outlasting the oral release; bursts and
  frication come from the acoustics, not scored events.

## CANDIDATE ideas for our planner
1. Nth-order critically damped target approximation (closed form, state
   continuity) instead of second-order motors; τ = move/10.
2. Targets with slopes for F0 (semitones) and pressure; compensate ≈ 5τ lag.
3. Separate approach and release τ (return gestures), per class.
4. Three context prototypes per consonant + place dominance, replacing
   pairwise coarticulation rules.
5. Plan acoustic landmarks (closure, velic, VOT) and solve gesture timing
   against the forward model.
6. Use the empirical timing laws (Kleiner, Fietkau, VOT polynomials, glottal
   leads, pressure ramps) as priors/features for the learned models.
