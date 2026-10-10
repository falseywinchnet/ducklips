# Continual imperfect holding throughout the speech engine

This implements the user's request for one mechanism that maintains the body
against passive relaxation, irregular overlapping muscle twitches and delayed
correction. It is integrated into `Synthesizer` under `living_control`, selected
in the CLI with `--living-control`. This path uses the current elastic-volume
`TissueSource`, the moving conservative tongue and the local airflow engine.
It does not run the earlier two-mass source. Human naturalness is unaccepted.

## Listening feedback, 10 October

The listener clarified that the reported silence was their muted playback,
then described the latest take as: "it sounds like someone sickly."
This applies to the round 003 waveform, freshly reproduced for the mono MP3
in `out/somatic-playback-004/`; the fresh raw waveform was verified identical.
This is evidence of a person-like but unhealthy-sounding impression, not an
acceptance of healthy natural phonation. The specific distinction between
hoarseness, weakness and shakiness has not yet been supplied.

Round 004 is an isolation experiment, generated with
`ducklips_somatic_study out/somatic-round-004 --no-control-variation`.
It sets only `SomaticControl::Settings::noise` to zero. Delayed compensation,
mean twitch dynamics, opposing tone, passive relaxation and finite breath
remain active. The tissue source's separately varying relaxation also remains
active, so this is not an all-randomness-disabled comparison. The normal
engine retains its requested imperfect holding mechanism. The exact render
receipt and complete enhanced inspection are in `somatic_control_rounds.json`.

The complete native enhanced field for round 004 still has broad sustained
texture and weak harmonic organization. The late vertical transients of
round 003 are absent. Its fold-area and flow traces retain large irregular
cycle-to-cycle changes despite much smoother holding trajectories. Thus the
added holding variation is not necessary for the sustained visual defect;
the post-release transients are affected by it. This does not yet establish
which physical assumption produces the perceived sickly quality. In
particular, tissue relaxation variability and the approximate source/tract
pressure coupling remain separate candidates. The listener has not assessed
round 004. Normal engine settings and their imperfect holding remain intact.

## What carries history

`SomaticControl` maintains 29 coordinates: tongue body, jaw, blade and root;
lip opening and protrusion; three constriction gestures with five coordinates
each; tongue drag; velum; glottal posture; breath pressure; pitch intent; and
effort. Fixed anatomical residuals and microphone gain are not muscle states.

Each coordinate has position, velocity, relaxing strain memory, two opposing
motor populations, overlapping twitch state, filtered sensing, delayed sensing
and bounded integral adaptation. A cue changes intention. Existing physical
and neural history continues across cues and audio callback boundaries.

The activation of a selected constriction is distinct from resting muscle
tone. An unused contact slot does not receive tonic recruitment: the tongue,
lips and jaw underneath it already have their own opposing tone. When a new
idle gesture acquires a contact location, its coordinate reference translates
while its small displacement, velocity, stress and neural history persist.
This prevents dormant, undefined contact locations from becoming active
because of a holding twitch.

The reduced mechanical equation in normalized coordinates is

```
q'' + 2 r q' + r²(q - g m) = r²(muscle + load)
m' = positive_rate(t) * (q - m)
```

Here `m` is stress memory and `g` is the relaxing fraction. At a held position
the relaxing part of restoring stress decreases. Neural drive supplies a
steady-load estimate plus delayed proportional and integral corrections.
Two populations of 64 irregularly firing motor units maintain opposing
effort. The population mean is carried by a smooth twitch response; a small
fraction of the finite population's deviation from its expected mean remains.
Randomness enters motor recruitment/firing and positive relaxation rates;
there is no random position offset, waveform multiplier or prescribed vocal
cycle. Related controls share a small effort variation and retain independent
motor histories. Random seeds permit exact replay.

Prototype defaults are 15% relaxing stiffness, 200 ms strain-memory time,
20 ms sensory filtering, 30 ms transport delay, 2% opposing tone, and 8% of
finite-pool fluctuations. The reused motor units have 15/35 ms twitch times
and irregular firing intervals. These are modeling choices, not an identified
human neuromuscular system. The same twitch distribution across different
muscles is a current simplification. The controller can compensate for small
loads; these checks do not establish anatomical correctness.

Breath has finite volume and compliance, depleted by measured glottal outflow.
Airway pressure applies small assumed generalized loads to tongue, jaw, lips
and velum. The four physiological coordinates are interpolated within each
1 ms control block. Old jitter, shimmer and periodic tremor are disabled in
this path. Glottal posture changes tissue equilibrium without teleporting
surface positions; pitch intent changes the compiled geometric prestress.

## Physical limits

- Most controls are reduced posture coordinates, not individual anatomical
  muscles. Their normalized loads and capacities are not forces in newtons.
- The source's backing/ends remain fixed. Laryngeal control senses its posture
  coordinate, not measured cartilage position or the tissue's cycle-averaged
  gap. The pitch-intent mapping is a prestress hypothesis, not an F0 servo or
  an established pitch calibration.
- The source computes quasi-steady wall pressure while the tract inlet carries
  airflow inertia. Their coupling is approximate. The coupled output shows
  much more broadband structure than the user's sustained-vowel reference;
  the mechanism causing the excess has not been isolated in this experiment.
- Respiratory compliance, coordinate scales, pressure-load projections and
  feedback gains need body-specific identification.
- Continuous neural maintenance does not establish human-sounding speech.

## Reproduction

Build `ducklips_somatic_study`, run it with a new output directory, then run
`tools/package_somatic_study.py` on that directory using the scientific Python
environment. The harness records a fresh 1.2 s /ah/ intention with pressure
and opening release at 0.75 s, all controller states and sample-level physical
traces. Packaging preserves a source/executable snapshot and hashes, applies
one constant listening gain, and generates the complete enhanced field.
There is no output fade, EQ, crop or added audio noise. Only the original
magma field, at all 684 frequency rows and every time cell, is inspected.

The full engine requires 48 kHz and 48-frame control blocks on this path.
This restriction is validated before the real-time reset/render functions.
The public controller supports whole-millisecond updates and configurable
noise, tone, relaxation, feedback, breath compliance and available volume.

## Verification and rounds

The focused checks cover sustained-load compensation, identical-seed replay,
independent motor histories, neutral-posture opposing tone, finite breath
accounting and fixed anatomy. With noise disabled, feedback reduced the
remaining tongue displacement under a fixed load from 0.188235 to 0.00641203
after four seconds. With finite-pool variation enabled, a fixed 120 Hz intent
produced a holding-coordinate range of 119.280–120.924 Hz over two seconds;
that is an internal intention coordinate, not measured acoustic F0.

The integrated test includes a new dorsal-contact cue and compares contiguous
rendering with 73-frame chunks that cross the 48-frame control boundaries.
The resulting audio is identical. Invalid settings are rejected before
entering the noexcept reset. Existing core synthesis, local airflow and tissue
source checks also run as focused regressions.

- Round 001: complete integrated audio, no tonic activation at neutral intent.
  Enhanced field inspected: extensive broadband texture, weak low harmonic
  organization, a clean physical release; not accepted as a human vowel.
- Round 002: opposing tone also recruited unused gesture slots. Geometry
  failed after 3552 frames because an undefined dorsal contact became active.
  Its full enhanced field confirms the early stop. This is a failed render;
  its apparent 57.9x rate is meaningless because the rest is failure silence.
- Round 003: tonic recruitment respects gesture activation and idle contact
  references are relocated without erasing physical history. All 57600 frames
  render without geometry/airflow failure, clipping or nonfinite samples.
  Full enhanced inspection still finds broad texture rather than clear human
  harmonic organization, and narrow late transients after pressure release.
  Their cause is not yet isolated. The exact receipt and visual findings are
  recorded in `somatic_control_rounds.json`. This is not a naturalness success.

Each round resides in `out/somatic-round-NNN/`. The newest full-engine case
uses the tissue source at 96 kHz and a 24-section oral / 9-section nasal tract
at 48 kHz. It is a diagnostic of the mechanism, not an accepted speech result.

## Efficiency boundary

The first complete integrated render took 0.217 s for 1.2 s of sound: 5.52x
realtime. Instrumentation reproduced exactly the same raw audio and measured
0.0062 s in the shared controller, 0.0909 s in tissue, 0.0917 s in airflow and
0.0209 s in geometry. Controller complexity is not the main bottleneck.

The representation to investigate below the 10x requirement is continuous
slowly varying geometry/force coefficients coupled to the existing acoustic
state, with justified multirate tissue updates. Recomputing tissue and airway
coefficients at every acoustic substep dominates this case. That change is
not yet implemented or validated; it must preserve contact, pressure feedback
and the enhanced structure of the resulting sound. Performance optimization
remains secondary to sound quality by the user's instruction.

Round 003 measures 0.21837 s for 1.2 s, or 5.495x realtime: 0.00626 s control,
0.09384 s tissue, 0.08731 s airflow and 0.02005 s geometry. Render timings
include diagnostic trace capture and timers but exclude CSV/WAV writes and
offline enhanced depiction. These are local observations, not a general
performance guarantee.
