# Papa timbre study: preserve the intelligible baseline

The user accepted round 009 as sounding like Papa, while describing a mechanical
sound and a bit of fax-machine accent. This is one accepted intelligible word;
naturalness is not accepted. The previous failed rounds remain failed.

## Fixed baseline

- Engine source commit: `920aaa6` (local air exchanges, conservative tongue).
- Audio: `out/visual-round-009/speech.wav`.
- Audio SHA-256: `029323121e3e1980ec3e7fa24d6e4741a96a6eaf0a532c4f2b70cdc6938dab84`.
- Enhanced field: `out/visual-round-009/enhanced.npz`.
- Body: `bodies/word-audition.body`; speaking rate 0.8; phones P AA0 P AA1.
- 48 kHz stereo, 72480 frames, one fixed playback gain 882.4414167632915.

The full baseline and human-reference enhanced fields were inspected without
resampling, recoloring or removing any time/frequency cells. The human is a
different body and pitch; its spectrum is not a template every ridge must match.

## Observations before intervention

The complete accepted field contains a strong descending harmonic family in
the second vowel, along with comparatively level resonance bands. Source traces
show a requested pitch rise to about 135 Hz and descent to about 92 Hz. Detected
cycles follow this trajectory fairly closely during sustained voicing. The
pitch detector holds its last estimate through silence; these held values are
not evidence of continuing phonation. There is no demonstrated controller
runaway. Per-cycle tension corrections still deserve scrutiny independently.

A supplemental three-cycle trace at 0.987–1.013 s shows glottal flow returning
to zero with a conspicuous corner in its derivative, followed by lip oscillation
while the glottis is closed. The configured fold contact smoothing area is
3e-8 m². Dividing that area by the locally observed area-change speed near
closure gives approximately 8–13 microseconds (10th–90th percentile in the
0.95–1.15 s window), below the 20.83-microsecond audio step. This is a local
scale estimate, not a measured physical closure duration or proof of aliasing.
The fold tissue uses two internal substeps, but the resulting opening is
passed to the local airflow update once per audio sample. Thus tissue substeps
do not themselves resolve the airflow seal. The configured one-strip model
also closes its entire 14 mm fold length together. These observations motivate
inspecting contact geometry and its airflow coupling, rather than smoothing
the rendered output.

The synthesis implementation and accepted body have been preserved. Diagnostic
body copies change only the explicitly recorded controls below. Observation
points retain the coupled folds and airway; a source tap is not an isolated
free-running glottis. In particular, airway pressure can appear in its output.

## Controlled rounds

**010 — coupled source observation.** Freshly rendered with exactly the
baseline body, plan, noise and fixed playback gain. Its physical trace, probe
and executable are byte-identical to round 009. The full 684 × 567 enhanced
field was inspected. The descending comb and considerable upper texture are
already present at the coupled source. Some narrow ridges continue through
the lip-closed interval at this internal observation point, while that interval
is largely quiet at the baseline microphone. This is internal motion, not
proof of an audible intersyllabic whistle. Airway feedback remains connected,
so the image does not assign every source ridge to the folds themselves.

**011 — stochastic jet pressures disabled.** Fresh full microphone rendering
with only `aspiration_gain=0` and `frication_gain=0` added to a diagnostic copy
of the baseline body. Mean-flow jet resistance, anatomy, plan, pitch controls,
jitter and shimmer remain active. This tests the generated pressure noise,
not all randomness or every possible turbulent effect. Its trace still has
both voiced flow regions. Full 684 × 567 inspection shows the descending comb,
upper bands and ringing remain visually very similar to baseline at the same
gain and color reference. Disabling these sources does not remove the targeted
visual features. This is not a claim that the waveforms or perceived timbres
are identical.

**012 — near-seal aperture rounding.** A separate baseline body copy changes
only `folds.contact_softness` from 3e-8 to 3e-7 m². Baseline aspiration and
frication are restored. This is a sensitivity experiment on the existing
contact law, not a fitted body, a new default, or an output filter. It retains
the word plan, 46 detected cycles and full duration. Supplemental traces show
the source derivative corners becoming rounded while the large open-phase
flow pulses remain close to baseline. The complete enhanced field shows a clear
reduction in dense upper texture, while both vowel regions, the principal
lower harmonic structure and descending pitch remain. Narrow resonances are
still visible. This is stronger causal evidence for the contact mechanism
than simply finding a source trace corner: modifying that mechanism changes
the targeted visual feature, while removing generated pressure noise did not.
The user subsequently reported: “it sounds softer but it still has that tinny
quality like a transistor hum.” Softening is audible; the target timbre defect
persists. This overrides any inference that the upper-field reduction alone
resolved the perceived problem.

## What this establishes, and the next physical question

The visible upper texture has a demonstrated contact-related contribution.
The engine radiates pressure from the time derivative of flow. The one-strip
near-seal corner makes an abrupt pressure edge; repeated edges excite the
same airway repeatedly. Altering the seal upstream of the coupled air state
rounds that edge and changes the full emitted field. It is not an output EQ
or a waveform replacement.

This does not prove that all remaining narrow bands are erroneous, that the
airway is adequately damped, or that this is the sole perceptual cause of the
reported fax accent. The large pitch descent is a separate planned motion;
it remains in the contact experiment. Realistic contact needs a defensible
relationship between fold geometry, progressively changing opening, air
momentum and tissue pressure. The existing pressure force uses signed strip
openings while flow uses a rounded minimum/soft seal; that constitutive
consistency is an unresolved next question. Do not promote an arbitrary wider
rounding area merely because it darkens the upper field.

The user-facing study sample is fresh round 012. The accepted body, dictionary,
engine and speaking-rate baseline remain unchanged.

## Limits

The enhanced field is the preserved nonlinear display, not calibrated power.
Image intensity and broad-band fractions are not perceptual quality scores.
No change in naturalness is claimed from traces or images. No recognizer score
is used to define this timbre problem. The baseline engine is retained while
its source/airway/radiation contributions are studied.

## Reproduction and retained evidence

All three diagnostic directories retain their executable, body, dictionary,
source snapshot, native waveform, fixed-gain waveform, physical trace, full
enhanced arrays/images and commands/hashes in `round.json`. The runner now
accepts `--tap` and `--playback-gain`; it rejects a declared gain that would
clip rather than silently renormalizing. Experimental body differences are
fully present in each frozen `body.body`.

```
# Same baseline options in each invocation of tools/render_visual_round.py:
--text 'Papa.' --rate .8 --playback-gain 882.4414167632915
# 010: --tap source; unchanged word-audition.body
# 011: --tap microphones; baseline + aspiration_gain=0, frication_gain=0
# 012: --tap microphones; baseline with contact_softness=3e-7
```

`reports/papa_timbre_study.json` records the exact round hashes and checks.
`out/visual-round-009/timbre-traces.png` and `cycle-detail.png` are supplemental
trace plots. `out/visual-round-012/contact-detail.png` overlays the contact
experiment at the same absolute sample times, with no phase alignment. These
plots do not replace the complete enhanced field of each fresh recording.

## Continuation: hum within the voiced vowels

The user clarified that the transistor-like quality is part of the voiced
sounds, not a separate sound between or after them. The audible softening of
012 therefore does not establish that contact-edge smoothing addresses this
remaining defect. The investigation continues from that body as a controlled
reference, without promoting it to the default.

**013 — lip observation.** Freshly generated from the unchanged 012 body at
its original common gain. The physical trace is byte-identical to 012. The
complete enhanced field retains the same voiced harmonic pattern and narrow
bands. The nasal trace is effectively zero; the scaled throat derivative RMS
is about 0.86% of the lip derivative before distance/delay propagation. This
makes microphone mixing or a separate nose/throat source an unlikely main
origin of the observed pattern. It is not an auditory equivalence claim.

**014 — rejected boundary-loss control.** An isolated source tree at
`out/papa-loss-engine-study` starts from commit 1dc517c and adds the six-state
thermal/viscous relaxation terms previously tried in failed round 002. The
reason for this repeat is the newly intelligible, correctly sealed baseline;
the old failed speech could not establish the effect on this voiced-timbre
problem. The main local engine and body are unchanged. Complete patch, build,
frozen sources and recording are retained locally.

The full enhanced field shows a conspicuous closure stripe at 0.807 s and
little change in the main voiced pattern. At that instant lip flow drops from
29.5 to 0.055 cm³/s in one sample while the represented lip aperture remains
about 0.011 cm². The implementation linearly interpolates divergent boundary
loss coefficients across geometry blocks; its treatment of near-closing
geometry and stored loss state needs examination. This is a hypothesis about
the new transient, not a proven diagnosis of the original hum.

The original playback gain would have clipped this new transient and was
rejected. The complete waveform was instead given one fixed gain
334.762270875, and its enhanced color reference was scaled by the same gain
ratio relative to 012. This preserves the physical amplitude comparison;
the listening files themselves do not have the same gain. The full source
hash and exact gain are in the round receipt. No click was cut or filtered.
The local-engine integrity suite passed, yet the speech change was rejected.
It also rendered at about 4.2× realtime, below the required 10× threshold.
This experiment is not merged into the working engine.

**015 — distributed fold geometry.** Return to the unchanged local airway and
the 012 body, adding only `folds.strips=16`. This exposes an existing physical
representation through body configuration: a bowed longitudinal vibration
shape with a posterior opening, instead of one rectangular aperture closing
along its full length together. The same moving lower/upper layers drive all
sampled strips; there are no added oscillators or scripted voice cycles.
Its shape and posterior opening distribution are model assumptions, not a
measured fit to a person. Its adduction coordinate treatment also remains a
reduced model rather than a complete three-dimensional tissue description.

The source trace changes over the whole pulse: a small posterior airflow
persists instead of a long exactly shut interval. The inspected complete enhanced field retains both voiced regions and the
descending harmonic family. Relative lower harmonic strengths and upper
texture change, but similar narrow resonances remain; this is not a visual
demonstration that the transistor-like quality has disappeared. This diagnostic
render takes 0.156 s for 1.51 s of audio (about 9.7× realtime), slightly below
the performance requirement. Performance work remains parked while timbre is
assessed; this is not a performance-qualified replacement. Default bodies
retain the existing one-strip geometry.

**016 — fresh voiced-level-matched audition.** Re-render the exact 015 physical
configuration, then apply one constant gain to the entire recording, chosen
to match 012's stereo RMS over voiced windows 0.72–0.79 and 0.95–1.15 s.
The declared gain is 1203.4128531012666; the color reference scales in the same
proportion to 0.24285731813851033, preserving physical-level comparisons with
the previous fields. This is RMS level matching, not a perceptual loudness
model. Its purpose is to keep a lower vowel level from confounding a timbre
comparison. The 015 physical run is freshly repeated, not reused as an old
audio deliverable. The final round receipt verifies byte-identical native audio and physical
traces, the exact waveform-to-full-field hash, and a PCM voiced RMS ratio
of 1.00000036 to 012. The complete 684 × 567 field was inspected: both vowel
regions remain, along with residual narrow resonances. Its appearance agrees
with 015 at the common physical scale. The diagnostic render was about 9.9×
realtime with full sample-rate traces and probes enabled.

The user's description of 012 as softer alone does not distinguish loudness
from timbre: the earlier common physical gain preserved actual level changes.
No naturalness acceptance is inferred for 015 or 016 from visual changes.


The retained repository change is only the transactional `folds.strips` body
configuration input, its documented geometry semantics and parser tests. The
existing fold/airway equations and all default bodies remain unchanged. The
updated stream test executable passes with zero failures, including invalid
integer rejection and existing chunk-invariance checks. The failed loss copy
remains isolated and is not a candidate default. A CTest regex invocation
selected no tests; validation used the actual test executable directly.


## Source investigation after the 016 rejection

User feedback on 016: “sounds worse, actually. the p is less distinct. and the
underlying behavior i mentioned remains.” The user suspects the source
excitation. Round 016 and its identical physical configuration 015 are rejected
for advancement: weaker P articulation and no relief of the embedded voiced
hum, despite the RMS-matched comparison. Return to the one-strip 012 body for
investigation. The accepted intelligibility count remains one word (009);
there is no accepted natural voice.

Source audit: the current finite-flow force uses the pressure difference
`lung - tract` as its lower-layer surface pressure before the kinetic-pressure
subtraction, and the upper-layer aerodynamic drive is zero. Those are different
quantities: a drop drives flow, while local pressure relative to the tissue
exterior drives its surface. The existing simplified model effectively takes
the downstream pressure as the surface-pressure reference. This is a concrete
omission to test, not an established cause of the listener's hum.

Relevant primary research: [Titze 2022, How can vocal folds oscillate with a
limited mucosal wave?](https://pmc.ncbi.nlm.nih.gov/articles/PMC9533394/) distinguishes
transglottal pressure from medial-surface pressure and describes the alternating
push from an inertive supraglottal load. It motivates checking this reference;
it does not validate our reduced pressure distribution or predict this word's
quality. [Perrine and Scherer 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10497319/)
compares low-dimensional source flow shapes with human recordings, including
rise/fall skew and closing slope. Our microphone recording alone does not give
an independently measured human glottal waveform.

**017 — instrumented one-strip baseline.** Fresh audio is byte-identical to
012. Added observation of signed lower/upper openings, rest opening, lung
pressure, and endpoint glottal flow; these do not alter dynamics. The complete
684 x 567 enhanced field was inspected before the next source experiment: two
voiced regions remain with a quiet closure, strong lower harmonics and residual
upper horizontal ridges. Supplemental whole-utterance and cycle-detail traces
show nearly symmetric source pulses, smooth mean flow, and closely agreeing
mean and endpoint flow. No large sample-to-sample source chatter is evident
in the inspected cycles. This does not exclude a subtle numerical defect.

**018 — downstream pressure on exposed fold tissue.** The isolated current
engine copy `out/papa-source-pressure-study` applies downstream gauge pressure
to exposed lower and upper fold surfaces. The original separated-jet drop
calculation remains; pressure is added only where the signed layer gap is
open. The new upper force uses the existing upper thickness and fold length.
No spring, contact, phonetic, lip, tract, timing or body parameter changes.
Default arguments keep isolated calibration and earlier callers unchanged.
This is a limited source-force experiment, not a complete description of
intraglottal pressure or a fluid/tissue energy-conservation proof.

The patch is retained in `source_pressure_experiment.patch`; it applies to the
instrumented main tree with `git apply`. Candidate equations are isolated
while the audible outcome is assessed. Both main and candidate local-engine
test executables pass. The candidate includes three pressure controls: equal
positive lung/tract pressure pushes exposed tissue despite zero pressure drop;
negative gauge pressure pulls exposed tissue inward; contact shields a closed
lower surface while downstream pressure acts on the open upper surface. These
tests establish the intended force behavior, not timbre.

018 is freshly rendered with the exact 012 body and speaking rate. Its one
constant playback gain is 571.3160472874877, selected from the same voiced
windows as 016. The final PCM RMS ratio to 012 in those windows is
1.0000002276. The entire sound has that same gain; no cuts, equalization,
compression or waveform replacement. The enhanced display reference is
0.11529566324322615 to preserve the common physical scale.

The source-cycle plot shows a steeper closing slope and a changed pulse shape;
52 detected cycles versus baseline 46 also means the loaded oscillator changed,
not just its output gain. Both release-window RMS levels remain within about
2 percent of 012 at the matched vowel level (0.690–0.707 and 0.890–0.909 s);
these windows are waveform checks and do not prove that either P is perceptually
clear. The traced render takes 0.156 s for 1.51 s of audio, about 9.7x realtime,
and is therefore not performance-qualified against the 10x target. Performance
work remains parked per the user's later instruction.

The complete 018 enhanced field was then inspected: both voiced regions and
the separating closure remain; no new conspicuous broadband click stripe
appears. Mid/upper voiced harmonic energy increases and lower harmonic
strengths change at the common physical display scale. Narrow bands remain.
This is evidence of a changed excitation and preserved visible word structure,
not evidence that the transistor quality has disappeared. Listening approval
is pending. The 017 full enhanced array is exactly equal to 012, confirming
that instrumentation did not change the baseline. No image was cropped,
desaturated or reduced; the full 684 x 567 cell images belong to their freshly
rendered WAVs and were inspected internally.
