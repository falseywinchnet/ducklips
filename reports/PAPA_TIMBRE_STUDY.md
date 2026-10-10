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
The recording has not received a human naturalness verdict.

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
