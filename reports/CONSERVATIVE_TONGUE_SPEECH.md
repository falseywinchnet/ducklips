# Conservative tongue in the speech engine

The volume-conserving material map is now connected to the canonical speech
path. It changes the actual airway used by the vocal-fold and acoustic solve,
including the mapped blade attachment. These recordings use
`--canonical-motion --conservative-tongue`; they are not the earlier compact
modal engine or a rendering of an isolated geometry demonstration.

Performance work remains parked. This is a functional and physical integration,
not a claim that speech is perfected. No NATO word has yet passed a documented
human pronunciation assessment.

## Geometry that reaches the acoustic solve

`ConservativeTongueSurface` constructs a dorsal lift from the current posture,
body gesture, native wall and commanded clearance. Unlike the earlier isolated
fixture's 0.9 cm contact patch and 1.5 cm shoulders, this integration uses the
existing tract model's 0.3 cm contact core and the incoming gesture width as
shoulder length. Shoulders are shortened to stay inside the resting body.
These are explicit geometric assumptions, not fitted human muscle data.

The existing divergence-free material flow moves the body. Its lower anchor is
`centre_y - radius/2`; the free gate starts at the lowest resting height of the
contact core. The command uses smoothstep activation, capped at contact. The
native map retains 32 implicit-midpoint substeps. Its conservative identity and
separate material validation are in [TONGUE_FLOW.md](TONGUE_FLOW.md).

The blade root is the image of its original body material point. Subsequent
allowed blade drag translates the same body and attachment. Airway rays are
intersected with the mapped circle itself: a 128-point contour brackets roots,
then Newton/bisection refines the actual surface. An inverse map detects wall
points inside the body. A failed map or ray solve stops synthesis and records
an error; no rigid-body fallback supplies replacement speech. At zero body
activation the resting circular body remains the common undeformed limit.
The existing vowel fit and empirical gap-to-area relation remain in use.

## A contact failure found during integration

The first connected alphabet render clipped 22 channel samples (11 stereo
frames). Geometry inspection also exposed an independent physical violation:
the fixed-length blade could pull the already closed dorsum through the wall.
In the combined closure fixture, full tip activation produced 0.88575 cm of
body/wall overlap. Passing the material solver was not sufficient evidence of
valid contact.

Blade drag is now limited along its requested direction by sampled body/wall
contact. The check visits native wall knots along the translation before
bisecting the first invalid interval. If the tip still lies beyond the fixed
2.2 cm blade reach, it is projected to the reachable position and its unmet
travel is reported. This preserves attachment and body volume. It does not
pretend the requested consonant was achieved.

The fixture now has sampled overlap 9.9998e-11 cm and a reported 1.22610 cm
unmet tip target. The contact-constrained alphabet before the spelling fix had
zero clipping, peak 0.690 without normalization, and maximum unmet tip travel
1.37333 cm. No gain change was used to remove clipping. This association does
not establish a complete acoustic cause for every previously clipped frame.

These contact checks sample 128 contour points. They are not a continuous
collision certificate, do not cover every articulator, and do not establish
physical validity across arbitrary bodies. Reported overlap concerns the
mapped body while its dorsal gesture is active. The fixed-length blade and
its unmet targets are major remaining model limitations.

## Word planning

The literal NATO spelling `X-ray` previously split into `X` and `ray`, bypassed
`XRAY` in the explicit dictionary and lost the opening EH vowel. The dictionary
now includes `X-RAY`; the planner preserves dictionary-recognized hyphenated
compounds as one word. Unknown compounds and repeated hyphens retain word
boundaries. This fixes the phoneme/nudge input, not the acoustic waveform.
All 26 written codewords now use explicit dictionary pronunciations.

## Current audio and enhanced fields

The final Quebec sample has 66,720 frames at 48 kHz (1.39 s), zero geometry
failures, zero nonfinite or clipped samples, zero sampled mapped-body overlap,
and zero constrained tip shortfall. Its largest moving-cavity mass residual is
2.106574280e-16 m³/s. It is byte-identical to the integration sample used to
compute the enhanced field.

The corrected full alphabet has 1,115,520 frames at 48 kHz (23.24 s), 2,375
cues, and no dictionary guesses. It has zero geometry-solve failures,
nonfinite samples, or clipped samples. Peak output is 0.690 without
normalization; maximum sampled body overlap is 1.00000008e-10 cm, maximum
unmet tip travel is 1.3733327 cm, and maximum cavity mass residual is
5.305204518e-16 m³/s. The launch binary and native source hashes were verified
against the final working tree. These numerical checks are not pronunciation
scores.

Local artifacts are under `out/` (not checked into Git):

- `canonical-conservative-quebec-current.wav`, `.csv`, `.log`.
- `canonical-conservative-nato-corrected.wav`, `.log`, and launch provenance.
- `canonical-conservative-quebec-enhanced.html`: full human and synthetic fields.
- `human-quebec-enhanced.{png,npz,json}` and
  `canonical-conservative-quebec-current-enhanced.{png,npz,json}`.

The human reference is the full `samples/derived/nato_q.wav` recording by
valeatory, with attribution and license in [samples/README.md](../samples/README.md).
It has 684×78 enhanced cells; the full synthetic sample has 684×131. Both retain
magma, 70 dB, and their saved amplitude references (human 0.37113784342106715;
synthetic 0.015974152821882408). References differ, so displayed brightness is
not a calibrated relative sound-level measurement. PNGs map every cell to one
row pixel and four column pixels, verified against the saved linear arrays.
There is no ordinary-STFT companion panel, grayscale conversion, or omitted
sample interval. The renderer now defaults to the full WAV; excerpts require
an explicit duration.

Viewer labels now come from the embedded recording metadata, preventing the
old hardcoded Alfa labels from misidentifying Quebec. The embedded PNG bytes
were checked against the full arrays. The in-app browser blocked the local
file protocol, so this change has no fresh browser interaction claim.

The earlier Papa–Quebec comparison reduced the largest modeled closing
cavity displacement from about 3166.5 to 612.1 mL/s; peak modeled cavity pressure
changed from 830.3 to 874.2 Pa. Those are consequences of different modeled
geometry, not evidence of calibrated human airflow or better pronunciation.

## Validation and remaining work

The final Release suite passed all 10 tests. New tests cover mapped attachment,
contact-limited reach, continuity at vanishing activation, dense-contour ray
comparisons, unsupported geometry, explicit render failure, chunk-independent
streaming, and NATO compound planning. Against an independently intersected
8,192-vertex contour at three activation levels and 351 airway rays, hit
classification agreed and the largest gap difference was 8.31e-6 cm.

The earlier full ASan/UBSan suite passed 9/9 before the contact constraint was
added. After that change, the affected tongue-surface and streaming tests
passed 2/2 under ASan/UBSan. The subsequent lexical change passed its separate
ASan/UBSan planner test. These are separate validation scopes, not a claimed
fresh full sanitizer run of the final tree.

Still unresolved: attainable combined body/blade motion, calibrated tissue
forces and timing, transverse anatomy and the physical gap-to-area law,
whole-articulator contact, and human pronunciation/realism. The next physical
work should explain why the requested tip placement is unreachable in this
geometry, rather than tune gain or timing around the missed constriction.

## Reproduce

```sh
cmake --build build
ctest --test-dir build --output-on-failure
./build/ducklips_say --canonical-motion --conservative-tongue \
  --body bodies/reference.body --dictionary samples/nato.dict \
  --text 'Quebec.' --out out/canonical-conservative-quebec-current.wav \
  --probe out/canonical-conservative-quebec-current.csv
python tools/enhanced_specgram/render.py out/canonical-conservative-quebec-current.wav \
  --out out/canonical-conservative-quebec-current-enhanced --channel 0 \
  --max-hz 8000 --display-reference .015974152821882408
```

Exact artifact hashes and final alphabet results are recorded in
[conservative_tongue_speech_reference.json](conservative_tongue_speech_reference.json).
