# Local engine recovery attempt — listening failed

The coupled gas-source solver had become the speech implementation. Its latest
Papa attempt stopped after 135 ms of completed preparation, with no complete
word. Its data is preserved at `out/papa-coupled-oral-render`; it is no longer
the path used to produce listening samples.

`LocalSpeech` is a new audio core selected by `--local-flow`. Both the text CLI
and the incremental nudge stream use it through `Synthesizer`. It retains the
existing body configuration, articulator trackers, conservative tongue,
pressure-driven two-mass fold mechanism, pitch control and microphone paths.
It replaces the modal/canonical acoustic evolution and separate mean-pressure
system with local air storage and branch momentum. It does **not** run the
recent continuous-glottis adaptive solver, and it is not acoustically equivalent
to that reference. No recorded waveform or scripted fold cycle drives it.

## Representation and update

Each air volume stores excess air charge; each passage stores momentum.
Neighboring volumes exchange air through a fixed local midpoint update. Its
coefficients are scalar closed forms. A forward sweep followed by a reverse
sweep gives a symmetric composition with a fixed number of operations. There
is no global acoustic matrix factorization, eigenbasis, Newton loop, rejected
timestep or convergence controller in this audio update.

The rejected recording used an instantaneous Bernoulli inlet feeding a tiny
compliance. The subsequent [failure investigation](LOCAL_SPEECH_FAILURE.md)
adds the missing inertia of air crossing the folds; this has not established
intelligible speech. The same air state supplies pressure back to
the moving folds. Nasal and side branches participate directly. Mouth closure
therefore stores pressure in the same system that subsequently produces its
release. Nodal air is retained at a seal; absent passage momentum is cleared.

The current engine uses circular Poiseuille drag, separation loss at the
selected constriction and external openings, and the passive exterior
resistance/inertance load. Separation loss is included in the local pressure
exchange, not applied after an independently driven flow. Smooth sections do
not each become a separate turbulent jet. Yielding walls now have a restoring
compliance as well as mass and resistance: steady pressure deforms them rather
than continuously losing air through them.

This is a reduced physical model. The continuous source's internal compressible
gas states and the reference's frequency-dependent thermal/viscous memories are
not ported. The airway still uses the empirical upper-gap section model and its
compliance floor. Moving-geometry work, disappearing-volume transport and
turbulence energy are not given the reference solver's complete work ledger.
The static conservation tests below do not establish those missing properties.

## Body representation

The local path queries the carried 128-point tongue contour directly rather
than refining every ray against an implicit mapped circle. The smooth surface
remains available in the older/reference paths. For the tested undeformed
circle, ray gaps differ by less than 20 micrometres; this is not an error bound
for every deformed surface.

Two parts of the material map now use their exact motion: stationary tissue
below the anchor or outside the attachment's entire swept support, and vertical
shear wholly above the deformation gate. These retain the exact deformation
gradient and unit area determinant. The remaining body construction still
uses the previous material-map iterations. It is not yet a fully persistent
mechanical body engine, and remains the main time cost on complex consonants.

## Validation and fresh recordings

`ducklips_local_speech_tests` checks analytic oscillator refinement, energy and
air conservation, passive drag, exact sealed flow, elastic-wall static balance,
polygon/smooth-circle agreement, bit-identical chunked nudge rendering, and no
render-time heap allocation. It also checks that finite speech did not depend
on invalid-sample substitution or clipping. Release and address/undefined
sanitizer runs pass. Existing stream and tongue surface tests pass; the material
map change also passes the tongue-flow and tongue-surface suites.

A separate uniform-tube check compares measured modes with the analytic
discrete tube. With 24 oral sections and one fixed local substep, the first four
frequency errors are approximately 0.043%, 0.145%, 0.327% and 0.596%; halving the
step improves them. This supports the selected speech-band representation in
that simple geometry. It does not validate a moving constriction or establish
perceptual equivalence to the 64-section engine.

The latest frozen run is `out/recovered-engine-20261010-final/run.json`.
It includes the executable, body, dictionary, complete source snapshot, exact
commands and hashes for freshly rendered Papa, November, Foxtrot and the NATO
alphabet. Every native render completed with zero nonfinite samples and zero
clipping at the declared recording scale. No old speech audio was reused.

The recordings use 128 Pa per native full-scale sample to preserve the large
consonant transients without clipping. Each listening copy uses one fixed gain
for the entire waveform, recorded in the manifest. There is no time stretch,
looping, denoising or dynamic gain correction. This headroom does not repair
the excessive transients. The original 2 Pa scale clipped on several words.
Human listening rejected all four recordings: Papa had high-pitched noise and
a cut-off opening followed by “pa”; November was “click-click-ember”; Foxtrot
was only noise; the alphabet was clicks and garbage. The first two associations
follow the presentation order. No pronunciation is approved for this engine.
The passing numerical checks above did not establish intelligibility. The exact
feedback and recording hashes are retained in `engine_recovery.json`.

See `engine_recovery.json` for measured render times. Preparation is separate.
Papa exceeds 10× realtime in the measured run; the full alphabet does not.
Measurements of body-construction time identify the remaining bottleneck;
they are not an exemption from the whole-engine performance requirement.

```sh
cmake --build build --target ducklips_say ducklips_stream ducklips_local_speech_tests
./build/ducklips_say --local-flow --conservative-tongue \
  --body bodies/reference.body --dictionary samples/nato.dict \
  --text 'Papa.' --out out/papa-new.wav

# Historical reproduction tool; the alphabet below failed human listening.
OPENBLAS_NUM_THREADS=1 /Users/joshuahkuttenkuler/miniforge3/bin/python \
  tools/render_local_engine.py --out out/a-new-recovery-run --alphabet
```

`--local-sections` and `--local-substeps` provide explicit refinement controls.
The text CLI's local defaults are 24 oral sections, 9 nasal sections and one
fixed substep; the streaming CLI uses the same configuration. The unflagged
CLI and the default library options remain the older modal path for now.
