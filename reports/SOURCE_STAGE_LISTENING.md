# Source isolation and cumulative output listening

The user requested the vocal folds and their interactions in stages to locate
the unwanted timbre. The earlier 019 tap retained the Papa controller and
loaded airway. These new recordings remove that confound explicitly.

All recordings are fresh, 1.3 s at 48 kHz, with the same one-strip body and
round 018 source-pressure equations. There is no text/word planner. Lung
pressure is 600 Pa, rest opening 3e-6 m² and tension Q=1 throughout. Jitter,
shimmer, tremor, aspiration noise and frication noise are disabled. The source
is self-oscillating; no pulse waveform is supplied. The start of each file
includes its physical onset from rest. The clip endpoint is a recording
boundary, not a simulated vocal release.

The free source uses the existing local inertive inlet into a constant-pressure
reservoir: compliance 1e100 m³/Pa and all airway branch conductances zero.
This is a zero-load diagnostic boundary, not a physical excised-larynx model.
The finite-flow, two-mass pressure equation, two tissue substeps and one local
air step are the same as in the loaded cases.

The coupled runs use the actual Synthesizer and NudgeStream, with two isolated
diagnostic additions: constant source controls bypass the word/pitch modulation,
and the initial tract shape can be AX or AA. The glottis retains airway pressure
feedback. The closures condition uses AA with lip-only activation/release cues
at .35/.50/.80/.95 s. There is no glottal abduction, pitch movement, breath
modulation or consonant noise added with those closures. It should not be
represented as a correctly coordinated /p/ or a new Papa success.

All source/interior/lip observation files use rho/(4 pi * .4 * 128) times the
local volume-flow derivative. Interior taps are at oral branches 5, 11 and 17
(after sections 6, 12 and 18 of 24). These are diagnostic audifications of flow
change, not direct internal microphone pressures or separable causal stages.
The whole airway remains coupled for every interior observation.

The cumulative external composites apply the existing AirPath to (1) lips,
(2) lips plus nostrils, (3) lips plus nostrils plus the existing tissue-weighted
wall radiation. The third reconstructs the engine's microphone audio exactly
in all three conditions, with zero float sample difference. Adding the nasal
contribution is byte-identical at PCM16 for the closed velum; this is not a
missing clip. Interior signals are never added together.

A single playback gain 315.98259903736243 gives a maximum peak of 0.7 across
the set. No per-clip normalization, EQ, fades or cuts. The two ordered listening
sequences concatenate whole clips with 350 ms of silence. Each complete unique
constituent receives the preserved enhanced analysis (128-sample hop, all 684
frequency rows through 8 kHz, original magma). Exact PCM duplicates share the
same field. These are internal depictions, not user-facing image comparisons.
The packaging process records completion; do not claim inspection before the
fields have completed and been viewed.

Audio, raw recordings, physical traces and the manifest are under
`out/source-stage-listening`. The instrumented source tree and executable are
under `out/source-stage-study`. It starts from 018's exact source snapshot and
applies `source_stage_instrumentation.patch`; `tools/source_stage_listening.cpp`
uses those diagnostic headers. It is intentionally not a main CMake target.
Its newly compiled synthesizer/nudge objects link to the frozen 018 library;
`libducklips-base.a` preserves that dependency. All live-Synthesizer nonfinite,
clipping, geometry and airflow audits pass. Detailed provenance is in
`source_stage_listening.json`.
