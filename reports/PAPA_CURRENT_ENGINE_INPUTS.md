# Fresh Papa inputs for the coupled engine

The next word input is freshly generated from `Papa.`, `samples/nato.dict`
and `bodies/reference.body`. This prepares an actual utterance for the current
coupled source/mouth path. It does not render the previous speech engine or
produce a new listening sample.

`tongue_reach_audit --source-cues FILE.json` now exports the planner's timed
physiological targets and phone intervals alongside its geometry exports.
The targets include lung pressure, glottal rest opening, intended pitch,
effort, loudness and velum opening. They are behavioral cues, not an output
waveform or an imposed fold cycle. Exporting them leaves the geometry,
airway and lip CSV bytes unchanged. The acoustic source still needs to
consume these controls with correct actuation and work accounting.

The dictionary gives `P AX0 P AA1`. The planned phones occupy 0.600–1.01572 s;
the complete input also includes preparation and a tail. There are 1,366
geometry frames through 1.365 s. Native region decomposition and independent
geometry/cap checks pass every exported frame. Each is compatible with the
current adapter's two oral regions, four lip regions and eight throat regions.
The smallest checked storage volume is 9.0819e-9 m3. No disappearing lip or
throat storage or mismatched oral-region ownership was found at these frames.

The lip outlet closes at frame 0.532 s and reopens at 0.670 s, then closes
again at 0.742 s and reopens at 0.795 s. These are sampled geometry results;
they do not establish topology invariance between all frames or acoustic
correctness through the releases.

The 77 exported source cues expose unfinished integration explicitly. The
breath prelude opens the velum to 8e-5 m2 at 0.100 s and commands closure at
0.500 and 0.550 s. The present coupled adapter has no nasal branch. Ignoring
that input would not reproduce the declared body/gesture sequence. Glottal
rest-opening and intended-pitch targets also need physical source actuation;
only time-varying lung pressure is currently exposed by that adapter.

A held native load at the word's 0.900 s vowel supplies a separate onset
diagnostic. With tension 1, rest opening 3e-6 m2 and 600 Pa, the independent
local search finds a decaying mode near 248.76 Hz (−197.21/s). At 1750 Pa it
finds a growing mode near 347.66 Hz (+34.63/s). A tension sensitivity check
at 600 Pa also finds decaying modes at the sampled tensions. These are
local searches from specified seeds, not exhaustive stability proofs or
simulations of muscle-driven onset. Neither pressure nor body defaults have
been changed to turn these diagnostics into a voice claim.

The [input receipt](papa_current_engine_inputs.json) binds the fresh exports,
full geometry preflight, source cues and load diagnostics. Acoustic rendering,
human intelligibility and realism remain unverified.

```sh
build/ducklips_tongue_reach_audit bodies/reference.body samples/nato.dict 'Papa.' \
  --airway out/papa-current-articulation/airway.csv \
  --lip-profile out/papa-current-articulation/lips.csv \
  --source-cues out/papa-current-articulation/source-cues.json \
  > out/papa-current-articulation/geometry.csv
build/ducklips_airspace_geometry_probe out/papa-current-articulation/geometry.csv \
  0 2 256 out/papa-current-articulation/volumes.csv \
  --frames out/papa-current-articulation/poses --core bodies/reference.body
OPENBLAS_NUM_THREADS=1 /Users/joshuahkuttenkuler/miniforge3/bin/python \
  tools/word_geometry_preflight.py out/papa-current-articulation --word Papa
```
