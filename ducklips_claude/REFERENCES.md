# References used (not included in this repository)

None of these files is committed here. Paths are where they lived on the
laptop (`file_manager/research/ducklips/`) and the Mac mini (`~/ducklips/`).

## Recordings and corpora

| Reference | Local path | Source and licence | Used for |
|---|---|---|---|
| **CMU ARCTIC, speaker bdl** (US English male): `wav/`, phone labels `lab/` (end times), prompts `arctic.data` | `voice/ref/arctic_bdl/` | Carnegie Mellon University, festvox.org/cmu_arctic. Free for any use under the light conditions of its `COPYING` (kept beside the data). | The phone-identification test (`voice/tools/phone_test.py`: kNN over bdl a0001–a0040, ceiling on a0041–a0060); vowel formant tuning (`vtune.py`, `vopt.py`); fricative, burst, nasal and tilt comparisons; cycle jitter (`tract/analysis/wmatch.py`); formant transitions (`pinch.py`); `render_arctic.py` renders the same prompts through our engines. |
| **CMU ARCTIC, speaker rms** | Mini only: `~/ducklips/voice/ref/arctic_rms/` (analysed into `voice/out/arctic_rms/`) | As above. | Fitting the learned vowel pitch and phone duration models (`voice/data/pitch_en.txt`, `durations_en.txt`; `voice/tools/pitch_fit.py`, `duration_fit.py`) and `identify_tune.py`. |
| **NATO phonetic alphabet reading** by "valeatory" | `voice/ref/nato/` (48 kHz, per-letter clips and the whole reading, with `nato_clips.json` provenance) and `voice/ref/nato16/` (16 kHz) | Wikimedia Commons, `File:NATO_Phonetic_Alphabet_reading.ogg`, **CC BY-SA 3.0**. First channel, resampled to 48 kHz; source sha256 15197e85…6f88. | NATO word comparisons (`nato_compare.py`), burst and aspiration levels, the room stage's level. |
| **The owner's own recordings** | `voice/ref/user_three.wav` ("hello / my name is papa / and I speak with a body …"), `voice/ref/user_nato.wav` (the NATO alphabet), `surface/out/ref/` and earlier `voice/listen/` comparisons | Recorded by the owner for this work. **Private; not to be published.** | Pitch contours and pace, H1−H2 and tilt, the room floor, cycle jitter (RAP 0.23–0.42%), NATO comparisons, four-view visual comparisons. |
| **VocalTractLab 2.3 examples** (`vtl2.3-example01…05-orig.wav` and `-synth.wav`) | Mini: `/tmp/es/vtl/` | Downloaded with the owner's permission from vocaltractlab.de (the examples page). The recordings belong to their authors. | The owner's reference for natural synthesis: visual comparison, cycle jitter and contour shape. |

## Data and models built into the code (cited, not copied audio)

| Reference | Where in the code |
|---|---|
| B. H. Story (2008), vowel area functions from MRI | `tract/data/story_2008_areas.txt` / `.json`; `tools/sustain` |
| Dang & Honda (1994), nasal tract areas (subject 4), and the sinus volumes and ostia | `tract/src/voice_body.cpp` (`kNoseAreas`, `kSinus*`) |
| Dang & Honda (1997), piriform fossa depth and zero | `tract/src/voice_body.cpp` (piriform branch) |
| VocalTractLab 2.4 (P. Birkholz): the triangular glottis's static fold parameters, the M01 speaker's pitch calibration, the glottal noise source, viscous duct loss | `tract/src/cord.cpp`, `tract/src/tube.cpp`, `tract/src/voice_body.cpp`; notes in `learned/vtl_*.md`. VocalTractLab is GPL: only its published mathematics and parameter values were studied; none of its code was run or copied. |
| Xu's target approximation model of pitch | `voice/src/planner.cpp` (pitch targets) |
| Astra's sheet-tube work (`~/ducklips_base`, read-only; repository commits up to 15797b5) | Ported in `tract/src/sheet_*`, `tract/src/variation.cpp` (from `variation.mjs` at 42c7337); see `learned/astra_lessons.md` |
| CMU Pronouncing Dictionary conventions (ARPAbet phones) | `voice/src/lexicon.cpp`, phone tables |
