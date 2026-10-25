# Analysis scripts (2026-10-13)

These scripts were written during the voice work and copied here from the
Mini's `/tmp` so they survive. Most expect the Mini's paths: `/tmp/es`,
`~/ducklips/voice/ref`, and the builds under `~/ducklips`. Each file's
docstring or header gives its usage.

## Kinds

- **Validated, used to steer.** Four-view images (`voice/tools/tf_view.py`
  over `voice/build/tf_analyze`), zoomed wideband
  (`voice/tools/spectrogram2.py`), and:
  - `slice.py`: the loudest vowels' narrowband spectra;
  - `tilt.py`: vowel tilt and stop bursts;
  - `bursts2.py`: release-aligned burst spectra;
  - `phone_spectrum.py` and `phone_level.py`: per-phone spectra and levels
    against bdl;
  - `vtune.py` and `vopt.py`: vowel formants against bdl, and the posture
    search;
  - `pervowel.py`: vowel voicing and pitch.
- **Tried and found not to separate real speech from ours.** Rejected for
  steering: `hear.py`, `motion.py`, `ripple.py`, `hnr.py`, `flux.py`,
  `vibrator.py`, `f0mod.py`.
- **Score post-processes** (the prototypes for the planner's motion and
  intonation): `stretch.py`, `contour.py`.
- **Drivers:**
  - `four.sh`: the owner's four test sentences through the segment engine;
  - `seg_search.sh` and `seg_search_p.sh`: three-seed identification;
  - `seg_try.sh`: one seed with phone rows;
  - `r29.sh` and `full29.sh`: the sheet engine's batch metrics.

See `NOTES.md`, Milestone 4, for what each one found.
