# ducklips_claude

The speech work done with Claude (Claude Code) between 2026-10-08 and
2026-10-14. It was kept in `file_manager/research/ducklips/` on the laptop and
synced to `~/ducklips/` on the Mac mini, and was never committed until now.

Not included: every audio file (`*.wav`), the reference recordings and corpora
(`voice/ref/`), generated renders and analysis caches (`out/`, `*.npy`,
`*.npz`), build products, and spectrogram images of the owner's recordings.
[REFERENCES.md](REFERENCES.md) lists every reference the work used, where it
came from, its licence, and where it lived locally.

## Layout

The folders reference each other by relative paths (for example `tract`
builds `../voice/src/room.cpp`), so they are kept together here.

| Folder | What it is |
|---|---|
| `tract/` | **The main line.** The segment-module engine: `Speaker` / `VoiceBody`, a `Tube` with elastic, mass and support pieces per node, the triangular glottis (`Cord`), trachea, nose, sinuses, piriform and under-tongue branches, motor holds with wander. Tools: `speak`, `stream`, `voices`, `talk.sh`, `talk_remote.sh`. Also the sheet engine (`sheet_*`, the owner-rejected "vibrator" branch). `NOTES.md` is the full lab notebook; `GOAL.md` is the standing goal and state; `analysis/` holds the measurement scripts (see `analysis/README.md`). |
| `voice/` | The text front end and planner (`say`: text → phones → timed articulator events, "the body language"), the room stage, the phone tables and learned pitch and duration models (`data/`), and evaluation tools (`tools/`: phone identification, ARCTIC rendering, four-view enhanced spectrograms). `NOTES.md` is its notebook. |
| `voice/enhanced_specgram/` | The owner's high-resolution spectrogram method (MIT, see its `LICENSE`), used for the four-view images. |
| `surface/`, `boundary/`, `src/`, `include/`, `learned/`, `tests/`, `tools/` | Earlier engines and notes from 2026-10-08 to 10-10: the modal "surface" network, the boundary lab, and the original DuckLips front-end engine. `README.frontend.md`, `DESIGN.md` and `FINDINGS.md` describe it. `learned/` summarizes what was studied from VocalTractLab and from Astra's work. |

## Build (on the Mini)

```sh
cmake -S voice -B voice/build -DCMAKE_BUILD_TYPE=Release && cmake --build voice/build
cmake -S tract -B tract/build -DCMAKE_BUILD_TYPE=Release && cmake --build tract/build
tract/tools/talk.sh "Hello, my name is Papa." out.wav
```

## Where it stands (2026-10-14)

See `tract/GOAL.md` (state) and the end of `tract/NOTES.md` (owner feedback
on m43 / m44 and the fixes: gliding intonation, tremor band, fold posture,
nasal sinuses, skin low-pass, sibilant lips, coda velar place).
