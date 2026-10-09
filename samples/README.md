# Human speech references

Only human spoken speech belongs in this reference collection. Generated
DuckLips output goes under `out/`, never into the reference set. No recordings
are used as synthesis units or copied into the generated waveform.

`manifest.json` pins the original 44.780-second NATO alphabet/digits recording
by **valeatory**, made 2014-02-28, to SHA-256. Source:
[Wikimedia Commons](https://commons.wikimedia.org/wiki/File:NATO_Phonetic_Alphabet_reading.ogg).
Selected license: [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/).
It is a contributor's human reading, not an official NATO voice or an
endorsement. The original is Vorbis, so high-frequency/low-level details can
include codec artifacts. It is one speaker, not a representative population.

`nato_segments.json` contains 26 source-time word boundaries supplied by
[emonster / SkyStamp](https://skystamp.app/license/), under the same license.
Its source SHA-256 agrees with the original recording. Credits and license
are in `licenses/`. Boundaries are externally supplied labels, not a new
listening validation or phone alignment by DuckLips.

Fetch and prepare with a Python having NumPy, SciPy, and SoundFile:

```sh
python tools/fetch_samples.py --prepare
```

The downloader has a 2 MB per-recording cap, verifies hashes, and refuses new
downloads below 1 GiB free space. `raw/` and `derived/` are ignored by Git.
Preparation selects microphone/channel 0, resamples to 48 kHz with a polyphase
anti-alias filter, and writes float WAV without level normalization. It also
extracts the 26 labelled words from the original decoded recording; those
modified audio files remain CC BY-SA 3.0. Receipts retain attribution, source
hashes, sampling details, and hashes of the resulting WAVs.

`nato.dict` is a working rhotic-English ARPAbet transcription guided by
[NATO's spelling alphabet](https://www.nato.int/en/news-and-events/articles/news/2017/12/21/nato-phonetic-alphabet-codes-and-signals).
It avoids crude letter-to-sound guesses. It is not a certified pronunciation
standard and should be checked by listening to each generated word. The
canonical written forms include **Alfa** and **Juliett**. Both `X-ray` and
`Xray` have explicit entries; known hyphenated compounds stay together in the
planner.

The pre-existing personal reference in `surface/out/ref/` stays local and
ignored. It has not been copied into this distributable reference corpus.
