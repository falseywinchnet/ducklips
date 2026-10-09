# NATO phonetic alphabet audio

## Attribution and license

**NATO Phonetic Alphabet reading.ogg**, recorded by **valeatory**, 28 February 2014.

- Original: https://commons.wikimedia.org/wiki/File:NATO_Phonetic_Alphabet_reading.ogg
- Creator: https://commons.wikimedia.org/wiki/User:Valeatory
- License: Creative Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0).
- License link: https://creativecommons.org/licenses/by-sa/3.0/
- Full legal terms: https://creativecommons.org/licenses/by-sa/3.0/legalcode

The source offers a choice of CC BY-SA 3.0 or GFDL. SkyStamp uses **CC BY-SA 3.0**.

Modifications by **emonster**: extracted the 26 letter words (Alfa through Zulu),
trimmed surrounding silence with short padding, and converted the stereo Ogg
Vorbis recording to mono AAC in individual M4A files. The spoken digits at the
end of the original are excluded. No speech synthesis or voice alteration is used.

The resulting `nato-a.m4a` through `nato-z.m4a` clips and emonster's audio edits
are distributed under **CC BY-SA 3.0**. You may copy, share, and adapt them,
including commercially, subject to that license's attribution and ShareAlike
conditions. Retain the title, creator, source, license link, and modification
notice when redistributing. Do not imply endorsement by the creator, Wikimedia,
NATO, or ICAO. No endorsement is claimed here.

SkyStamp's proprietary app license and its restrictions do **not** apply to these
audio assets. No additional legal restrictions or DRM are applied to the downloadable
audio pack. This license notice covers the audio; it does not relicense SkyStamp's
application code or unrelated assets.

## In-app access

Credits and clickable source/license links appear in:

- Quiz → NATO Phonetic Alphabet → Information.
- Settings → Company Info → License Agreement → Phonetic alphabet audio credits.

Both locations link to the website's audio attribution section. The app bundles
all 26 clips, this notice, and the full license for offline use; it contains no
audio download archive or export control.

## Website download

The audio attribution section at https://skystamp.app/license/#phonetic-audio
(and https://skystamp.app/ja/license/#phonetic-audio in Japanese) offers
`SkyStamp-NATO-audio-CC-BY-SA-3.0.zip` without an account. The pack contains all
26 clips, this notice, the full license, and clip timings. The section is also
included on the other localized license pages with explicitly marked English copy.

The generator writes the archive to `website/artwork/phonetic-alphabet/`.
The website asset pipeline copies it to `public/downloads/`, preserving it across
site rebuilds. Download and reuse are website features, separate from quiz playback.

## Quiz behavior

Ten unique letters are sampled from A–Z. Sessions alternate five letter-to-word
questions with five audio-to-letter questions, each with four distinct choices.
The answer's pronunciation becomes available after answering a letter-to-word
question; listening questions offer playback from the start. Playback is initiated
by the user, can be replayed, and stops on leaving the question, backgrounding,
audio interruption, or headphone disconnection. Questions, scoring, results, and
audio work offline. Quiz interface copy is localized in all 23 supported app languages. Code words
are invariant across app languages, including
the aviation spellings **Alfa**, **Juliett**, and **X-ray**.

## Source and reproducibility

- Original recording: `Scripts/Audio/source/NATO_Phonetic_Alphabet_reading.ogg`.
- Original checksum, letter mapping, and padded cut boundaries:
  `Scripts/Audio/phonetic-alphabet.json`.
- Generator: `Scripts/Audio/prepare_phonetic_audio.py` (Python 3 and FFmpeg).
- Runtime assets: `SkyStamp/SkyStamp/Resources/PhoneticAlphabet/`.
- Full license in the repository: `Docs/licenses/CC-BY-SA-3.0.txt`.

Run from the repository root:

```sh
python3 Scripts/Audio/prepare_phonetic_audio.py --ffmpeg /path/to/ffmpeg
```

The generator rejects a changed source checksum rather than applying stale cuts.
If the recording or timings change, verify every word and boundary again before
shipping. Publish the matching website archive when updating the distributed clips.
