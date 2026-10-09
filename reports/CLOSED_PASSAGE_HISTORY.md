# Acoustic history at complete contact

The user's alphabet assessment identified clicks in Foxtrot and November.
Tracing the exact heard rendering exposed a separate acoustic boundary defect
from the [single mean-cavity pressure problem](CAVITY_TOPOLOGY.md).

An acoustic passage with zero area has zero inverse inertance G and carries
no volume flow. The canonical integrator nevertheless continued integrating
its momentum from the pressure difference across it. That momentum was
unobservable while the passage was closed (`Q = G momentum = 0`). When the
passage reopened, it became a large flow. Its viscous/radiation relaxation
history also remained attached to an air path that no longer existed.

The canonical speech path now removes momentum and relaxation state for a
fully closed passage, after each substep and at the geometry endpoint. It
preserves nodal pressure and the history of every open passage. This changes
the actual speech engine; there is no waveform edit, gain change, word-specific
condition or adjustment to the planner's gestures.

## A regression that fails before the fix

Start at rest with a sealed section in a short tube. Apply a 500 Pa series
pressure source across the absent passage for 50 ms. No air moves and the
reported stored energy remains zero. Turn the source off, then reopen over
1 ms. Both air spaces remain unpressurized, so there is no physical source
state to release.

| Measured after reopening | Before | Corrected |
|---|---:|---:|
| Maximum stored energy | 0.116659618 J | 0 J |
| Maximum lip flow acceleration | 149.955962 m³/s² | 0 m³/s² |

The same pressure source acting through an open passage still stores
1.656497682e−5 J, identical before and after. Its state persists when the
source stops. The regression therefore checks an absent-passage boundary,
not blanket muting of pressure sources or acoustic state.

The test was compiled against the old native library and failed as shown,
then passed against the corrected library. It is now
`ducklips_canonical_closure_tests` in CTest. The full Release suite passed
12/12, including host-buffer independence and live nudge streaming.

## Connection to the heard audio

The diagnostic rerender `out/nato-click-audit.wav` is byte-identical to
`out/nato-local-material.wav`, which the user assessed. Its SHA-256 is
`9836e9a7b5b135b2dd030a68d24f483c0aa9d9d4078694bdeb6b18da7aa389b8`.
Both use the canonical network and conservative material tongue.

In Foxtrot the largest left-channel sample jump is 0.767327264 at 5.0660625 s.
The corresponding lip-radiation candidate is −6.153250 m³/s² at 5.065146 s,
around the opening /t/. In November, the strongest candidate occurs as the
initial /n/ releases into /ow/: lip acceleration −0.9479232 m³/s² at
11.846167 s, followed by the largest microphone sample jump at 11.847 s.
These timestamps identify waveform candidates; the user did not give exact
click timestamps. Propagation from lips to microphone accounts for a delay.

The new alphabet uses the same body, dictionary, cues, seed and unnormalized
output. Its 23,239 geometry-probe steps match the prior rendering exactly
in constriction area, attachment coordinates/motion, tip miss and body overlap.
The full run has zero geometry failures, nonfinite samples and clipping.
ASan/UBSan passed the affected streaming and new closure tests, 2/2.

| Whole word, same interval | Before | Corrected |
|---|---:|---:|
| Foxtrot peak | 0.529955 | 0.140136 |
| Foxtrot largest left sample jump | 0.767327 | 0.237927 |
| Foxtrot RMS | 0.015874 | 0.013427 |
| November peak | 0.141429 | 0.051293 |
| November largest left sample jump | 0.029986 | 0.021027 |
| November RMS | 0.011768 | 0.011459 |

At the original Foxtrot /t/ event (5.060–5.075 s), the largest sample jump
falls from 0.767327 to 0.131492. The largest remaining Foxtrot jump moves to
4.9620417 s, near /k/ → /s/. At the November /n/ release (11.840–11.860 s),
the jump falls from 0.029986 to 0.007163. These reductions support a causal
contribution from the closed-passage state, but do not establish that all
audible clicks have gone. The remaining Foxtrot event needs further work.

The excerpts preserve the exact full-sequence samples: Foxtrot 4.55–5.65 s
and November 11.70–12.65 s, including adjacent quiet time. They have no fades,
normalization or waveform repair. Files are
`out/foxtrot-closed-contact-context.wav`,
`out/november-closed-contact-context.wav`, and the complete
`out/nato-closed-contact.wav`.

Full enhanced fields for each before/after word use the same saved reference
0.015974152821882408, magma and 70 dB. Foxtrot has 684×104 cells; November
684×90. PNGs contain every cell, one pixel per row and four per column; exact
pixel comparison against the saved NPZ fields passed for all four. No
conventional STFT companion is displayed. Their filenames are
`out/{foxtrot,november}-{before-contact,closed-contact}-enhanced.{png,npz,json}`.

Source/binary/audio hashes, exact render command, intervals and measurements
are in [closed_passage_reference.json](closed_passage_reference.json).
The old listening assessment is not a score for the new recording.

## Remaining physical limits

This is an exact-contact boundary rule. It is not a validated continuum
model of fluid being squeezed out during closing contact, nor a smooth
near-contact viscothermal limit. An exploratory tiny-positive-area sweep
still shows large energy changes when the loss-memory coefficients change
with geometry. Those states need a physically justified moving-geometry
interpretation; resetting only fully closed passages does not solve that
problem.

The production renderer also still uses one mean-pressure cavity. Separating
its pressure storage, validating the gap-to-area law and improving source
and articulation realism remain necessary. A reduction in waveform spikes
does not establish correct pronunciation or human realism.
