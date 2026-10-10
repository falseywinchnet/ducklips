# Below-mouth feedback stage — 2026-10-10

The user's saved six-section oscillatory tube is preserved verbatim in
`web/tube/presets/user-vocal-tube-2026-10-10.json`. The default now extends it
into a nine-section below-mouth system. Mouth, sinuses, under-tongue branches,
teeth and lips are not implemented in this change.

## Configuration

The added upstream sequence reflects the downstream material order, with
half-length copies to introduce feedback without doubling the tube length.
Total length is 263.5 mm. The table gives nominal multipliers, before the
whole-body elasticity control and continual variation.

| Section | Length mm | Mass | Elasticity | Support |
|---|---:|---:|---:|---:|
| Inlet expansion | 50 | 1 | 1 | 1 |
| Lower elastic adapter | 17.5 | 0.2 | 3 | 0.2 |
| Lower firm support | 22.5 | 4.7 | 0.2 | 3.2 |
| Lower heavy tissue | 17.5 | 3.5 | 1.5 | 0.2 |
| Vocal cord | 16 | 1.85 | 1.5 | 0.2 |
| Upper heavy tissue | 35 | 3.5 | 1.5 | 0.2 |
| Upper firm support | 45 | 4.7 | 0.2 | 3.2 |
| Upper elastic adapter | 35 | 0.2 | 3 | 0.2 |
| Throat expansion | 25 | 1 | 1 | 1 |

Main diameter is 5.372888 mm, both expansions 10.745777 mm, and the cord's
resting opening 1.184192 mm. Thus the user's 2x expansion refers to diameter,
not volume or area. No inserted blocks are present. The insertion UI is gone;
legacy saved blocks remain readable.

The material sheet has 60 axial x 8 angular positions and 480 computed modes;
60 modes are axisymmetric. The circumferential symmetry of ring controls and
uniform inlet allows an exact single-sector runtime. The full eight-sector
path is retained and tested. Nominal dry mode range is 52.269–1062.426 Hz.

Independent seeded material states combine smoothed noise, sawtooth relaxation
and signed bubble envelopes. Elasticity stays within +/-10% of each section's
set value. The whole-body multiplier spans 1–4. Variation changes the physical
response cache, damping, and tissue impedance/delay; it is not audio noise or
a prescribed glottal waveform. See `web/tube/MODEL.md` for the approximation
and constitutive-model limits.

## Fresh visual experiments

All audio below was newly rendered with the sheet engine. Every directory
contains its exact input design, raw float WAV, source snapshot, trace and
receipt. Full enhanced fields use the same 48 kHz input, 0–8 kHz range, hop256,
original magma map, and every computed cell. The field is not calibrated
acoustic power; output levels below come from the waveform itself.

| Local output directory | Change | RMS | Peak | Visual observation |
|---|---|---:|---:|---|
| `out/below-mouth-reference` | Exact user reference | 0.0236784 | 0.0394957 | Strong, mostly stationary bands |
| `out/below-mouth-mirrored-01` | Nine sections at 570 Pa, variation enabled | 0.0009957 | 0.0036096 | Continuous bending bands, substantially quieter output |
| `out/below-mouth-pressure-850` | Raise inlet to 850 Pa | 0.0079483 | 0.0311207 | Continuous bands retained, stronger output |
| `out/below-mouth-cord-350um` | Tighten resting cord radius to 0.35 mm at 570 Pa | 0.0162809 | 0.0868104 | Broad noisy interval followed by a more banded state; rejected as default |
| `out/below-mouth-final` | Final preset; exact mixed-boundary damping update | 0.0082439 | 0.0333692 | Six-second field inspected separately below |

The final choice keeps the user's resting cord opening and uses 850 Pa. The
final audio contains six seconds after a one-second pressure startup; there
is no normalization, EQ, pitch generator or added audio waveform. This is
material/feedback development, not a demonstrated natural human voice.

The reference replay at `out/below-mouth-reference-replay` is byte-identical
to the freshly rendered reference above. Its depiction is therefore exactly
the same field; no independent changed sound is being substituted.

## Representation and validation

The added sections initially fell below 10x realtime (9.3x browser; roughly
8.6x Node). Air transport was the expensive part. When all air paths have
one-sample delay, their histories now share a time-major plane, removing
repeated per-link circular indexing while retaining every path, direction
and FIR. Longer/fractional paths keep the general representation. Full
available histories transfer on geometry changes. This is tested for
sample-exact equivalence through both directions of that transition.

Final six-second core rendering measured **10.995x realtime**. The browser
startup core benchmark measured **11.3x** on both localhost and the LAN URL
on this Mac mini. These are synthesis-core measurements, excluding geometry
compilation, UI drawing and offline enhanced rendering; the MacBook Neo's
actual speed has not been measured here. Full asymmetric operation has not
been optimized to that throughput.

Validation covers modal completeness and residuals, exact original fixture
regressions, independently recomputed perturbed eigensystems, seeded bounded
variation, pressure-only edits preserving histories, topology changes,
1–4 elasticity running without faults, release after removing pressure,
rotational equivalence, and sample-exact worklet/LAN block-size parity.

Browser verification covers nine handles, support arrow-key adjustment and
undo (1.18 -> 0.98 -> 1.18 mm), preset switching to six reference handles,
1–4 whole-body elasticity, nominal mode display, and LAN Worker startup.
LAN page and module routes return HTTP 200. The LAN server listens on
0.0.0.0:8765; tested at http://192.168.1.18:8765/ from this host. This does not
establish reachability from the other physical laptop if its network differs.

The browser screenshot is `workshop.png`. Enhanced fields are inspection
artifacts in the local output directories, not embedded replacement UI panels.

The end-to-end LAN take is saved at
`out/flat-tube-web/2026-10-10T19-13-38.156Z-21333b/`. Its complete enhanced field
was inspected: sustained, gradually bending bands with uneven texture, without
the broad noisy interval of the tighter-cord experiment. Its source snapshot
matches all nine current browser/engine modules, the raw audio hash matches,
and its capture records no design edits during the two-second take. The page
returned a working `Full enhanced field` link after processing. Browser test
playback was stopped and the main workshop left ready with monitor at 25%.

All **27 tests passed**; see `tests.txt`. `evidence.json` contains the exact
recording hashes, source hashes and measured core timing receipts.

The complete final six-second field was inspected at native cell resolution:
684 rows x 1126 time columns, every sample from 0 through 288000 included.
Bands persist across the duration, bend gradually and develop uneven texture;
there is no broad noisy opening like the tighter-cord candidate. Several strong
low bands remain nearly stationary. This is evidence of an evolving sustained
regime, not proof of organic human phonation. The enhanced pipeline took
1158.27 seconds offline; its registration work grows faster than linearly with
image duration. The image hash provenance matches the final WAV exactly.

The six-second recording and full field are also available over LAN at
http://192.168.1.18:8765/renders/below-mouth-20261010/index.html . These served
files are hard links to the exact final artifacts, avoiding duplicate storage.
