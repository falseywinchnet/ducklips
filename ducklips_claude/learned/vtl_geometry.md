# VocalTractLab 2.4 — 3D vocal tract geometry and area function (study notes)

Source: VTL 2.4 `Sources/` (GPL; mathematics only): VocalTract.cpp/.h,
Surface.cpp, Geometry.cpp, Splines.cpp, AnatomyParams.cpp. Units cm, cm².
x anterior, y up, z lateral; only z ≤ 0 is built and mirrored. Anatomy values
are the built-in defaults (VocalTract.cpp:102–174); speaker files override.
Pipeline `calculateAll()` (2371): calcSurfaces → calcCenterLine →
calcCrossSections → crossSectionsToTubeSections → getTube() (7553).

## Parameters (VocalTract.h:160; ranges VocalTract.cpp:155–173)

| Param | Meaning | min..max (neutral) |
|---|---|---|
| HX | narrow↔wide larynx blend | 0..1 (1.0) |
| HY | hyoid/larynx height, cm | −6..−3.5 (−4.75) |
| JX | horizontal jaw shift, cm | −0.5..0 (0) |
| JA | jaw angle, deg | −7..0 (−2) |
| LP | lip protrusion (along lip-corner path) | −1..1 (−0.07) |
| LD | vertical lip distance, cm (<0 compression) | −2..4 (0.95) |
| VS | velum shape high↔mid | 0..1 (0) |
| VO | velic opening, cm² | −0.1..1.5 (−0.1) |
| TCX, TCY | tongue-body ellipse centre | −3..4, −3..1 (−0.4, −1.46) |
| TTX, TTY | tongue-tip circle centre | 1.5..5.5, −3..2.5 (3.5, −1.0) |
| TBX, TBY | blade Bézier control point | −3..4, −3..5 (2.0, 0.5) |
| TRX, TRY | root Bézier control point | −4..2, −6..0 (0, 0) |
| TS1, TS2, TS3 | tongue-side elevation root/back/tip | 0..1, 0..1, −1..1 (0) |

## Articulator surfaces

All articulators are rib grids (numRibs × numRibPoints, quads split in two
triangles; Surface.cpp:197–208), ribs mostly rational quadratic Béziers
resampled at equal arc length (Splines.cpp:200–270, 1025–1046).
- Larynx: narrow/wide 8-point outlines blended by HX, sheared along the
  oblique pharynx wall (fulcrum (−2.372, 1.214), a = −98°), shifted by HY.
- Pharynx back wall: 3 ribs, half-depth 3.4/2 → 3.8/2.
- Velum: LOW/MID/HIGH 6-rib surfaces blended by VS and VO; acoustic port area
  is `max(VO, 0)` directly.
- Hard palate: 9 ribs from per-rib palate height and angle; blended into velum
  for x < 1.2 cm.
- Mandible/lower teeth: same construction mirrored; rigid rotation by JA about
  (−6.5, 2.0) after rest offset (0, −1.2) + JX. Teeth: 25 ribs with 0.15 cm
  interdental grooves.
- Lips: corner path lerp(narrow, wide, LD) rotated by **half** of JA (lips
  follow the jaw at 50%); `h = max(0.5·LD, −0.05)` is half the lip gap; JX < 0
  raises both lips (labiodentals); lip edges clamped between incisor tips.
- Cheeks: fill strips at the teeth line — **no buccal space** in the model.
- Tongue (`calcTongue` 3806–4249): midsagittal contour of root Bézier (hyoid →
  (TRX, TRY) → tangent on body ellipse, r = 1.8) + ellipse arc + blade Bézier
  (→ (TBX, TBY) → tangent on tip circle r = 0.2) + tip arc; 33 dynamic + 4
  static ribs × 11 points. Root automatic: `TRX = 0.938|C − H| − 5.11`,
  `TRY = 0.831·TCX − 3.03`. Tongue sides: elevation `min(p/0.3, 1)` cm for
  p ≥ 0, `max(−0.5·p/−0.15, −0.5)` cm for p < 0, knots along the contour,
  coronal rib a cubic smoothstep over half-width 1.75 cm; tonsil bump.

## Centre line and cross-sections

- Rough centre line: vertical behind the tongue, quarter ellipse around it,
  horizontal to the lips; 129 points; ceiling/floor intersections; midpoints
  smoothed over 2 cm arc length; normals de-crossed by bisection.
- Each cut plane intersects all triangles; profiles sampled at 96 × 0.0729 cm
  laterally; upper profile from upper cover/teeth/lip/uvula, lower from lower
  cover/teeth/lip/epiglottis; teeth rule removes the vestibule beside the
  incisors; lip closure collapses the profile; mirrored.
- Tongue merged as `lower = max(lower, tongue)`; lateral ends trimmed while
  lower > upper − 0.1 cm → a tongue within 1 mm of the hull everywhere leaves
  area 0: **that is how stops close**. Lateral channels: one lumen per side kept
  (gap ≥ 0.01 cm after skipping < 0.1 cm), so /l/ keeps side channels with the
  midline closed.
- `A = Σ ½(gap_i + gap_{i+1})Δ`, perimeter from both profiles (computed but
  the acoustics use `2√(πA)`).
- **Minimum-area floors**: incisor zone [incisor − 0.5, + 0.3] cm ≥ 0.15 cm²;
  tongue-side floor `0.15(p − 0.2)/0.2` up to 0.15 (p > 0.2, fricative
  groove), `0.20(p + 0.2)/(−0.2)` up to 0.20 (p < −0.2, lateral); |p| ≤ 0.2 →
  0, full closure only there. TS2 applies to the back, TS3 from tip − 2 cm to
  the lips. The /s/ groove is a parameter, not geometry.
- Tube sections: 18 pharynx (glottis → nasal port) + 32 mouth; lengths vary
  with shape; each section takes the **minimum** area over its span and that
  section's articulator tag; `MIN_AREA_CM2` 1e-4.

## Other cavities
- No separate sublingual or buccal branch.
- Piriform fossa: 5-section tapered side tube, `maxArea = 2V/L`, L 2.5 cm,
  V 1.5 cm³ (defaults), coupled at pharynx section 4.
- Nose: 19 × 0.57 cm (scaled), Dang & Honda subject 4 areas
  {1.30,1.25,1.30,1.20,1.60,2.50,3.60,2.25,2.50,1.75,2.50,2.10,1.50,1.50,0.75,
  1.25,1.60,1.65,1.40}; first 4 sections ramp `A_port + i²(A_4 − A_port)/16`.
- Sinuses: 4 Helmholtz resonators, V 11.3/6.8/33.0/6.2 cm³, neck length
  0.3/0.3/0.45/1.0 cm, neck area 0.185/0.185/0.145/0.11 cm², at nose sections
  8, 9, 11, 12.
- Trachea 23 × 1 cm, 2.5 cm² (first two 4.0, 3.0); glottis 2 × 0.3 cm.

## Shapes and constraints
- Speaker file `<shapes>`; consonants `tt-`, `tb-`, `ll-` stored as C(a),
  C(i), C(u); interpolation as in vtl_gestures.md.
- Tongue limited against the hull (pushed 0.3 cm out; ellipse moved toward an
  anchor until contact), rib points clamped to [lower + 0.02, upper − 0.02];
  order and box constraints for tip, blade, root; epiglottis rotates to clear
  the tongue root; HY below the rear jaw rib.

## CANDIDATE ideas for our engine
1. Context prototypes for consonants (a/i/u) + downstream dominance.
2. Minimum-area floors tied to tongue-side bracing (fricative 0.15, lateral
   0.20, incisor 0.15 cm²); closure only when bracing is off.
3. Section area = minimum over its span; grid anchored at anatomical
   landmarks (nasal port).
4. Contact as a geometric clamp with a 1 mm closure threshold.
5. Articulator tags per section to place and shape noise sources.
6. Fixed side branches (piriform, 4 sinuses, 19-section nose with quadratic
   port ramp).
7. Ellipse + tip circle tongue with tangent Béziers and hull-aware limits.
8. Lips at half the jaw rotation; JX < 0 raises both lips.
