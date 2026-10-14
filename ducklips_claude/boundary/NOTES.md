# Boundary-surface sibilant prototype: notes

Status: research probe, 2026-10-08. It is not a candidate engine and has no
DuckLips integration.

Environment: Apple M4 Mac mini, one thread, Apple clang -O2, no fast-math.
The internal rate is 192 kHz (4x), patches are 4 mm (2.5 mm in the duct), and
c = 350 m/s. Reproduce on the Mini with `./build/boundary_lab validate`,
`run DIR [--face-cm 2 --full-kernel 1]` and `stability`. Outputs land in
`build/out/`, because the Mini's `out/` is erased by every `m4build` sync.

## Formulation

The model treats walls as planar quads split into patches (strips along u,
joined at their edges along v). The pressure on each patch comes from a
time-domain Kirchhoff collocation over every patch it sees:
2R sum A [cos b (p/r^2 + p'/(c r)) - rho a/r]/(4 pi) at the retarded time.
Delays are linearly interpolated, and p' is a central difference. Soft
tissue is a strip membrane with its own modes (Jacobi), driven by p, and
its acceleration is re-emitted. The microphone is another receiver, with
factor 1. The jet is +/- noise pressure on the two faces of the top row of
the lower incisors, which acts as a Curle dipole on the tooth. Control: a
1-D FDTD tube built from the same numbers with the same jet noise and a
baffled R-L radiation load.

## Results

| Test | Result |
|---|---|
| Closed rigid box 4x1x1 cm, R=1, 2.5 mm | modes 4348/8719/13102 Hz vs 4375/8750/13125 (**MEASURED**, <1%) |
| Same box: losses | B 480-980 Hz where the continuum has none. Finer patches: B grows (609-979 Hz). **MEASURED** numerical dissipation |
| Same box, 5 mm patches or 8x oversampling | NaN: late-time MOT instability. It appears when c dt / patch is below ~0.7 (**MEASURED**) |
| Open duct 2 cm, culled kernel | 2859 -> 2824 Hz as patches shrink from 2.5 to 1.25 mm; reference 3529-3733, 1-D 3557. Converges to a 20% error, B ~1700 Hz (**MEASURED**) |
| Open duct + 2 cm baffle, culled | identical (2859 Hz): under the visibility rule the interior can never hear the baffle |
| Open duct + baffle, full signed kernel | 3627 (2 cm), 3445 (4 cm), 3662 Hz (1.67 mm). B 609-727 vs 275 for 1-D (**MEASURED**) |
| /s/, /sh/ culled, stable | /s/: peak 5.8-6.2 kHz, B 5.6 kHz. /sh/: peak 4.8-5.2 kHz, B 3.7 kHz. /sh/-/s/ = -1.5 dB. The spectrum is featureless (`out/boundary_s_sh.png`) |
| /s/, /sh/ full kernel (with/without face, 2-4 mm teeth, 2.5 mm patches) | NaN in 5 of 6 variants; the sixth grows by +880 dB |
| full kernel + two-point averaging | stable, but B 5-6 kHz: no resonance left. Duct B 3082 Hz (4x) / 1547 Hz (8x) |
| 1-D tube control, same geometry | /s/ 3773 + 10236 Hz, /sh/ 1969 + 7295 + 9416 Hz, /sh/-/s/ +0.3 dB. B 230-700 Hz |
| DuckLips engine (baritone) | /s/ -9.8 dB re /a/, peak 3.5 kHz, centroid 7.2 kHz; /sh/ -14.0 dB, peak 2.5 kHz, centroid 5.3 kHz; /sh/-/s/ -4.2 dB; RTF 0.065 for the whole voice |
| Membranes (59-141 modes, all below 610 Hz) | change the output by -40 to -45 dB and double the cost. Soft tissue is mass-controlled at 2-8 kHz |
| Reciprocity, tongue -> lower lip (culled) | mean 5.3 dB magnitude error (worst 15.6), 12.5 deg phase, 1-10 kHz |
| Cost | rigid RTF 2.6 (/s/, 103 patches) and 5.9 (/sh/, 155); 9.7 with membranes; 145 for the 576-patch duct. 1-D control 0.008. Cost grows as patches^2 x oversampling |

## Failed or rejected attempts

1. **REJECTED: "air only links surfaces that see each other."** With a
   line-of-sight rule, an open end never gets its pressure-release
   reflection. A baffle cannot help, because the interior cannot see it.
   The open-duct resonance converges to a value 20% low. Correct radiation
   needs the full signed Kirchhoff kernel, including couplings from behind
   and through solids, plus a closing surface (the face).
2. **REJECTED (for this geometry): explicit full-kernel marching.** The
   full kernel is right for the baffled duct (+/-3%). In the mouth it
   diverges: thin teeth slabs and non-convex corners give strong back-face
   near-field couplings with sub-sample delays, which are clamped to one
   sample. This is the classic explicit-MOT late-time instability.
3. **REJECTED: two-point time averaging as the stabilizer.** It is stable,
   but the dissipation is O(rate), not O(hops). It erases every resonance.
4. **HYPOTHESIS (untested): an implicit near-field step.** Solve the
   sub-sample couplings simultaneously each step, which is the "little
   feedback", or use Galerkin/Burton-Miller. That is the standard route to
   stable, low-loss TDBEM. It costs a sparse solve per sample on top of
   P^2 taps.
5. **The membrane eigen-decomposition does not matter for sibilants.** The
   modes lie below 610 Hz, and the tissue impedance is about 1000 times
   rho c at 4 kHz. The walls act as rigid reflectors.
6. Setup: `m4build` deletes the remote `out/`, and a DuckLips run there was
   lost mid-session. Outputs now go to `build/out/`.

## House style (author review only)

There is no `auto`, no lambdas and no structured bindings. `step()` does not
allocate. Known deviations: ternaries in diagnostic `printf` arguments, and
compound conditions in a few `for` headers. `inject_acceleration_at`
reallocates the active tap lists; it is a setup operation.
