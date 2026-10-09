# Frozen surface port operator

The detailed speech airway is too slow to be the accepted renderer. In the
current frozen full-body probe it advances 4,988 directed paths at 192 kHz,
allocates 319,232 delay samples, and runs at about 0.20× realtime. Approximately
66% of the measured step time is in `gather()` (geometry/admittance updates,
propagation, and incoming-wave sums); the rest is predominantly scattering.
A faster representation must preserve the cavity response before it replaces
this engine.

`SurfaceNetwork::linear_snapshot()` and `tools/surface_ports.py` provide an
independent frequency-domain evaluation of the **same** network. It solves for
563 junction pressures rather than explicitly advancing all the stored delay
samples. This is an exact elimination of the frozen paths, not a modal truncation.
The matrix depends on frequency and still contains every path's delay/filter.
Its 563 pressure unknowns are therefore **not** a claim that the time-domain
system has only 563 states, or that the rendering speed target has been met.

## Elimination from pressure and flow

At a junction, an incident pressure wave `a` and outgoing wave `b` satisfy
`p = a + b` and outward volume flow `u = Y(b-a)`, where `Y` is the path's
characteristic admittance in m³/(Pa s). These are the standard acoustic
[scattering-junction relations](https://dsprelated.com/freebooks/pasp/Lossless_Scattering.html).

For a paired path between pieces i and j, let `h_ij(z)` and `h_ji(z)` be its
forward and reverse transfer functions, including propagation and loss. Then

```text
a_i = h_ji (p_j - a_j)
a_j = h_ij (p_i - a_i)
D   = 1 - h_ij h_ji

u_i = Y [(1 + h_ij h_ji) p_i - 2 h_ji p_j] / D
u_j = Y [(1 + h_ij h_ji) p_j - 2 h_ij p_i] / D
```

Each paired path thus contributes a two-by-two admittance block. A yielding
wall or radiation termination with reflection filter R contributes the shunt
admittance `Y(1-R)/(1+R)`. Summing these contributions gives `A(z) p = u_source`.
The observed radiated signal uses the same discrete flow derivative as the
engine: `sample_rate * (1-z^-1) * sum(radiation_admittance * aperture_pressure)`.
There is no spectral equalizer or fitted response in this construction.

The exported paths reproduce the actual discrete implementations:

- Axial: integer delay, first-order fractional all-pass, and boundary-loss shelf.
  Their cascaded numerator and denominator are exported without approximation.
- Transverse and floor paths: their actual linearly interpolated delay.
  Fixed power normalization converts back to the same pressure-wave relation.
- Walls: mass/spring/resistance reflection filters, including their direct term.
- Lips and nostrils: the actual radiation reflection filters and admittances.

Zero-admittance paths contribute no flow. The generic Python assembly can omit
completely disconnected pressure unknowns; none were omitted in the six
reported fixtures. It explicitly refuses a path-elimination pole rather than
silently regularizing it. A positive Laplace contour avoids lossless path poles
when comparing finite impulse records.

## Verification

`surface_response` exports a prepared frozen snapshot and independently advances
an impulse through the C++ time-domain implementation. The full-body fixture
copies the real `SurfaceSynthesizer` network at its initial schwa posture,
including the bent airway, transverse paths, nasal tract, side cavities, and
floor. A separate floor fixture opens two connections to exercise a genuine
loop. The lateral source excites one wall piece instead of a symmetric glottal
input, so the checks do not rely solely on plane-wave excitation.

The six 0.25-second probes use a 1e-10 m³/s one-sample flow input. They compare
202 frequency points from 48 Hz through 12 kHz, including complex pressure
responses and the radiated flow derivative. Both evaluations use
`z = exp((40 + j*2*pi*f)/sample_rate)`: multiplying the impulse by `exp(-40*t)`
matches that contour exactly and suppresses finite-record tails. It is an
analysis window, not acoustic damping added to the engine.

| Fixture / input | Pieces | Directed paths | Worst relative complex error |
|---|---:|---:|---:|
| Straight tube / glottis | 52 | 102 | 3.78e-10 |
| Solid transverse tube / glottis | 402 | 3,616 | 1.26e-9 |
| Bent transverse tube / glottis | 402 | 3,616 | 5.38e-10 |
| Bent tube with two floor openings / glottis | 411 | 3,760 | 4.93e-10 |
| Full body / glottis | 563 | 4,988 | 3.16e-6 |
| Full body / lateral piece | 563 | 4,988 | 3.76e-6 |

The largest discrepancies occur in the weak response at the almost-closed
nasal outlet. The relative-error denominator is floored at 1e-6 of that
channel's peak response, explicitly in the analysis script. This is agreement
with the existing model, not validation of its anatomy against a human.

All sampled path and termination real-power margins are nonnegative, and the
largest normalized linear-system residual is below 1.4e-12. Reversing the
full-body glottis/lateral source-observer pair agrees to 5.1e-14 relative error
in the direct frequency solution. The full-body 202-frequency calculation took
about 0.135 seconds here; this is an offline diagnostic duration, **not** a
speech rendering realtime factor. Timing includes Python sparse solves and
excludes the C++ snapshot preparation. C++ step timing includes per-sample
clock instrumentation; its 0.20× is not a new production speech benchmark.

Measured settings, binary hash, platform, response/snapshot hashes, per-channel
errors and timing are recorded in
[`reports/frozen_surface_ports.json`](../reports/frozen_surface_ports.json).
The recorded Git head precedes these uncommitted-at-measurement changes; the
binary hash identifies the executable actually measured. Release and
AddressSanitizer/UndefinedBehaviorSanitizer CTest runs each pass all three
suites (surface snapshot guards, numerical core, and streaming/configuration).

## Reproduce

```sh
cmake -S surface -B build/surface -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/surface
ctest --test-dir build/surface --output-on-failure
mkdir -p out/port_response
for fixture in tube solid bent floor body; do
  ./build/surface/surface_response --fixture "$fixture" \
    --out "out/port_response/$fixture" --duration 0.25
done
./build/surface/surface_response --fixture body --source lateral \
  --out out/port_response/body-lateral --duration 0.25
python tools/surface_ports.py \
  out/port_response/tube out/port_response/solid out/port_response/bent \
  out/port_response/floor out/port_response/body out/port_response/body-lateral \
  --out reports/frozen_surface_ports.json
```

The `.f64` files contain native-endian interleaved doubles, with pressure
channels and a final radiated-flow-derivative channel specified by the JSON.
Generated arrays stay outside Git. `--body bodies/reference.body` can configure
the full-body fixture. Running the analysis gate fails on complex-response
error above 1e-5, linear residual above 1e-9, or negative scaled passive-power
margin below -1e-10. It does not treat a failed reference match as a usable
reduction.

## Boundary and next work

The snapshot rejects moving/unsettled geometry, nonzero mean-flow kinetic
resistance, and membrane edge tension. Its current domain is the frozen,
zero-mean-flow linearization. It does not cover source nonlinearities, turbulent
series-pressure drives, articulator work, or time-varying state transfer.
No source or acoustic rendering algorithm was replaced by this change.

The subsequent [passive dynamic reduction](PASSIVE_REDUCTION.md) provides a
501-state frozen representation with all-port verification and a native modal
probe. Its coupled-fold case clears 10×, but the general multiport case still
fails. Carrying acoustic state through moving postures remains open. The current compact tree modes are still missing floor and
transverse physics; this reference prevents accepting them merely for being
fast. The full detailed-physics 10× gate and human-realism requirements remain
open.
