# ducklips — tube workshop

A speech-engine experiment built around connected tubes, each carrying a periodic
2D material sheet. The side cutaway is its editing view.
Open the two-wall cutaway, move its supports, and inspect
pressure-driven motion and outlet sound. Materials have mass, elasticity,
supported modes and separate solid propagation. Air travels bidirectionally
through an invisible crisscross network with FIR attenuation.

```sh
node web/tube/server.mjs
```

Open **http://127.0.0.1:8765/** on the host Mac. On another device on the same
LAN, use the LAN address printed by the server (currently
**http://192.168.1.18:8765/** on this Mac mini), or
**http://Joshuahs-Mac-mini.local:8765/**. Keep the host Mac awake and the server
running. No dependencies or build step are required.

LAN HTTP uses the same synthesis renderer in a Worker with scheduled Web Audio
buffers; localhost retains AudioWorklet. No certificate installation is needed.
LAN playback keeps about 0.2 seconds of audio buffered. To restrict access back
to this computer, start with `HOST=127.0.0.1 node web/tube/server.mjs`.

- **Move supports:** one diamond per section controls a ring around the tube.
  Drag down to tighten, up to expand, and horizontally to move the support.
  Both cutaway walls follow that ring. Arrow keys also work.
- **Tube configuration:** choose Open, Balanced, or Rounded to try different
  mouth and teeth/lip segment arrangements. The accepted nine-section below-mouth
  stage feeds the mouth, with sinus and under-tongue feedback branches. The
  below-mouth system and your exact six-section reference remain separate presets.
- **Choose a tube:** select Below-mouth, Mouth, Teeth and lips, Sinus pocket, or
  Under-tongue pocket to edit its supports. Purple markers locate the mouth
  branches. Pressure and tissue motion return through every connection.
- Select a band to adjust material, mass, elasticity or support stiffness.
  The fifth section of Below-mouth is the cord. The default network has
  **832 computed sheet modes**, including all 480 in the original source stage.
  Each stage's modes update when its geometry, material or support changes. The frequency-ordered mode selector previews the complete
  unrolled sheet, including its joined circumferential edges.
- **Whole-body elasticity:** 1–4 multiplies the material settings. Independent
  noise, sawtooth relaxation and bubble pulses vary every section by at most
  ±10%; this changes material response, not an imposed voice waveform.
- **Start pressure** connects the reservoir and outlet monitor. Set pressure,
  slight pressure wander and monitor gain in the sidebar. Playback gain is separate
  from the raw recording.
- **Record 2 seconds** captures the current engine and automatically computes
  its complete enhanced depiction locally. Save/open keeps the entire connected system as JSON, including every branch.
- Connections and mode previews are optional inspection overlays.

This is an experimental instrument, **not an accepted human voice**. Its present
acoustic mesh, modal tissue and contact approximation have important physical
limits, especially glottal jet formation. Read the [model and implementation
notes](web/tube/MODEL.md), [connected-tube experiments](reports/connected-tubes/README.md),
and [connected-sheet evidence](reports/flat-tube/SHEET.md).

```sh
node --test web/tube/*.test.mjs
```

The previous engine was rejected and is no longer the project entry point.
Its last pre-restart state is preserved at `archive/pre-flat-tube-20261010`
(`36d670d`), with its documentation in [README.legacy.md](README.legacy.md).
The new web engine imports none of its synthesis code. Historical sources and
recordings remain recoverable; their listening results do not score this tube.
