# ducklips — tube workshop

A speech-engine experiment built around a connected 2D material sheet that
wraps around a directly editable tube. The side cutaway is its editing view.
Open the two-wall cutaway, insert volume blocks, move its supports, and inspect
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
- **Insert volume:** choose a block volume and click in the tube. Drag it to
  relocate the volume; select it and press Delete to remove it.
- Select a band to adjust material, mass, elasticity or support stiffness.
  The second band is the initially tight vocal band; the first is its return
  chamber. All **384 supported sheet modes** update when geometry, material or
  support changes. The frequency-ordered mode selector previews the complete
  unrolled sheet, including its joined circumferential edges.
- **Start pressure** connects the reservoir and outlet monitor. Set pressure,
  slight pressure wander and monitor level in the sidebar.
- **Record 2 seconds** captures the current engine and automatically computes
  its complete enhanced depiction locally. Save/open keeps tube designs as JSON.
- Connections and mode previews are optional inspection overlays.

This is an experimental instrument, **not an accepted human voice**. Its present
acoustic mesh, modal tissue and contact approximation have important physical
limits, especially glottal jet formation. Read the [model and implementation
notes](web/tube/MODEL.md) and [connected-sheet evidence](reports/flat-tube/SHEET.md).

```sh
node --test web/tube/*.test.mjs
```

The previous engine was rejected and is no longer the project entry point.
Its last pre-restart state is preserved at `archive/pre-flat-tube-20261010`
(`36d670d`), with its documentation in [README.legacy.md](README.legacy.md).
The new web engine imports none of its synthesis code. Historical sources and
recordings remain recoverable; their listening results do not score this tube.
