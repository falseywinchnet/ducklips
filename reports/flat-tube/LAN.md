# LAN access — 2026-10-10

Requested access from the MacBook Neo over the same LAN. The server now listens
on all IPv4 interfaces on port 8765. Both http://192.168.1.18:8765/ and
http://Joshuahs-Mac-mini.local:8765/ returned HTTP 200 from this host. The
listener was verified as *:8765. No firewall settings were changed.

An HTTP LAN origin cannot use the previous AudioWorklet-only startup. Extracted
the existing renderer into audio-renderer.mjs and retained a thin AudioWorklet
wrapper for localhost/secure origins. The LAN path runs that same renderer in
a module Worker and schedules four 2048-frame PCM buffers with Web Audio.
It does not use ScriptProcessorNode, certificates, insecure-origin browser
flags, or a different voice source. Physical engine.mjs is unchanged.

Also replaced the page's use of secure-context-only structuredClone for plain
JSON designs, so recording/edit provenance works over HTTP. Recording receipts
now freeze the shared renderer and both audio transport paths.

Verification:

- Twelve tests passed, including exact equality of monitored PCM and raw
  captures across 128-frame and 2048-frame output blocks, and raw capture while
  monitor level is zero.
- A direct comparison against the preceding commit's AudioWorklet produced
  bit-identical PCM after extraction into the shared renderer.
- Browser opened the actual LAN origin, started pressure, ran the engine and
  captured a fresh two-second WAV without errors/warnings. The same-origin
  enhanced upload was accepted and generated a new receipt with six source
  files. The browser core benchmark displayed 10.7x.
- The host was tested through its LAN URL; no claim is made that the MacBook
  Neo itself was remotely operated or its Wi-Fi route independently tested.

The LAN queue adds about 0.2 seconds of playback latency. It remains susceptible
to sufficiently long UI-thread scheduling stalls. The host must remain awake
with the server running. Device-local saved geometry is associated with the
browser origin; localhost and the LAN URL have separate saved designs.
