#!/bin/bash
# Open-loop pitch of a held /a/ against lung pressure: lung scales 0.6..1.6
# of 650 Pa, ear off, dry; median F0 over 0.9-2.1 s.
# usage: fp.sh [extra options]
for s in 0.6 0.8 1.0 1.3 1.6; do
  GAIN=0.0005 ~/ducklips/tract/build/sheet_speak /tmp/es/held_aa.score /tmp/es/fp_$s.wav --pitch-feedback 0 --room 0 --lung-scale $s "$@" >/dev/null 2>&1
  ~/miniforge3/bin/python3 - /tmp/es/fp_$s.wav $s <<'PY' 2>/dev/null
import sys, numpy
sys.path.insert(0, "/tmp")
from vw import load, pitch_at
x = load(sys.argv[1]); f = [pitch_at(x, int(t * 16000))[0] for t in numpy.arange(0.9, 2.1, 0.02)]
f = [v for v in f if v > 0]; lev = 20 * numpy.log10(numpy.std(x[int(0.9*16000):int(2.1*16000)]) + 1e-9)
print("lungs %4.0f Pa  F0 %6.1f Hz  (n %d)  level %.1f dB" % (650 * float(sys.argv[2]), numpy.median(f) if f else 0, len(f), lev))
PY
done
