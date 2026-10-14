#!/bin/bash
# 10-40 Hz modulation and jitter of sentences 1+2 for option sets (each set: label=options).
cd /tmp/es; S=~/ducklips/tract/build/sheet_speak
for pair in "$@"; do
  name=${pair%%=*}; opts=${pair#*=}
  $S v1u.score vi1.wav $opts >/dev/null 2>&1; $S v2u.score vi2.wav $opts >/dev/null 2>&1
  ~/miniforge3/bin/python3 -c "
import numpy
from scipy.io import wavfile
a=wavfile.read('/tmp/es/vi1.wav')[1].astype(float); b=wavfile.read('/tmp/es/vi2.wav')[1].astype(float)
wavfile.write('/tmp/es/vi.wav',48000,numpy.concatenate([a,numpy.zeros(24000),b]).astype(numpy.int16))"
  ~/miniforge3/bin/python3 /tmp/vibrator.py "$name" /tmp/es/vi.wav 2>/dev/null | sed "s/stretches/[$opts] stretches/"
done
rm -f vi1.wav vi2.wav vi.wav
