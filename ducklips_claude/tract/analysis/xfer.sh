#!/bin/bash
# Sheet transfer vs ideal for a held score under option sets.
# usage: xfer.sh SCORE AREAS "options" ...
cd /tmp/es; score=$1; areas=$2; shift 2
for o in "$@"; do
  TRACE=/tmp/es/tr_xf.txt TRACE_EVERY=1 GAIN=0.0005 ${SHEET_BIN:-~/ducklips/tract/build/sheet_speak} $score /tmp/es/hx.wav --room 0 $o >/dev/null 2>&1
  echo "== [$o]"; ~/miniforge3/bin/python3 /tmp/sheet_transfer.py tr_xf.txt hx.wav $areas 2>/dev/null | tail -9
done
rm -f /tmp/es/tr_xf.txt
