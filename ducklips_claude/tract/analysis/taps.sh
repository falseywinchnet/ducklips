#!/bin/bash
cd /tmp/es/newbody
for v in "" "FRICATION=0" "ASPIRATION=0" "FRICATION=0 ASPIRATION=0" "BRANCHES=" ; do
  env ROOM=0 KINEMATIC=1 VOICEPOST=0.015,0.01 $v ~/ducklips/tract/build/speak k1.score tt.wav 73021 2>/dev/null
  echo "[$v] $(~/miniforge3/bin/python3 /tmp/thumps.py tt.wav | head -1)"
done
