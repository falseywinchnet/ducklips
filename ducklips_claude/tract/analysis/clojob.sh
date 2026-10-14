#!/bin/bash
# one render + cycle dump: clojob.sh DIR KEY "ENV..."
d=$1; envs=$2; k=$3
env ROOM=0 CYCLES=$d/$k.cyc $envs ~/ducklips/tract/build/speak $d/$k.score $d/$k.wav 1 >/dev/null 2>&1
DUMP=$d/${k}_wm.txt ~/miniforge3/bin/python3 /tmp/wmatch.py $d/$k.wav >/dev/null 2>&1
