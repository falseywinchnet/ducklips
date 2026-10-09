#!/bin/sh
# Renders a vowel/fricative calibration sequence and prints the rms of each
# segment's steady middle. Run from research/ducklips after building.
# usage: tools/measure_segments.sh [voice]
set -eu
voice=${1:-baritone}
./build/ducklips_say --voice "$voice" --phonemes "AA1:0.4 S:0.25 AA1:0.4 SH:0.25 AA1:0.4 F:0.2 AA1:0.4 HH:0.15 AA1:0.4" \
    --out out/calibration_"$voice".wav | grep -E "peak"
for window in "0.70 0.95 AA" "1.04 1.21 S" "1.69 1.86 SH" "2.34 2.46 F" "2.55 2.85 AA" "2.94 3.02 HH"; do
    set -- $window
    printf '%-3s ' "$3"
    ./build/ducklips_analyze out/calibration_"$voice".wav --from "$1" --to "$2" | head -1 | sed 's/.*peak/peak/'
done
