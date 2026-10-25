#!/bin/bash
# A voice profile, fitted to one reference speaker in rounds (owner,
# 2026-10-18: "use profiles and then to refine the profiles ... so that you can
# fit genders and depths"). Each round:
#   1. the engine's tree (profile_fit.py: source, noise, variation, tract,
#      tongue, cycle);
#   2. the vowel postures to the speaker's own formants (organ_vopt3.py with
#      TARGETS, the VocalTractLab-style tongue's geometry);
#   3. the planner (profile_fit.py: plan_pitch, plan_timing, plan_consonants);
#   4. the phones' mean durations to the speaker's labels (duration_calibrate.py,
#      closed form, planning only).
# Every stage starts from the last stage's best; a stage that improves nothing
# hands its input on. Renders through build_dev/speak.
# usage: profile_cycle.sh NAME REFDIR TRAIN HELDOUT TARGETS.json PARAMS ORGANS PLAN ROUNDS GENERATIONS [PROBE_EXTRA]
set -u
NAME=$1; REFDIR=$2; TRAIN=$3; HELDOUT=$4; TARGETS=$5; PARAMS=$6; ORGANS=$7; PLAN=$8; ROUNDS=$9; GENERATIONS=${10}
EXTRA_PROBE="anat.vtl=1 ${11:-}"
DURATION_IDS=${DURATION_IDS:-a0001-a0100}
HERE=$(cd "$(dirname "$0")" && pwd)
PY=$HOME/miniforge3/bin/python3
export SPEAK_BIN=$HOME/ducklips/tract/build_dev/speak
# One thread per worker: nine workers each spinning up the numerical
# libraries' threads on every core made a generation take 7 minutes.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
OUT=/tmp/es/profile/$NAME
mkdir -p "$OUT"
cp "$PARAMS" "$OUT/current.params"; cp "$ORGANS" "$OUT/current.organs"; cp "$PLAN" "$OUT/current.plan"
cd "$HERE"
for round in $(seq 1 "$ROUNDS"); do
  echo "=== $NAME round $round: tree $(date)"
  stage=${NAME}_r${round}_tree
  $PY profile_fit.py "$stage" "$REFDIR" "$TRAIN" "$HELDOUT" "$OUT/current.params" "$OUT/current.organs" "$OUT/current.plan" \
    "$GENERATIONS" source,noise,variation,tract,tongue,cycle 2>/dev/null | grep -v Warn
  [ -f "/tmp/es/profile/$stage/best.params" ] && cp "/tmp/es/profile/$stage/best.params" "$OUT/current.params"
  echo "=== $NAME round $round: vowels $(date)"
  (cd "$HERE" && PROBE_EXTRA="$EXTRA_PROBE" TARGETS="$TARGETS" SPEAK="$SPEAK_BIN" VOWEL_PITCH="${VOWEL_PITCH:-115}" \
     SPEAK_ENV="KINEMATIC=1 PARAMS=$OUT/current.params" MAXFEV=200 \
     $PY organ_vopt3.py "$OUT/current.organs" "$OUT/current.organs" "$OUT/vowels_r$round.organs" 6 \
     aa ae ao ah eh ey ow ih uh ax iy uw 2>/dev/null | grep -v Warn)
  [ -s "$OUT/vowels_r$round.organs" ] && cp "$OUT/vowels_r$round.organs" "$OUT/current.organs"
  echo "=== $NAME round $round: plan $(date)"
  stage=${NAME}_r${round}_plan
  $PY profile_fit.py "$stage" "$REFDIR" "$TRAIN" "$HELDOUT" "$OUT/current.params" "$OUT/current.organs" "$OUT/current.plan" \
    "$GENERATIONS" plan_pitch,plan_timing,plan_consonants 2>/dev/null | grep -v Warn
  [ -f "/tmp/es/profile/$stage/best.plan" ] && cp "/tmp/es/profile/$stage/best.plan" "$OUT/current.plan"
  echo "=== $NAME round $round: durations $(date)"
  $PY duration_calibrate.py "$REFDIR" "$DURATION_IDS" "$OUT/current.organs" "$OUT/current.plan" "$OUT/durations_r$round.plan" 3 \
    2>/dev/null | grep -v Warn | head -3
  [ -s "$OUT/durations_r$round.plan" ] && cp "$OUT/durations_r$round.plan" "$OUT/current.plan"
done
echo "=== $NAME done $(date)"
