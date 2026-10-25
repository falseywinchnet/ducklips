#!/bin/bash
# A voice profile in rounds, the class-credited way (2026-10-18); profile_cycle4
# plus the plan_motion and plan_nasal groups, line-buffered logs and a closed-
# form pitch stage. Each round:
#   1. class_fit.py: the vowels' organ postures, the consonants' plan shapes
#      and the frication tree settings, each ranked by its own phone classes'
#      cost; its best held-out mean (the start counts) goes on;
#   2. profile_fit.py on the settings every class shares, in two stages
#      (the source: source, noise, variation, cycle, tract, shape; the plan: plan_pitch,
#      plan_breath, plan_motion, plan_nasal), each kept only if its held-out
#      combined beats its start's (one render seed);
#   3. duration_calibrate.py: the phones' durations to the speaker's aligned
#      labels (closed form, planning only);
#   4. source_calibrate.py: the F0 percentiles' median, spread and low tail,
#      and the phase scatter (cycle lag, aspiration gain) (closed form, two
#      rendered passes; not for NATO).
# (profile_cycle3's formant refit is gone: matching LPC formants took bdl's
# /iy/ to a 0.036 cm^2 constriction and held-out 4.9 dB to 0.7.)
# Each stage's errors go to $OUT/r<round>_<stage>.err.
# usage: profile_cycle6.sh NAME REFDIR TRAIN HELDOUT PARAMS ORGANS PLAN ROUNDS CLASS_GENERATIONS TREE_GENERATIONS
set -u
NAME=$1; REFDIR=$2; TRAIN=$3; HELDOUT=$4; PARAMS=$5; ORGANS=$6; PLAN=$7; ROUNDS=$8; CLASS_GENERATIONS=$9; TREE_GENERATIONS=${10}
# Calibration sets spread across the corpus's recording sessions (bdl's F0
# on a0201-a0300 sits 0.7-1.6 st above a0101-a0120); the source set misses
# the held-out a0101-a0120.
DURATION_IDS=${DURATION_IDS:-a0001-a0593/3}
SOURCE_IDS=${SOURCE_IDS:-a0021-a0593/25}
HERE=$(cd "$(dirname "$0")" && pwd)
PY=$HOME/miniforge3/bin/python3
export SPEAK_BIN=$HOME/ducklips/tract/build_dev/speak
export SAY_BIN=$HOME/ducklips/voice/build_dev/say
# A refining step size for the shared stages (profile_fit's default 0.15 is
# for rough starts).
export PROFILE_SIGMA=${PROFILE_SIGMA:-0.08}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
OUT=/tmp/es/profile/$NAME
mkdir -p "$OUT"
cp "$PARAMS" "$OUT/current.params"; cp "$ORGANS" "$OUT/current.organs"; cp "$PLAN" "$OUT/current.plan"
cd "$HERE"
# The held-out combined from a "... heldout <report>" line.
heldout_of() {
  sed -nE 's/.*heldout .*combined (-?[0-9.]+).*/\1/p' | head -1
}
# A profile_fit stage on groups every class shares, kept only if its
# held-out combined beats its start's: shared_stage LABEL GROUPS.
shared_stage() {
  local label=$1 groups=$2
  echo "=== $NAME round $round: $label $(date)"
  local stage=${NAME}_r${round}_$label
  $PY profile_fit.py "$stage" "$REFDIR" "$TRAIN" "$HELDOUT" "$OUT/current.params" "$OUT/current.organs" "$OUT/current.plan" \
    "$TREE_GENERATIONS" "$groups" 2>"$OUT/r${round}_$label.err" | grep --line-buffered -v Warn | tee "$OUT/r${round}_$label.log"
  if [ -f "/tmp/es/profile/$stage/best.params" ]; then
    local before after
    before=$(grep "^start train" "$OUT/r${round}_$label.log" | heldout_of)
    after=$(head -1 "/tmp/es/profile/$stage/best.params" | heldout_of)
    echo "$label: held-out before $before after $after"
    if $PY -c "import sys; sys.exit(0 if float('$after') > float('$before') else 1)" 2>/dev/null; then
      cp "/tmp/es/profile/$stage/best.params" "$OUT/current.params"
      [ -f "/tmp/es/profile/$stage/best.plan" ] && cp "/tmp/es/profile/$stage/best.plan" "$OUT/current.plan"
      echo "$label: kept"
    else
      echo "$label: rejected"
    fi
  fi
}
# The closed-form stages: durations to the aligned labels, then the source (not
# for NATO). calibrate ROUND.
calibrate() {
  local round=$1
  echo "=== $NAME round $round: durations $(date)"
  $PY duration_calibrate.py "$REFDIR" "$DURATION_IDS" "$OUT/current.organs" "$OUT/current.plan" "$OUT/durations_r$round.plan" 3 \
    2>"$OUT/r${round}_durations.err" | grep --line-buffered -v Warn | head -3
  [ -s "$OUT/durations_r$round.plan" ] && cp "$OUT/durations_r$round.plan" "$OUT/current.plan"
  if [ "$REFDIR" != nato ]; then
    echo "=== $NAME round $round: source $(date)"
    $PY source_calibrate.py "$REFDIR" "$SOURCE_IDS" "$OUT/current.params" "$OUT/current.organs" "$OUT/current.plan" \
      "$OUT/source_r$round.params" "$OUT/source_r$round.plan" 2 2>"$OUT/r${round}_source.err" | grep --line-buffered -v Warn
    if [ -s "$OUT/source_r$round.plan" ] && [ -s "$OUT/source_r$round.params" ]; then
      cp "$OUT/source_r$round.plan" "$OUT/current.plan"
      cp "$OUT/source_r$round.params" "$OUT/current.params"
    fi
  fi
}
# (Calibrated once before the first round: rms5's start sat at held-out
# timing 1.1 through a whole round on old pooled durations.)
calibrate 0
for round in $(seq 1 "$ROUNDS"); do
  echo "=== $NAME round $round: classes $(date)"
  stage=${NAME}_r${round}_class
  $PY class_fit.py "$stage" "$REFDIR" "$TRAIN" "$HELDOUT" "$OUT/current.params" "$OUT/current.organs" "$OUT/current.plan" \
    "$CLASS_GENERATIONS" 2>"$OUT/r${round}_class.err" | grep --line-buffered -v Warn
  if [ -s "/tmp/es/class/$stage/best_organs.txt" ]; then
    cp "/tmp/es/class/$stage/best_organs.txt" "$OUT/current.organs"
    cp "/tmp/es/class/$stage/best_plan.params" "$OUT/current.plan"
    cp "/tmp/es/class/$stage/best_voice.params" "$OUT/current.params"
  fi
  shared_stage tree source,noise,variation,cycle,tract,shape
  shared_stage plan plan_pitch,plan_breath,plan_motion,plan_nasal
  calibrate $round
done
echo "=== $NAME done $(date)"
