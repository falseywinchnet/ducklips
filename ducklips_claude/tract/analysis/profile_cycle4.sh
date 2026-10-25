#!/bin/bash
# A voice profile in rounds, the class-credited way (2026-10-18). Each round:
#   1. class_fit.py: the vowels' organ postures, the consonants' plan shapes
#      and the frication tree settings, each ranked by its own phone classes'
#      cost; its best held-out mean (the start counts) goes on;
#   2. profile_fit.py on the source, noise, variation, cycle, plan_pitch and plan_breath groups (the
#      settings every class shares), kept only if its held-out combined beats
#      its start's (one render seed; v4 measure);
#   3. duration_calibrate.py: the phones' mean durations to the speaker's
#      labels (closed form, planning only).
# (profile_cycle3's formant refit is gone: matching LPC formants took bdl's
# /iy/ to a 0.036 cm^2 constriction and held-out 4.9 dB to 0.7.)
# Each stage's errors go to $OUT/r<round>_<stage>.err.
# usage: profile_cycle4.sh NAME REFDIR TRAIN HELDOUT PARAMS ORGANS PLAN ROUNDS CLASS_GENERATIONS TREE_GENERATIONS
set -u
NAME=$1; REFDIR=$2; TRAIN=$3; HELDOUT=$4; PARAMS=$5; ORGANS=$6; PLAN=$7; ROUNDS=$8; CLASS_GENERATIONS=$9; TREE_GENERATIONS=${10}
DURATION_IDS=${DURATION_IDS:-a0001-a0100}
HERE=$(cd "$(dirname "$0")" && pwd)
PY=$HOME/miniforge3/bin/python3
export SPEAK_BIN=$HOME/ducklips/tract/build_dev/speak
export SAY_BIN=$HOME/ducklips/voice/build_dev/say
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
OUT=/tmp/es/profile/$NAME
mkdir -p "$OUT"
cp "$PARAMS" "$OUT/current.params"; cp "$ORGANS" "$OUT/current.organs"; cp "$PLAN" "$OUT/current.plan"
cd "$HERE"
# The held-out combined from a "... heldout <report>" line.
heldout_of() {
  sed -nE 's/.*heldout .*combined (-?[0-9.]+).*/\1/p' | head -1
}
for round in $(seq 1 "$ROUNDS"); do
  echo "=== $NAME round $round: classes $(date)"
  stage=${NAME}_r${round}_class
  $PY class_fit.py "$stage" "$REFDIR" "$TRAIN" "$HELDOUT" "$OUT/current.params" "$OUT/current.organs" "$OUT/current.plan" \
    "$CLASS_GENERATIONS" 2>"$OUT/r${round}_class.err" | grep -v Warn
  if [ -s "/tmp/es/class/$stage/best_organs.txt" ]; then
    cp "/tmp/es/class/$stage/best_organs.txt" "$OUT/current.organs"
    cp "/tmp/es/class/$stage/best_plan.params" "$OUT/current.plan"
    cp "/tmp/es/class/$stage/best_voice.params" "$OUT/current.params"
  fi
  echo "=== $NAME round $round: tree $(date)"
  stage=${NAME}_r${round}_tree
  $PY profile_fit.py "$stage" "$REFDIR" "$TRAIN" "$HELDOUT" "$OUT/current.params" "$OUT/current.organs" "$OUT/current.plan" \
    "$TREE_GENERATIONS" source,noise,variation,cycle,plan_pitch,plan_breath 2>"$OUT/r${round}_tree.err" | grep -v Warn | tee "$OUT/r${round}_tree.log"
  if [ -f "/tmp/es/profile/$stage/best.params" ]; then
    before=$(grep "^start train" "$OUT/r${round}_tree.log" | heldout_of)
    after=$(head -1 "/tmp/es/profile/$stage/best.params" | heldout_of)
    echo "tree: held-out before $before after $after"
    if $PY -c "import sys; sys.exit(0 if float('$after') > float('$before') else 1)" 2>/dev/null; then
      cp "/tmp/es/profile/$stage/best.params" "$OUT/current.params"
      [ -f "/tmp/es/profile/$stage/best.plan" ] && cp "/tmp/es/profile/$stage/best.plan" "$OUT/current.plan"
      echo "tree: kept"
    else
      echo "tree: rejected"
    fi
  fi
  echo "=== $NAME round $round: durations $(date)"
  $PY duration_calibrate.py "$REFDIR" "$DURATION_IDS" "$OUT/current.organs" "$OUT/current.plan" "$OUT/durations_r$round.plan" 3 \
    2>"$OUT/r${round}_durations.err" | grep -v Warn | head -3
  [ -s "$OUT/durations_r$round.plan" ] && cp "$OUT/durations_r$round.plan" "$OUT/current.plan"
done
echo "=== $NAME done $(date)"
