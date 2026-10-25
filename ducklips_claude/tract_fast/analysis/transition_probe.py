"""One transition's timing grids alone, from the true earlier phones (the
first vowel is recovered exactly in these words): which leaving and arrival
the grids choose, and each rung's cost at the truth and at the choice.

usage: transition_probe.py WORD INDEX   (POOL workers; the model from the
environment: PITCH_EAR, ENGINE_GLIDE_PHASE)"""
import os
import sys
from multiprocessing import Pool

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import driver_ledger
import posture_recover
import sequence_recover
import tube_recover


def main() -> None:
    word = sequence_recover.build_word(sys.argv[1].split(","))
    sequence_recover.configure_targets(word)
    index = int(sys.argv[2])
    before = [numpy.array(word.true_tubes[k]) for k in range(index)]
    if index in word.stops:
        guess = numpy.array(before[sequence_recover.vowel_before(index, word.stops)])
        guess[word.stops[index]] = word.narrowings.get(index, sequence_recover.STOP_AREA)
    else:
        guess = posture_recover.nodes(word.hints[word.names[index]], word.count)[0]
    print("%s, transition %d; pitch ear %s, glide phase %s" % (
        " ".join(word.specs), index, posture_recover.PITCH_EAR_SETTING,
        os.environ.get("ENGINE_GLIDE_PHASE", "1")))
    with Pool(driver_ledger.POOL) as pool:
        arrivals, leavings = sequence_recover.refine_transition(
            pool, before, guess, word.hinted_arrivals.copy(), word.hinted_leavings.copy(), word.length, index,
            word.flow_out, word.true_arrivals, word.true_leavings,
            word.hinted_velum if word.nasals else None)
    tubes = numpy.array(before + [guess] * (word.phones - index))
    begin, span = sequence_recover.window(word.hinted_leavings[index - 1] - 0.02, word.hinted_arrivals[index] + 0.03)
    for band in (1500.0, 0.0):
        target = tube_recover.Target(word.flow_out, begin, span, band, 1.0, True)
        found = sequence_recover.timing_cost((tubes, arrivals, leavings, None, word.length, target))
        truth_arrivals = arrivals.copy()
        truth_leavings = leavings.copy()
        truth_arrivals[index] = word.true_arrivals[index]
        truth_leavings[index - 1] = word.true_leavings[index - 1]
        true_cost = sequence_recover.timing_cost((tubes, truth_arrivals, truth_leavings, None, word.length, target))
        print("band %.0f Hz: cost at the choice %.4g, at the true times %.4g" % (band, found, true_cost))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
