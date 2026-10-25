"""More final-joint rounds from a finished word's result, with the round
trip after each (a run's final joint restarts Levenberg-Marquardt from the
parameters every round, so continuing from the saved result is what a run
with more rounds would have done).

usage: joint_continue.py WORD ROUNDS [MIN_GAIN_DB]
(WORD as sequence_recover's; reads out/reversible/sequence_<label>.npz,
writes sequence_<label>_continued.npz)"""
import os
import sys
from multiprocessing import Pool

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import driver_ledger
import posture_recover
import reverse_core
import sequence_recover
import tube_recover


def round_trip(word: sequence_recover.Word, parameters: numpy.ndarray) -> float:
    """The recovered sequence's sound against the true one's (dB), as
    sequence_recover.main measures it."""
    tubes, arrivals, leavings = word.shape.split(parameters)
    velum = word.shape.velum(parameters)
    times, frames = sequence_recover.keyframes(tubes, arrivals, leavings)
    port = sequence_recover.velum_keyframes(arrivals, leavings, velum) if word.nasals else None
    _, again = posture_recover.render((tubes[0], word.length), word.seconds, keyframes=(times, frames), seal=word.seal,
                                      velum=port, jet=word.jet,
                                      jet_gate=sequence_recover.jet_keyframes(arrivals, leavings, word.fricatives),
                                      folds=sequence_recover.fold_keyframes(arrivals, leavings, word.voiceless))
    clip = slice(word.start, word.start + word.size)
    return float(reverse_core.snr(word.flow_out[clip], again[clip]))


def main() -> None:
    word = sequence_recover.build_word(sys.argv[1].split(","))
    sequence_recover.configure_targets(word)
    rounds = int(sys.argv[2])
    least_gain = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    label = "_".join(spec.replace("@", "-").replace("^", "v").replace("~", "n").replace("!", "x")
                     for spec in word.specs).lower()
    kept = numpy.load(os.path.join(sequence_recover.OUT, "sequence_%s.npz" % label))
    parameters = word.shape.join(kept["tubes"], kept["arrivals"], kept["leavings"], kept["velum"])
    print("from the saved result: equation error %.1f dB, round trip %.1f dB" % (
        sequence_recover.level(word.shape.residual(parameters, word.whole), word.reference),
        round_trip(word, parameters)))
    sys.stdout.flush()
    with Pool(driver_ledger.POOL) as pool:
        for index in range(rounds):
            with driver_ledger.stage("continued joint %d" % index):
                parameters, used, gain = tube_recover.stage_gain(pool, word.shape, parameters, word.whole, 30)
            tubes, arrivals, leavings = word.shape.split(parameters)
            rms = [tube_recover.tube_rms((tubes[k], word.length), (word.true_tubes[k], word.length))
                   for k in range(word.phones)]
            print("round %d: %d evaluations, gain %.2f dB; equation error %.1f dB; round trip %.1f dB; tubes rms max "
                  "%.4f" % (index, used, gain, sequence_recover.level(word.shape.residual(parameters, word.whole),
                                                                       word.reference),
                            round_trip(word, parameters), max(rms)))
            sys.stdout.flush()
            if gain < least_gain:
                break
    tubes, arrivals, leavings = word.shape.split(parameters)
    numpy.savez(os.path.join(sequence_recover.OUT, "sequence_%s_continued.npz" % label), tubes=tubes,
                arrivals=arrivals, leavings=leavings, velum=word.shape.velum(parameters))
    driver_ledger.LEDGER.report(depth=1)


if __name__ == "__main__":
    main()
