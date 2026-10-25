"""A recovered word's answer polished by block rounds (sequence_recover's
FAST_ALGO blocks: the times alone, then the tubes and velum alone, round
after round while a round gains): from out/reversible/sequence_<label>.npz,
the same known answer rebuilt from its specs; the polished answer saved
beside it (_polished) with its round trip and wave file.

usage (on the Mini): POOL=8 RENDERER=incremental polish_word.py PHONE,PHONE,...
(the specs the word was recovered with)"""
import os
import sys
from multiprocessing import Pool
from typing import List

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import driver_ledger
import posture_recover
import reverse_core
import sequence_recover as sr


def label_of(specs: List[str]) -> str:
    """The saved files' label (sequence_recover.main's)."""
    return "_".join(spec.replace("@", "-").replace("^", "v").replace("~", "n").replace("!", "x")
                    for spec in specs).lower()


def round_trip(word: sr.Word, parameters: numpy.ndarray, wave: str) -> float:
    """The answer's sound against the true one's (dB), its wave file kept."""
    found_tubes, found_arrivals, found_leavings = word.shape.split(parameters)
    found_velum = word.shape.velum(parameters)
    times, frames = sr.keyframes(found_tubes, found_arrivals, found_leavings)
    port = sr.velum_keyframes(found_arrivals, found_leavings, found_velum) if word.nasals else None
    _, again = posture_recover.render((found_tubes[0], word.length), word.seconds, keyframes=(times, frames),
                                      seal=word.seal, velum=port, jet=word.jet, wave=wave,
                                      jet_gate=sr.jet_keyframes(found_arrivals, found_leavings, word.fricatives),
                                      folds=sr.fold_keyframes(found_arrivals, found_leavings, word.voiceless))
    stretch = slice(word.start, word.start + word.size)
    return reverse_core.snr(word.flow_out[stretch], again[stretch])


def main() -> None:
    specs = sys.argv[1].split(",")
    word = sr.build_word(specs)
    sr.configure_targets(word)
    label = label_of(specs)
    saved = numpy.load(os.path.join(sr.OUT, "sequence_%s.npz" % label))
    parameters = word.shape.join(saved["tubes"], saved["arrivals"], saved["leavings"], saved["velum"])
    print("%s: from the saved answer, round trip %.1f dB" % (
        " ".join(specs), round_trip(word, parameters, os.path.join(sr.OUT, "sequence_%s_before_polish.wav" % label))))
    sys.stdout.flush()
    report = sr.BlockReport(word.shape, word.length, word.true_tubes, word.true_arrivals, word.true_leavings,
                            word.nasals)
    with Pool(driver_ledger.POOL) as pool:
        parameters = sr.block_rounds(pool, word.shape, parameters, word.whole, word.reference, report)
    found_tubes, found_arrivals, found_leavings = word.shape.split(parameters)
    print("polished: round trip %.1f dB" % round_trip(
        word, parameters, os.path.join(sr.OUT, "sequence_%s_recovered_polished.wav" % label)))
    numpy.savez(os.path.join(sr.OUT, "sequence_%s_polished.npz" % label), tubes=found_tubes, arrivals=found_arrivals,
                leavings=found_leavings, velum=word.shape.velum(parameters), length=word.length,
                true_tubes=word.true_tubes, true_arrivals=word.true_arrivals, true_leavings=word.true_leavings,
                true_velum=word.true_velum)


if __name__ == "__main__":
    main()
