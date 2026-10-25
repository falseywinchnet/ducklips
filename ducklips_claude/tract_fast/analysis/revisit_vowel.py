"""A recovered word's vowel searched again in its sequence as now known (the
other phones and the times as recovered), from the hint's survey and its
basins as at first; then the whole word's joint and block rounds.
(Alfa's last vowel came out 0.89 cm^2 off from a first search whose
context was the hinted times and a prefix not yet joined; at the recovered
times the true vowel scored -83.7 dB in its window against -74.9 for the
found one.)

usage (on the Mini): POOL=8 RENDERER=incremental revisit_vowel.py PHONE,PHONE,... INDEX[,INDEX...]
(from the latest saved answer: _revisited if there is one, else the
recovered one; each vowel searched again in turn, then the joints and
block rounds once)"""
import os
import sys
from multiprocessing import Pool
from typing import List

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import driver_ledger
import tube_recover
import sequence_recover as sr
import polish_word


def main() -> None:
    specs: List[str] = sys.argv[1].split(",")
    indices = [int(value) for value in sys.argv[2].split(",")]
    word = sr.build_word(specs)
    sr.configure_targets(word)
    label = polish_word.label_of(specs)
    latest = os.path.join(sr.OUT, "sequence_%s_revisited.npz" % label)
    if not os.path.exists(latest):
        latest = os.path.join(sr.OUT, "sequence_%s.npz" % label)
    saved = numpy.load(latest)
    tubes = [numpy.array(tube, float) for tube in saved["tubes"]]
    arrivals = numpy.array(saved["arrivals"], float)
    leavings = numpy.array(saved["leavings"], float)
    velum = numpy.array(saved["velum"], float)
    hints = sr.table(os.path.join(sr.DATA, "organs_vtl_tongue.txt"))
    port_areas = velum if word.nasals else None
    print("%s: from %s (round trip %.1f dB)" % (" ".join(specs), os.path.basename(latest), polish_word.round_trip(
        word, word.shape.join(numpy.array(tubes), arrivals, leavings, velum),
        os.path.join(sr.OUT, "sequence_%s_before_revisit.wav" % label))))
    sys.stdout.flush()
    with Pool(driver_ledger.POOL) as pool:
        for index in indices:
            name = specs[index].split("@")[0]
            begin, span = sr.window(arrivals[index] + sr.SETTLE, leavings[index] - sr.CLEAR)
            leading = {k: (word.stops[k], word.narrowings.get(k, sr.STOP_AREA)) for k in range(index)
                       if k in word.stops}
            if len(leading) == index and index > 0:
                # Only consonants before it: they are made from each candidate
                # (as the first search makes a word's leading consonants).
                context = sr.PrefixContext([], arrivals, leavings, word.length, word.seal, port_areas, word.jet,
                                           word.voiceless, leading, word.fricatives)
            else:
                context = sr.PrefixContext(tubes[:index], arrivals, leavings, word.length, word.seal, port_areas,
                                           word.jet, word.voiceless, fricatives=word.fricatives)
            print("phone %d (%s) searched again in its recovered sequence, window %.3f-%.3f s" % (
                index, name, begin / sr.STEP_RATE, (begin + span) / sr.STEP_RATE))
            sys.stdout.flush()
            value, found = tube_recover.recover(pool, hints[name], word.flow_out, begin, span,
                                                (word.true_tubes[index], word.length), word.count, context)
            print("%s again: %.1f dB, tube rms %.4f" % (
                name, value, tube_recover.tube_rms(found, (word.true_tubes[index], word.length))))
            sys.stdout.flush()
            tubes[index] = numpy.array(found[0], float)
            # Phones made from this vowel follow it.
            for k in range(len(specs)):
                if k in word.stops and sr.vowel_before(k, word.stops) == index:
                    tubes[k] = tubes[index].copy()
                    tubes[k][word.stops[k]] = word.narrowings.get(k, sr.STOP_AREA)
        parameters = word.shape.join(numpy.array(tubes), arrivals, leavings, velum)
        print("before the joint: %.1f dB" % sr.level(word.shape.residual(parameters, word.whole), word.reference))
        for rung in range(4):
            parameters, used = tube_recover.stage(pool, word.shape, parameters, word.whole, 30)
            print("joint %d: %d evaluations; %.1f dB" % (
                rung, used, sr.level(word.shape.residual(parameters, word.whole), word.reference)))
            sys.stdout.flush()
        report = sr.BlockReport(word.shape, word.length, word.true_tubes, word.true_arrivals, word.true_leavings,
                                word.nasals)
        parameters = sr.block_rounds(pool, word.shape, parameters, word.whole, word.reference, report)
    found_tubes, found_arrivals, found_leavings = word.shape.split(parameters)
    print("revisited: round trip %.1f dB" % polish_word.round_trip(
        word, parameters, os.path.join(sr.OUT, "sequence_%s_recovered_revisited.wav" % label)))
    numpy.savez(os.path.join(sr.OUT, "sequence_%s_revisited.npz" % label), tubes=found_tubes,
                arrivals=found_arrivals, leavings=found_leavings, velum=word.shape.velum(parameters),
                length=word.length, true_tubes=word.true_tubes, true_arrivals=word.true_arrivals,
                true_leavings=word.true_leavings, true_velum=word.true_velum)


if __name__ == "__main__":
    main()
