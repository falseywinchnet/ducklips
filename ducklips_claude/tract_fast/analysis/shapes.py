"""Shapes from sound: each vowel's posture shared across a voice's words.

The owner's goal (2026-10-24): the invertible engine is the instrument that
recovers speech geometry for the fast engine; the postures must be
recovered, not given. Steps 10-15 (NOTES_pyramid.md) measured the limits on
one word under the free time field: the equation error is a plateau far
from the truth, a held vowel's sound pins 2-4 combinations of the body's 9
parameters, and from the hints (about 1 cm2 rms off) the articulation did
not come back. Two things differ here:
  - the timing is the traverse family (puppet.TraverseProblem): a midpoint a
    transition, a lag and a move time an articulator; it cannot trade
    freely against shape;
  - a vowel's posture is one set of body parameters shared by every word
    and every transition it appears in, heard through many moves, not one
    hold.

Unknowns: the voice's (voice.VoiceProblem, traverse mode, groups) plus,
per vowel, the free body fields (posture_recover.FIELDS less FROZEN). Truth:
the words' postures from profile_bdl5 (as build_word makes them). Start:
truth + FRACTION (hint - truth), the hint table organs_vtl_tongue. A stop's
tube is its vowel's with the seal node closed (pyramid_word.targets).

Usage: shapes.py WORD;WORD;...  FRACTION  [traverse-known]
"""

import concurrent.futures
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy

import fields_word
import posture_recover
import puppet
import pyramid_word
import sequence_recover
import voice

FROZEN = ()
# No field frozen. Step 13's near-nulls are combinations (/ah/'s lip
# aperture against width, 4e-8; protrusion 2.2e-7), not single fields:
# frozen at the start, protrusion and width capped the reachable level
# (IY-AA-UW and Papa, the rest at the truth) at -38.7 dB from 20% of the way
# to the hints and -22.7 dB from the hints; width alone -36.5 / -21.9,
# protrusion alone -39.9 / -23.5; none frozen -250.1 (step 60). The damped
# solve handles a near-null combination; freezing a heard field does not.
FREE = tuple(field for field in posture_recover.FIELDS if field not in FROZEN)
MEAN_MOVE = puppet.DURATION_START
GAUGE = "pharynx"   # the articulator whose lag is the time reference (held at zero)
# A posture's difference step at a rung is the same fraction of its field's
# spread over the vowel table as the rung's time secant is of the plan's mean
# move (10 ms of 65 ms: 15%; one step: 0.016%): one resolution for shape and
# time, no constant of its own.


def vowel_names(contexts: Sequence[pyramid_word.Context]) -> List[str]:
    """The vowels the words contain (phones that are not stops or
    narrowings), in first-seen order."""
    names = []
    for context in contexts:
        word = context.word
        for k in range(word.phones):
            if k not in word.stops and word.names[k] not in names:
                names.append(word.names[k])
    return names


def field_spreads(hints: Dict[str, Dict[str, float]]) -> numpy.ndarray:
    """Each free field's spread (max - min) over the hint table's vowels."""
    spreads = numpy.empty(len(FREE))
    for position, field in enumerate(FREE):
        values = numpy.array([hints[name][field] for name in hints])
        spreads[position] = float(values.max() - values.min())
    return spreads


class TubeCache:
    """posture_recover.nodes_full per (posture, node count), kept: a residual
    asks for the same few vowels many times. Full precision from the
    process's one kept probe (step 62: four-decimal areas were 1e-4 cm2
    quanta, the size of the finest posture steps; and a probe process a
    posture forked a threaded worker, which macOS aborted)."""

    def __init__(self) -> None:
        self.tubes: Dict[Tuple, numpy.ndarray] = {}

    def tube(self, posture: Dict[str, float], count: int) -> numpy.ndarray:
        key = (tuple(posture[field] for field in posture_recover.FIELDS), count)
        found = self.tubes.get(key)
        if found is None:
            found = posture_recover.nodes_full(posture, count)[0]
            self.tubes[key] = found
        return found


class ShapesProblem(voice.VoiceProblem):
    """The voice's vector followed by each vowel's free body fields."""

    def __init__(self, contexts: Sequence[pyramid_word.Context], secant: float, cached: bool,
                 fraction: float) -> None:
        voice.VoiceProblem.__init__(self, contexts, secant, "traverse", cached, True)
        self.vowels = vowel_names(contexts)
        word = contexts[0].word
        self.truth_postures = {name: dict(word.profile[name]) for name in self.vowels}
        self.hint_postures = {name: dict(word.hints[name]) for name in self.vowels}
        self.fraction = fraction
        self.start_postures = {}
        for name in self.vowels:
            posture = {}
            for field in posture_recover.FIELDS:
                truth = self.truth_postures[name][field]
                posture[field] = truth + fraction * (self.hint_postures[name][field] - truth)
            self.start_postures[name] = posture
        self.spreads = field_spreads(word.hints)
        self.posture_offset = voice.VoiceProblem.size(self)
        self.tube_cache = TubeCache()
        # Per word, its vowels' indices into self.vowels.
        self.word_vowels = []
        for context in contexts:
            names = [context.word.names[k] for k in range(context.word.phones) if k not in context.word.stops]
            self.word_vowels.append(sorted(set(self.vowels.index(name) for name in names)))
        self.steps = self.steps_for(secant)

    def size(self) -> int:
        return voice.VoiceProblem.size(self) + len(getattr(self, "vowels", [])) * len(FREE)

    def postures_of(self, vector: numpy.ndarray) -> Dict[str, Dict[str, float]]:
        postures = {}
        for index, name in enumerate(self.vowels):
            posture = dict(self.start_postures[name])
            first = self.posture_offset + index * len(FREE)
            for position, field in enumerate(FREE):
                posture[field] = float(vector[first + position])
            postures[name] = posture
        return postures

    def tubes_for(self, index: int, postures: Dict[str, Dict[str, float]]) -> numpy.ndarray:
        """Word `index`'s phones' tubes from the postures (a stop's from its
        vowel, its seal node closed or narrowed, as build_word makes them)."""
        context = self.words[index].context
        word = context.word
        tubes = numpy.zeros_like(word.true_tubes)
        for k in range(word.phones):
            if k not in word.stops:
                tubes[k] = self.tube_cache.tube(postures[word.names[k]], context.count)
        for k in sorted(word.stops):
            tubes[k] = tubes[sequence_recover.vowel_before(k, word.stops)]
            tubes[k, word.stops[k]] = word.narrowings.get(k, sequence_recover.STOP_AREA)
        return tubes

    def word_vector(self, vector: numpy.ndarray, index: int) -> numpy.ndarray:
        own = voice.VoiceProblem.word_vector(self, vector, index)
        word = self.words[index]
        context = word.context
        word.rows = pyramid_word.targets(context, self.tubes_for(index, self.postures_of(vector)),
                                         context.word.true_velum)
        return own

    def columns_of(self, index: int, active: numpy.ndarray) -> numpy.ndarray:
        columns = voice.VoiceProblem.columns_of(self, index, active)
        extra = []
        for vowel in self.word_vowels[index]:
            first = self.posture_offset + vowel * len(FREE)
            extra.extend(range(first, first + len(FREE)))
        extra = numpy.array(extra, dtype=numpy.int64)
        return numpy.concatenate([columns, extra[active[extra]]])

    def steps_for(self, secant: float) -> numpy.ndarray:
        steps = numpy.full(self.size(), secant)
        relative = secant / MEAN_MOVE
        for index in range(len(self.vowels)):
            first = self.posture_offset + index * len(FREE)
            steps[first:first + len(FREE)] = relative * self.spreads
        return steps

    def start_vector(self) -> numpy.ndarray:
        vector = numpy.zeros(self.size())
        vector[:self.posture_offset] = voice.VoiceProblem.start_vector(self)[:self.posture_offset]
        for index, name in enumerate(self.vowels):
            first = self.posture_offset + index * len(FREE)
            vector[first:first + len(FREE)] = [self.start_postures[name][field] for field in FREE]
        return vector

    def true_vector(self) -> numpy.ndarray:
        vector = numpy.zeros(self.size())
        vector[:self.posture_offset] = voice.VoiceProblem.true_vector(self)[:self.posture_offset]
        for index, name in enumerate(self.vowels):
            first = self.posture_offset + index * len(FREE)
            vector[first:first + len(FREE)] = [self.truth_postures[name][field] for field in FREE]
        return vector


def posture_report(problem: ShapesProblem, vector: numpy.ndarray, label: str) -> None:
    """Each vowel's tube against the truth's (cm2 rms over the nodes of the
    first word that has it) and its body fields' errors."""
    postures = problem.postures_of(vector)
    parts = []
    for name in problem.vowels:
        index = next(i for i, vowels in enumerate(problem.word_vowels) if problem.vowels.index(name) in vowels)
        count = problem.words[index].context.count
        found = problem.tube_cache.tube(postures[name], count)
        want = problem.tube_cache.tube(problem.truth_postures[name], count)
        parts.append("%s %.3f cm2" % (name, float(numpy.sqrt(numpy.mean((found - want) ** 2)))))
    print("%s postures (tube rms off the truth): %s" % (label, "; ".join(parts)), flush=True)


def build_contexts(word_specs: Sequence[Sequence[str]]) -> List[pyramid_word.Context]:
    """voice.build_contexts' words (traverse truth), their true tubes
    rebuilt from the true postures at full precision first, so the truth
    is exactly in the family the solve reads through TubeCache."""
    cache = TubeCache()
    contexts = []
    for specs in word_specs:
        context = pyramid_word.Context(list(specs), ["areas"])
        word = context.word
        tubes = numpy.zeros_like(word.true_tubes)
        for k in range(word.phones):
            if k not in word.stops:
                tubes[k] = cache.tube(word.profile[word.names[k]], context.count)
        for k in sorted(word.stops):
            tubes[k] = tubes[sequence_recover.vowel_before(k, word.stops)]
            tubes[k, word.stops[k]] = word.narrowings.get(k, sequence_recover.STOP_AREA)
        word.true_tubes = tubes
        puppet.retime_truth_traverse(context, voice.true_lags_of(context), voice.true_extras_of(context, "traverse"))
        contexts.append(context)
    return contexts


def worker_start(word_specs: Sequence[Sequence[str]], secant: float, cached: bool, fraction: float) -> None:
    """The pool's problem: voice.word_block and voice.word_residual read voice._voice."""
    contexts = build_contexts(word_specs)
    voice._voice = ShapesProblem(contexts, secant, cached, fraction)


def solve(word_specs: Sequence[Sequence[str]], fraction: float, traverse_known: bool,
          cached: bool) -> Tuple[ShapesProblem, numpy.ndarray]:
    contexts = build_contexts(word_specs)
    # The main process only reports: no residual caches (2.1 GB warm, step 63).
    problem = ShapesProblem(contexts, puppet.SECANTS[0], False, fraction)
    vector = problem.start_vector()
    truth = problem.true_vector()
    frozen = numpy.zeros(problem.size(), dtype=bool)
    if traverse_known:
        vector[int(problem.offsets[-1]):problem.posture_offset] = truth[int(problem.offsets[-1]):problem.posture_offset]
        frozen[int(problem.offsets[-1]):problem.posture_offset] = True
    elif GAUGE in problem.shared:
        # The gauge: c added to every lag and taken from every midpoint
        # changes nothing heard. Free, the solve drifted along it (every raw
        # lag 3-7 ms high together, Alfa's events 38 ms off, step 63); one
        # articulator's lag is held at zero as the reference. A convention,
        # not a fitted value.
        reference = int(problem.offsets[-1]) + problem.shared.index(GAUGE)
        vector[reference] = 0.0
        frozen[reference] = True
    print("shapes: %d words, %d vowels (%s), %d unknowns (%d midpoints, %d traverse%s, %d posture); start %.0f%% "
          "of the way to the hints; the truth in the family %.1f dB" % (
              len(problem.words), len(problem.vowels), " ".join(problem.vowels), problem.size() - int(frozen.sum()),
              int(problem.offsets[-1]), problem.posture_offset - int(problem.offsets[-1]),
              " held at the truth" if traverse_known else "", len(problem.vowels) * len(FREE), 100.0 * fraction,
              pyramid_word.level_db(problem.residual(truth))), flush=True)
    posture_report(problem, vector, "start")
    workers = fields_word.POOL
    problem.workers = workers
    heard = numpy.zeros(problem.size(), dtype=bool)
    # Futures executors: a worker that dies raises BrokenProcessPool here
    # instead of leaving its task waited for (multiprocessing.Pool hung
    # 19.5 h); pinned, so each worker caches only its own words.
    with voice.PinnedPool(workers, worker_start, (word_specs, puppet.SECANTS[0], cached, fraction)) as pool:
        current = voice.residual_parallel(pool, problem, vector)
        voice.report(problem, vector, pyramid_word.level_db(current), "start", current)
        for secant in puppet.SECANTS:
            problem.secant = secant
            problem.steps = problem.steps_for(secant)
            print("secant %.3g s (postures %.3g of each field's spread)" % (secant, secant / MEAN_MOVE), flush=True)
            vector, heard = voice.levenberg_broyden(pool, problem, vector, frozen, fields_word.ITERATIONS)
            posture_report(problem, vector, "secant %.3g s" % secant)
    voice.lag_report(problem, vector, heard)
    return problem, vector


def main() -> None:
    word_specs = [item.split(",") for item in sys.argv[1].split(";")]
    fraction = float(sys.argv[2])
    options = sys.argv[3:]
    problem, vector = solve(word_specs, fraction, "traverse-known" in options, "plain" not in options)
    numpy.save(os.path.join(sequence_recover.SCRATCH, "shapes_%03d.npy" % int(round(100 * fraction))), vector)


if __name__ == "__main__":
    main()
