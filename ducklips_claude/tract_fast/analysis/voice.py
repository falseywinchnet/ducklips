"""The voice: per-articulator traverse shared across words (Phase D).

Several words of one synthetic speaker, each with its own event times
(puppet.PuppetProblem), and one traverse shared by articulator: node k of
the tube, the velum, the folds, each a lag and (with `extras`) an extra
move time. The equation errors of all the words are fitted at once; a
traverse the words' transitions do not move is unheard and stays at its
start. Cross-validation: the traverse from all words but one, the held-out
word's events then fitted with it fixed, its equation error compared
against no traverse and the true one.

Usage: voice.py WORD;WORD;...  [extras] [holdout INDEX]
  each WORD a comma-separated spec list as puppet.py takes.
"""

import multiprocessing
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy

import fields_word
import puppet
import pyramid_word
import sequence_recover

SECANTS = puppet.SECANTS


GROUPS = puppet.GROUPS   # the articulator groups (one definition, in puppet.py)
MODES = ("lags", "extras", "traverse")
# What the voice shares per articulator: its lag alone; its lag and an extra
# move time on the word's own (leaving, arrival) events; or its lag and its
# own move time, the words keeping one midpoint a transition
# (puppet.TraverseProblem, step 57).


def articulator_names(context: pyramid_word.Context, groups: bool = False) -> List[str]:
    """Each command column's articulator, shared across words by name. Named
    from the lips: the words' node counts differ (44 and 46, set by the
    first vowel's length), so a node counted from the glottis is a different
    place in each word (the voice's truth was then not in its family:
    -23.9 dB, step 51). With `groups`, one name an articulator group."""
    names = []
    for node in range(context.count):
        distance = context.count - 1 - node
        if groups:
            for name, first, last in GROUPS:
                if first <= distance < last:
                    names.append(name)
                    break
        else:
            names.append("lips-%d" % distance)
    if context.has_velum:
        names.append("velum")
    if context.has_folds:
        names.extend(["folds", "folds", "folds"])
    return names


def true_lags_of(context: pyramid_word.Context) -> numpy.ndarray:
    return puppet.friction_truth(context)


def true_extras_of(context: pyramid_word.Context, mode: str) -> numpy.ndarray:
    """The truth's second shared quantity a column: the extra move times
    (extras), the move times (traverse), or nothing (lags)."""
    columns = context.count + (1 if context.has_velum else 0) + (3 if context.has_folds else 0)
    if mode == "extras":
        return puppet.friction_truth_extras(context)
    if mode == "traverse":
        return puppet.traverse_truth_durations(context)
    return numpy.zeros(columns)


def second_block(mode: str) -> bool:
    """Whether the shared block carries a second quantity after the lags."""
    return mode != "lags"


class VoiceProblem:
    """vector = [word 0's events, word 1's events, ..., the shared lags (one
    an articulator name, in `shared` order), then (extras, traverse) the
    shared extra move times or move times in the same order]. In traverse
    mode a word's events are its transitions' midpoints."""

    def __init__(self, contexts: Sequence[pyramid_word.Context], secant: float, mode: str = "lags",
                 cached: bool = False, groups: bool = False) -> None:
        if mode not in MODES:
            raise ValueError("mode %r is not one of %s" % (mode, MODES))
        self.mode = mode
        self.groups = groups
        self.words: List[puppet.PuppetProblem] = []
        self.shared: List[str] = []
        self.maps: List[numpy.ndarray] = []   # per word, each column's index into the shared lags
        for context in contexts:
            names = articulator_names(context, groups)
            for name in names:
                if name not in self.shared:
                    self.shared.append(name)
            self.maps.append(numpy.array([self.shared.index(name) for name in names], dtype=numpy.int64))
            if mode == "traverse":
                word = puppet.TraverseProblem(context, secant=secant)
            else:
                word = puppet.PuppetProblem(context, secant=secant)
            if cached:
                puppet.attach_cache(word)
            self.words.append(word)
        self.offsets = numpy.cumsum([0] + [word.events() for word in self.words])
        self.secant = secant
        self.used = 1
        self.steps = numpy.full(self.size(), secant)

    def size(self) -> int:
        return int(self.offsets[-1]) + (2 if second_block(self.mode) else 1) * len(self.shared)

    def word_vector(self, vector: numpy.ndarray, index: int) -> numpy.ndarray:
        """Word `index`'s own vector (its events; its traverse set from the shared)."""
        word = self.words[index]
        word.lags = self.shared_lags(vector)[self.maps[index]]
        if self.mode == "extras":
            word.extras = self.shared_extras(vector)[self.maps[index]]
        elif self.mode == "traverse":
            word.durations = self.shared_extras(vector)[self.maps[index]]
        return vector[self.offsets[index]:self.offsets[index + 1]]

    def shared_lags(self, vector: numpy.ndarray) -> numpy.ndarray:
        return vector[self.offsets[-1]:self.offsets[-1] + len(self.shared)]

    def shared_extras(self, vector: numpy.ndarray) -> numpy.ndarray:
        return vector[self.offsets[-1] + len(self.shared):]

    def residual(self, vector: numpy.ndarray) -> numpy.ndarray:
        parts = []
        for index, word in enumerate(self.words):
            parts.append(word.residual(self.word_vector(vector, index)))
        return numpy.concatenate(parts)

    def blocks(self) -> List[slice]:
        """Each word's stretch of the residual."""
        found = []
        start = 0
        for word in self.words:
            size = word.context.word.size
            found.append(slice(start, start + size))
            start += size
        return found

    def start_vector(self) -> numpy.ndarray:
        """The hinted events, the lags at zero; in traverse mode the hinted
        midpoints and every move time at the plan's mean move."""
        vector = numpy.zeros(self.size())
        for index, word in enumerate(self.words):
            context_word = word.context.word
            if self.mode == "traverse":
                own = word.vector_of_midpoints(puppet.hinted_midpoints(context_word))
            else:
                own = word.vector_of(context_word.hinted_arrivals, context_word.hinted_leavings)
            vector[self.offsets[index]:self.offsets[index + 1]] = own[:word.events()]
        if self.mode == "traverse":
            self.shared_extras(vector)[:] = puppet.DURATION_START
        return vector

    def true_vector(self) -> numpy.ndarray:
        vector = numpy.zeros(self.size())
        for index, word in enumerate(self.words):
            context_word = word.context.word
            if self.mode == "traverse":
                own = word.vector_of_midpoints(puppet.true_midpoints(context_word))
            else:
                own = word.vector_of(context_word.true_arrivals, context_word.true_leavings)
            vector[self.offsets[index]:self.offsets[index + 1]] = own[:word.events()]
            self.shared_lags(vector)[self.maps[index]] = true_lags_of(word.context)
            if second_block(self.mode):
                self.shared_extras(vector)[self.maps[index]] = true_extras_of(word.context, self.mode)
        return vector


_voice: Optional[VoiceProblem] = None


def worker_start(word_specs: Sequence[Sequence[str]], secant: float, mode: str, cached: bool = False,
                 groups: bool = False) -> None:
    global _voice
    contexts = build_contexts(word_specs, mode)
    _voice = VoiceProblem(contexts, secant, mode, cached, groups)


def build_contexts(word_specs: Sequence[Sequence[str]], mode: str) -> List[pyramid_word.Context]:
    """The words, their truths retimed with the pattern's lags (and, with
    extras, its extra move times)."""
    contexts = []
    for specs in word_specs:
        context = pyramid_word.Context(list(specs), ["areas"])
        if mode == "traverse":
            puppet.retime_truth_traverse(context, true_lags_of(context), true_extras_of(context, mode))
        else:
            puppet.retime_truth(context, true_lags_of(context), true_extras_of(context, mode))
        contexts.append(context)
    return contexts


def worker_columns(task: Tuple[numpy.ndarray, List[int], numpy.ndarray]) -> numpy.ndarray:
    """Central-difference columns; an event's column touches its word's
    block only (the other blocks are not recomputed)."""
    # Three fields from a main process started before the steps travelled
    # with the task (it then reads the pool's own steps).
    base, indices, base_residual = task[0], task[1], task[2]
    steps = task[3] if len(task) > 3 else _voice.steps
    blocks = _voice.blocks()
    columns = numpy.zeros((len(base_residual), len(indices)))
    for position, index in enumerate(indices):
        step = steps[index]
        ahead = base.copy()
        ahead[index] += step
        behind = base.copy()
        behind[index] -= step
        owner = int(numpy.searchsorted(_voice.offsets, index, side="right") - 1)
        if index < _voice.offsets[-1]:
            word = _voice.words[owner]
            block = blocks[owner]
            columns[block, position] = (word.residual(_voice.word_vector(ahead, owner))
                                        - word.residual(_voice.word_vector(behind, owner))) / (2.0 * step)
        else:
            columns[:, position] = (_voice.residual(ahead) - _voice.residual(behind)) / (2.0 * step)
    return columns


def jacobian(pool, problem: VoiceProblem, vector: numpy.ndarray, current: numpy.ndarray,
             frozen: numpy.ndarray) -> numpy.ndarray:
    active = numpy.flatnonzero(~frozen)
    chunks = [chunk for chunk in numpy.array_split(active, fields_word.POOL * 2) if len(chunk) > 0]
    tasks = [(vector, list(chunk), current, problem.steps) for chunk in chunks]
    matrix = numpy.zeros((len(current), len(vector)))
    if tasks:
        matrix[:, active] = numpy.concatenate(pool.map(worker_columns, tasks), axis=1)
    return matrix


def report(problem: VoiceProblem, vector: numpy.ndarray, cost: float, label: str) -> None:
    truth = problem.true_vector()
    parts = []
    for index, word in enumerate(problem.words):
        own_truth = problem.word_vector(truth, index).copy()
        # word_vector sets the word's traverse from the vector it reads: the
        # vector's own last, so the block's level is the vector's.
        own = problem.word_vector(vector, index)
        block_cost = pyramid_word.level_db(word.residual(own))
        parts.append("%s %.1f dB, events %.3f ms max off" % (
            "/".join(word.context.word.names), block_cost, 1000.0 * numpy.abs(own - own_truth).max()))
    traverse = vector[problem.offsets[-1]:]
    true_traverse = truth[problem.offsets[-1]:]
    print("%s: all %.1f dB; %s; traverse off the truth %.3f ms max" % (
        label, cost, "; ".join(parts), 1000.0 * numpy.abs(traverse - true_traverse).max()), flush=True)


def lag_report(problem: VoiceProblem, vector: numpy.ndarray, heard: numpy.ndarray) -> None:
    """The traverse against the truth, raw and in the pharynx's gauge
    (adding c to every lag and taking it from every event changes nothing,
    step 54)."""
    truth = problem.true_vector()
    if "pharynx" in problem.shared:
        lags = problem.shared_lags(vector)
        reference = lags[problem.shared.index("pharynx")]
        true_lags = problem.shared_lags(truth)
        true_reference = true_lags[problem.shared.index("pharynx")]
        print("lags relative to the pharynx (ms; truth): %s" % " ".join(
            "%s:%+.2f(%+.0f)" % (name, 1000.0 * (lags[j] - reference), 1000.0 * (true_lags[j] - true_reference))
            for j, name in enumerate(problem.shared) if name != "pharynx"), flush=True)
    kinds = [("lags", problem.shared_lags)]
    if second_block(problem.mode):
        kinds.append(("move times" if problem.mode == "traverse" else "extras", problem.shared_extras))
    for title, part in kinds:
        values = part(vector)
        true_values = part(truth)
        part_heard = part(heard.astype(numpy.float64)) > 0.5
        print("%s (ms; truth): %s" % (title, " ".join(
            "%s:%s(%+.0f)" % (name, "unheard" if not flag else "%+.2f" % (1000.0 * values[j]), 1000.0 * true_values[j])
            for j, (name, flag) in enumerate(zip(problem.shared, part_heard))
            if flag or true_values[j] != 0.0)), flush=True)


def levenberg(pool, problem: VoiceProblem, vector: numpy.ndarray, frozen: numpy.ndarray,
              iterations: int, quiet: bool = False) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """fields_word.levenberg's step (damped Gauss-Newton in the secant's
    units, Nielsen's damping, the agreement gate) over the voice; returns
    the vector and which unknowns the sound heard (nonzero columns). Quiet:
    no report a step (the inner solves of varpro)."""
    current = problem.residual(vector)
    cost = pyramid_word.level_db(current)
    heard = numpy.zeros(len(vector), dtype=bool)
    for iteration in range(iterations):
        matrix = jacobian(pool, problem, vector, current, frozen)
        matrix[:, ~numpy.all(numpy.isfinite(matrix), axis=0)] = 0.0
        norms = numpy.sqrt(numpy.sum(matrix ** 2, axis=0))
        heard = norms > 0.0
        scale = problem.steps.copy()
        mean = float(numpy.mean((norms[heard] * scale[heard]) ** 2)) if numpy.any(heard) else 0.0
        if mean <= 0.0:
            break
        scale /= numpy.sqrt(mean)
        scaled = matrix[:, heard] * scale[heard][None, :]
        left, singular, rows = numpy.linalg.svd(scaled, full_matrices=False)
        projected = left.T @ (-current)
        damping = pyramid_word.DAMPING_START
        first = damping
        growth = 2.0
        accepted = False
        while not accepted and damping <= pyramid_word.DAMPING_LIMIT * first:
            z = rows.T @ (singular / (singular ** 2 + damping) * projected)
            step = numpy.zeros(len(vector))
            step[heard] = scale[heard] * z
            trial = problem.residual(vector + step)
            predicted = float(numpy.sum(current ** 2) - numpy.sum((current + matrix @ step) ** 2))
            actual = float(numpy.sum(current ** 2) - numpy.sum(trial ** 2))
            ratio = actual / predicted if predicted > 0.0 else 0.0
            if actual > 0.0 and ratio >= fields_word.AGREEMENT:
                vector, current = vector + step, trial
                cost = pyramid_word.level_db(current)
                damping *= max(1.0 / 3.0, 1.0 - (2.0 * ratio - 1.0) ** 3)
                accepted = True
            else:
                damping *= growth
                growth *= 2.0
        if not quiet:
            report(problem, vector, cost, "step %d%s" % (iteration, "" if accepted else ": no damping lowered it"))
        if not accepted:
            break
    return vector, heard


def solve(word_specs: Sequence[Sequence[str]], mode: str = "lags",
          fixed_traverse: Optional[numpy.ndarray] = None, label: str = "voice",
          cached: bool = False, groups: bool = False) -> Tuple[VoiceProblem, numpy.ndarray]:
    """The voice solved over the secant rungs. With `fixed_traverse` (the
    shared block as the vector holds it) the traverse is held there and
    only the events move."""
    contexts = build_contexts(word_specs, mode)
    problem = VoiceProblem(contexts, SECANTS[0], mode, cached, groups)
    vector = problem.start_vector()
    frozen = numpy.zeros(problem.size(), dtype=bool)
    if fixed_traverse is not None:
        vector[problem.offsets[-1]:] = fixed_traverse
        frozen[problem.offsets[-1]:] = True
    truth = problem.true_vector()
    print("%s: %d words, %d unknowns (%d events, %d shared); the truth in the family %.1f dB" % (
        label, len(problem.words), problem.size(), int(problem.offsets[-1]), problem.size() - int(problem.offsets[-1]),
        pyramid_word.level_db(problem.residual(truth))), flush=True)
    report(problem, vector, pyramid_word.level_db(problem.residual(vector)), "%s start" % label)
    heard = numpy.zeros(problem.size(), dtype=bool)
    for secant in SECANTS:
        problem.secant = secant
        problem.steps = numpy.full(problem.size(), secant)
        print("%s secant %.3g s" % (label, secant), flush=True)
        with multiprocessing.Pool(fields_word.POOL, initializer=worker_start,
                                  initargs=(word_specs, secant, mode, cached, groups)) as pool:
            vector, heard = levenberg(pool, problem, vector, frozen, fields_word.ITERATIONS)
    lag_report(problem, vector, heard)
    return problem, vector


# ---- variable projection: the events solved inside, the traverse outside ----

INNER_RUNGS_COLD = SECANTS                      # from the hint: every rung
INNER_RUNGS_WARM = (0.0001, SECANTS[-1])        # from a solved neighbour: the fine rungs
OUTER_RUNGS = (0.001, 0.0001)                   # s; the traverse's secants
OUTER_ITERATIONS = 6


def inner_solve(pool, problem: VoiceProblem, vector: numpy.ndarray, rungs: Sequence[float],
                iterations: int) -> numpy.ndarray:
    """Every word's events solved with the traverse held where the vector
    has it (the shared block frozen), over the given secant rungs."""
    frozen = numpy.zeros(problem.size(), dtype=bool)
    frozen[problem.offsets[-1]:] = True
    for secant in rungs:
        problem.secant = secant
        problem.steps = numpy.full(problem.size(), secant)
        vector, heard = levenberg(pool, problem, vector, frozen, iterations, quiet=True)
    return vector


def projected(pool, problem: VoiceProblem, vector: numpy.ndarray, traverse: numpy.ndarray,
              rungs: Sequence[float]) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """The vector with the given traverse and its events re-solved from the
    vector's, and the residual there: the projected (outer) function."""
    trial = vector.copy()
    trial[problem.offsets[-1]:] = traverse
    trial = inner_solve(pool, problem, trial, rungs, fields_word.ITERATIONS)
    return trial, problem.residual(trial)


def varpro(word_specs: Sequence[Sequence[str]], mode: str, cached: bool, groups: bool,
           label: str = "varpro") -> Tuple[VoiceProblem, numpy.ndarray]:
    """Golub-Pereyra variable projection over the shared traverse: the
    outer unknowns are the lags (and extras); each outer evaluation solves
    every word's events exactly (the event-only fit, which reaches -50 to
    -90 dB where the joint fit stalls at -24 to -36, steps 44, 52, 53). The
    outer Jacobian by central differences of the projected function, one
    inner solve a column side; damped Gauss-Newton with Nielsen's damping
    and the agreement gate, as in levenberg."""
    contexts = build_contexts(word_specs, mode)
    problem = VoiceProblem(contexts, SECANTS[0], mode, cached, groups)
    shared = slice(int(problem.offsets[-1]), problem.size())
    vector = problem.start_vector()
    truth = problem.true_vector()
    print("%s: %d words, %d shared unknowns over %d events; the truth in the family %.1f dB" % (
        label, len(problem.words), problem.size() - shared.start, shared.start,
        pyramid_word.level_db(problem.residual(truth))), flush=True)
    with multiprocessing.Pool(fields_word.POOL, initializer=worker_start,
                              initargs=(word_specs, SECANTS[0], mode, cached, groups)) as pool:
        vector = inner_solve(pool, problem, vector, INNER_RUNGS_COLD, fields_word.ITERATIONS)
        current = problem.residual(vector)
        cost = pyramid_word.level_db(current)
        report(problem, vector, cost, "%s: events solved at the start's traverse" % label)
        heard = numpy.zeros(problem.size(), dtype=bool)
        for secant in OUTER_RUNGS:
            print("%s outer secant %.3g s" % (label, secant), flush=True)
            for iteration in range(OUTER_ITERATIONS):
                traverse = vector[shared].copy()
                width = len(traverse)
                matrix = numpy.zeros((len(current), width))
                for column in range(width):
                    ahead = traverse.copy()
                    ahead[column] += secant
                    behind = traverse.copy()
                    behind[column] -= secant
                    _, ahead_residual = projected(pool, problem, vector, ahead, INNER_RUNGS_WARM)
                    _, behind_residual = projected(pool, problem, vector, behind, INNER_RUNGS_WARM)
                    matrix[:, column] = (ahead_residual - behind_residual) / (2.0 * secant)
                norms = numpy.sqrt(numpy.sum(matrix ** 2, axis=0))
                live = norms > 0.0
                heard[shared] = live
                if not numpy.any(live):
                    print("%s: no shared unknown is heard" % label, flush=True)
                    break
                scale = numpy.full(width, secant)
                scale /= numpy.sqrt(float(numpy.mean((norms[live] * scale[live]) ** 2)))
                scaled = matrix[:, live] * scale[live][None, :]
                left, singular, rows = numpy.linalg.svd(scaled, full_matrices=False)
                projected_residual = left.T @ (-current)
                damping = pyramid_word.DAMPING_START
                first = damping
                growth = 2.0
                accepted = False
                while not accepted and damping <= pyramid_word.DAMPING_LIMIT * first:
                    z = rows.T @ (singular / (singular ** 2 + damping) * projected_residual)
                    step = numpy.zeros(width)
                    step[live] = scale[live] * z
                    trial, trial_residual = projected(pool, problem, vector, traverse + step, INNER_RUNGS_WARM)
                    predicted = float(numpy.sum(current ** 2) - numpy.sum((current + matrix @ step) ** 2))
                    actual = float(numpy.sum(current ** 2) - numpy.sum(trial_residual ** 2))
                    ratio = actual / predicted if predicted > 0.0 else 0.0
                    if actual > 0.0 and ratio >= fields_word.AGREEMENT:
                        vector, current = trial, trial_residual
                        cost = pyramid_word.level_db(current)
                        damping *= max(1.0 / 3.0, 1.0 - (2.0 * ratio - 1.0) ** 3)
                        accepted = True
                    else:
                        damping *= growth
                        growth *= 2.0
                report(problem, vector, cost, "%s outer step %d%s" % (
                    label, iteration, "" if accepted else ": no damping lowered it"))
                lag_report(problem, vector, heard)
                if not accepted:
                    break
    return problem, vector


ALTERNATE_ROUNDS = 3


def alternate(word_specs: Sequence[Sequence[str]], mode: str, cached: bool, groups: bool,
              start_traverse: Optional[numpy.ndarray]) -> Tuple[VoiceProblem, numpy.ndarray]:
    """Block coordinate descent from a traverse: every word's events solved
    with the traverse held (the event-only fit, which reaches µs), then the
    traverse solved with the events held, in rounds. (The joint LM from six
    words left the traverse 1.5 ms off and the events up to 12 ms off,
    step 54.)"""
    contexts = build_contexts(word_specs, mode)
    problem = VoiceProblem(contexts, SECANTS[0], mode, cached, groups)
    shared = int(problem.offsets[-1])
    vector = problem.start_vector()
    if start_traverse is not None:
        vector[shared:] = start_traverse
    print("alternate: %d words, %d events, %d shared; the truth in the family %.1f dB" % (
        len(problem.words), shared, problem.size() - shared, pyramid_word.level_db(problem.residual(problem.true_vector()))),
        flush=True)
    events_frozen = numpy.zeros(problem.size(), dtype=bool)
    events_frozen[:shared] = True
    heard = numpy.zeros(problem.size(), dtype=bool)
    with multiprocessing.Pool(fields_word.POOL, initializer=worker_start,
                              initargs=(word_specs, SECANTS[0], mode, cached, groups)) as pool:
        for round_index in range(ALTERNATE_ROUNDS):
            rungs = INNER_RUNGS_COLD if round_index == 0 else INNER_RUNGS_WARM
            vector = inner_solve(pool, problem, vector, rungs, fields_word.ITERATIONS)
            report(problem, vector, pyramid_word.level_db(problem.residual(vector)),
                   "round %d, events solved" % round_index)
            for secant in OUTER_RUNGS:
                problem.secant = secant
                problem.steps = numpy.full(problem.size(), secant)
                vector, heard = levenberg(pool, problem, vector, events_frozen, fields_word.ITERATIONS, quiet=True)
            report(problem, vector, pyramid_word.level_db(problem.residual(vector)),
                   "round %d, traverse solved" % round_index)
            lag_report(problem, vector, heard)
    return problem, vector


def main() -> None:
    word_specs = [item.split(",") for item in sys.argv[1].split(";")]
    options = sys.argv[2:]
    mode = "traverse" if "traverse" in options else ("extras" if "extras" in options else "lags")
    cached = "plain" not in options   # the resumable residual: bit-identical on all six words (step 50)
    groups = "groups" in options
    if "alternate" in options:
        path = "%s/voice_traverse.npy" % sequence_recover.SCRATCH
        start = numpy.load(path) if os.path.exists(path) else None
        print("alternate from %s" % ("the saved traverse %s" % path if start is not None else "zero"), flush=True)
        problem, vector = alternate(word_specs, mode, cached, groups, start)
        numpy.save("%s/voice_alternate.npy" % sequence_recover.SCRATCH, vector)
        return
    if "varpro" in options:
        problem, vector = varpro(word_specs, mode, cached, groups)
        numpy.save("%s/voice_varpro_traverse.npy" % sequence_recover.SCRATCH, vector[problem.offsets[-1]:])
        return
    holdout = int(options[options.index("holdout") + 1]) if "holdout" in options else None
    if holdout is None:
        problem, vector = solve(word_specs, mode, cached=cached, groups=groups)
        numpy.save("%s/voice_traverse.npy" % sequence_recover.SCRATCH, vector[problem.offsets[-1]:])
        return
    # Cross-validation: the traverse from the other words, the held-out
    # word's events fitted under it, against none and the true one.
    training = [specs for index, specs in enumerate(word_specs) if index != holdout]
    trained, trained_vector = solve(training, mode, label="training", cached=cached, groups=groups)
    held = [word_specs[holdout]]
    held_problem = VoiceProblem(build_contexts(held, mode), SECANTS[0], mode, False, groups)
    # The learned traverse in the held-out word's shared block (a full-size
    # vector, the shared accessors slice it); an articulator the training
    # words never had stays at zero and is named.
    learned_full = numpy.zeros(held_problem.size())
    parts = [(trained.shared_lags, held_problem.shared_lags)]
    if second_block(mode):
        parts.append((trained.shared_extras, held_problem.shared_extras))
    unseen = []
    for trained_part, held_part in parts:
        for j, name in enumerate(held_problem.shared):
            if name in trained.shared:
                held_part(learned_full)[j] = trained_part(trained_vector)[trained.shared.index(name)]
            else:
                unseen.append(name)
    if unseen:
        print("articulators the training words lack (traverse zero): %s" % ", ".join(sorted(set(unseen))), flush=True)
    learned = learned_full[held_problem.offsets[-1]:]
    truth = held_problem.true_vector()[held_problem.offsets[-1]:]
    for label, traverse in (("held-out, no traverse", numpy.zeros(len(learned))), ("held-out, learned traverse", learned),
                            ("held-out, true traverse", truth)):
        solve(held, mode, fixed_traverse=traverse, label=label, cached=cached, groups=groups)


if __name__ == "__main__":
    main()
