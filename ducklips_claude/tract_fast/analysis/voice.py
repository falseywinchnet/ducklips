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

import concurrent.futures
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
FAST = os.environ.get("VOICE_SLOW", "0") != "1"   # the fast path (step 59); VOICE_SLOW=1 the step-58 solver


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
        self.workers = fields_word.POOL   # the pool's size, set by the solver that makes the pool
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
        # Bounded: a subclass may append its own unknowns after the shared block.
        return vector[self.offsets[-1] + len(self.shared):self.offsets[-1] + 2 * len(self.shared)]

    def residual(self, vector: numpy.ndarray) -> numpy.ndarray:
        parts = []
        for index, word in enumerate(self.words):
            parts.append(word.residual(self.word_vector(vector, index)))
        return numpy.concatenate(parts)

    def columns_of(self, index: int, active: numpy.ndarray) -> numpy.ndarray:
        """The columns word `index`'s rows depend on: its own events and the
        shared quantities of the articulators it has, those not frozen."""
        own = numpy.arange(int(self.offsets[index]), int(self.offsets[index + 1]))
        used = numpy.unique(self.maps[index])
        shared = int(self.offsets[-1]) + used
        if second_block(self.mode):
            shared = numpy.concatenate([shared, shared + len(self.shared)])
        columns = numpy.concatenate([own, shared])
        return columns[active[columns]]

    def steps_for(self, secant: float) -> numpy.ndarray:
        """Every unknown's difference step at a rung: all times here."""
        return numpy.full(self.size(), secant)

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


def report(problem: VoiceProblem, vector: numpy.ndarray, cost: float, label: str,
           current: Optional[numpy.ndarray] = None) -> None:
    """The step's line: all words' level, each word's level and events off
    the truth, the traverse off the truth. With `current` (the vector's
    residual) the word levels are read from it, not re-evaluated (nine
    sequential residuals a step in the main process, ~2 s, step 59)."""
    truth = problem.true_vector()
    blocks = problem.blocks()
    parts = []
    for index, word in enumerate(problem.words):
        own_truth = problem.word_vector(truth, index).copy()
        # word_vector sets the word's traverse from the vector it reads: the
        # vector's own last, so the block's level is the vector's.
        own = problem.word_vector(vector, index)
        if current is not None:
            block_cost = pyramid_word.level_db(current[blocks[index]])
        else:
            block_cost = pyramid_word.level_db(word.residual(own))
        parts.append("%s %.1f dB, events %.3f ms max off" % (
            "/".join(word.context.word.names), block_cost, 1000.0 * numpy.abs(own - own_truth).max()))
    shared_end = int(problem.offsets[-1]) + (2 if second_block(problem.mode) else 1) * len(problem.shared)
    traverse = vector[problem.offsets[-1]:shared_end]
    true_traverse = truth[problem.offsets[-1]:shared_end]
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


# ---- the fast path: one task a word, Broyden updates between Jacobians (step 59) ----

def word_block(task: Tuple[numpy.ndarray, int, numpy.ndarray, numpy.ndarray]) -> Tuple[int, numpy.ndarray, numpy.ndarray]:
    """One word's rows of the Jacobian for the task's columns, central
    differences over the task's steps, every residual resumed in this
    process's cache for the word (warmed at the base first). Returns the
    word, its columns, the block (rows x columns)."""
    vector, index, steps, columns = task
    word = _voice.words[index]
    block = numpy.zeros((word.context.word.size, len(columns)))
    word.residual(_voice.word_vector(vector, index))
    for position, column in enumerate(columns):
        step = steps[column]
        ahead = vector.copy()
        ahead[column] += step
        behind = vector.copy()
        behind[column] -= step
        block[:, position] = (word.residual(_voice.word_vector(ahead, index))
                              - word.residual(_voice.word_vector(behind, index))) / (2.0 * step)
    return index, columns, block


def word_residual(task: Tuple[numpy.ndarray, int]) -> Tuple[int, numpy.ndarray]:
    vector, index = task
    return index, _voice.words[index].residual(_voice.word_vector(vector, index))


def task_order(problem: VoiceProblem) -> List[int]:
    """The words longest first (a pool finishes when its longest task does)."""
    sizes = numpy.array([word.context.word.size for word in problem.words])
    return [int(index) for index in numpy.argsort(-sizes, kind="stable")]


class PinnedPool:
    """Single-process executors, one a slot, and a fixed plan of which slot
    computes which word's columns: a worker then caches only its own words.
    (With any task to any worker every worker cached all nine words, 2.1 GB
    each warm; a 4-worker pod with 8 GB lost a worker to the memory limit,
    and 6 workers on the Mini's 16 GB went to swap, step 63.) The plan is
    McNaughton's wrap-around rule: the slots filled in turn to an equal share
    of the Jacobian's work (columns times residual length), a word split
    where a slot fills, so the wall time stays the balanced one."""

    def __init__(self, slots: int, initializer, initargs: Tuple) -> None:
        self.slots = slots
        self.executors = [concurrent.futures.ProcessPoolExecutor(1, initializer=initializer, initargs=initargs)
                          for _ in range(slots)]
        self.plans: Dict[bytes, Tuple[List[List[Tuple[int, numpy.ndarray]]], Dict[int, int]]] = {}

    def __enter__(self) -> "PinnedPool":
        return self

    def __exit__(self, kind, value, trace) -> None:
        for executor in self.executors:
            executor.shutdown(wait=True, cancel_futures=True)

    def plan(self, problem: "VoiceProblem", frozen: numpy.ndarray) -> Tuple[List[List[Tuple[int, numpy.ndarray]]], Dict[int, int]]:
        """Per slot its (word, columns) shares, and per word the slot that
        evaluates its residuals (the one holding most of its columns)."""
        key = frozen.tobytes()
        found = self.plans.get(key)
        if found is not None:
            return found
        active = ~frozen
        items = [(index, problem.columns_of(index, active)) for index in task_order(problem)]
        sizes = {index: problem.words[index].context.word.size for index, _ in items}
        total = float(sum(len(columns) * sizes[index] for index, columns in items))
        target = total / self.slots
        shares: List[List[Tuple[int, numpy.ndarray]]] = [[] for _ in range(self.slots)]
        held: Dict[int, Dict[int, int]] = {}
        slot = 0
        filled = 0.0
        for index, columns in items:
            size = sizes[index]
            remaining = columns
            while len(remaining) > 0:
                if slot == self.slots - 1:
                    take = len(remaining)
                else:
                    take = int(numpy.clip(round((target - filled) / size), 1, len(remaining)))
                shares[slot].append((index, remaining[:take]))
                counts = held.setdefault(index, {})
                counts[slot] = counts.get(slot, 0) + take
                remaining = remaining[take:]
                filled += take * size
                if slot < self.slots - 1 and filled >= target - 0.5 * size:
                    slot += 1
                    filled = 0.0
        owner = {}
        for index, _ in items:
            counts = held.get(index, {0: 0})
            owner[index] = max(counts, key=counts.__getitem__)
        found = (shares, owner)
        self.plans[key] = found
        return found


def run_tasks(pool, function, tasks: List) -> List:
    """`function`'s results over `tasks`, in completion order, from a
    multiprocessing.Pool or a concurrent.futures executor (whose dead
    worker raises BrokenProcessPool, where the Pool waited forever)."""
    if hasattr(pool, "imap_unordered"):
        return list(pool.imap_unordered(function, tasks))
    futures = [pool.submit(function, task) for task in tasks]
    return [future.result() for future in concurrent.futures.as_completed(futures)]


def residual_parallel(pool, problem: VoiceProblem, vector: numpy.ndarray) -> numpy.ndarray:
    """The voice's residual, one task a word (on a PinnedPool, in the slot
    that holds the word's cache)."""
    blocks = problem.blocks()
    current = numpy.empty(blocks[-1].stop)
    tasks = [(vector, i) for i in task_order(problem)]
    if isinstance(pool, PinnedPool):
        owner = pool.plan(problem, numpy.zeros(problem.size(), dtype=bool))[1]
        futures = [pool.executors[owner[task[1]]].submit(word_residual, task) for task in tasks]
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    else:
        results = run_tasks(pool, word_residual, tasks)
    for index, part in results:
        current[blocks[index]] = part
    return current


TASKS_A_WORKER = 3   # tasks a worker on average: enough to even out cores of two speeds (4 performance, 6 efficiency)


def jacobian_tasks(problem: VoiceProblem, vector: numpy.ndarray, frozen: numpy.ndarray, workers: int) -> List:
    """A Jacobian's tasks: each word's columns (its own events and the
    shared columns it uses) cut into chunks, the longest words' first, so
    that the slowest task no longer sets the wall time (one task a word,
    Delta's 22 columns on an efficiency core did, step 62). Each chunk costs
    one warm-up residual more."""
    active = ~frozen
    per_word = [(index, problem.columns_of(index, active)) for index in task_order(problem)]
    total = sum(len(columns) for _, columns in per_word)
    chunk = max(1, int(numpy.ceil(total / float(TASKS_A_WORKER * workers))))
    tasks = []
    for index, columns in per_word:
        for first in range(0, len(columns), chunk):
            tasks.append((vector, index, problem.steps, columns[first:first + chunk]))
    return tasks


def jacobian_by_word(pool, problem: VoiceProblem, vector: numpy.ndarray, frozen: numpy.ndarray) -> numpy.ndarray:
    """The voice's Jacobian in word-and-column chunks: each word's rows for
    its own events and the shared columns it uses (zero elsewhere). (Spread
    by column instead, every worker evaluated all nine words for each shared
    column, 82% of a step's residuals.)"""
    blocks = problem.blocks()
    matrix = numpy.zeros((blocks[-1].stop, problem.size()))
    if isinstance(pool, PinnedPool):
        shares = pool.plan(problem, frozen)[0]
        futures = []
        for slot, share in enumerate(shares):
            for index, columns in share:
                futures.append(pool.executors[slot].submit(word_block, (vector, index, problem.steps, columns)))
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    else:
        results = run_tasks(pool, word_block, jacobian_tasks(problem, vector, frozen, problem.workers))
    for index, columns, block in results:
        matrix[blocks[index], columns] = block
    return matrix


def levenberg_broyden(pool, problem: VoiceProblem, vector: numpy.ndarray, frozen: numpy.ndarray,
                      iterations: int, quiet: bool = False) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """levenberg's step (damped Gauss-Newton in the secant's units, Nielsen's
    damping, the agreement gate), with the Jacobian computed by word and,
    after an accepted step, corrected by Broyden's rank-one update with the
    residual's change the step itself measured: J += (dr - J s) s' / s's.
    A full Jacobian again only when the gate rejects a step on an updated
    one. Returns the vector and which unknowns the sound heard."""
    current = residual_parallel(pool, problem, vector)
    cost = pyramid_word.level_db(current)
    matrix = jacobian_by_word(pool, problem, vector, frozen)
    fresh = True
    jacobians = 1
    residuals = 1
    heard = numpy.zeros(len(vector), dtype=bool)
    accepted_steps = 0
    damping = pyramid_word.DAMPING_START
    while accepted_steps < iterations:
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
        first = pyramid_word.DAMPING_START
        growth = 2.0
        accepted = False
        refreshed = False
        while not accepted and damping <= pyramid_word.DAMPING_LIMIT * first:
            z = rows.T @ (singular / (singular ** 2 + damping) * projected)
            step = numpy.zeros(len(vector))
            step[heard] = scale[heard] * z
            trial = residual_parallel(pool, problem, vector + step)
            residuals += 1
            predicted = float(numpy.sum(current ** 2) - numpy.sum((current + matrix @ step) ** 2))
            actual = float(numpy.sum(current ** 2) - numpy.sum(trial ** 2))
            ratio = actual / predicted if predicted > 0.0 else 0.0
            if actual > 0.0 and ratio >= fields_word.AGREEMENT:
                if numpy.all(numpy.isfinite(trial)):
                    change = trial - current
                    matrix += numpy.outer(change - matrix @ step, step) / float(step @ step)
                    fresh = False
                vector, current = vector + step, trial
                cost = pyramid_word.level_db(current)
                damping *= max(1.0 / 3.0, 1.0 - (2.0 * ratio - 1.0) ** 3)
                accepted = True
            elif not fresh:
                # The updated Jacobian disagreed with the sound: a full one,
                # the same damping, the step taken again.
                matrix = jacobian_by_word(pool, problem, vector, frozen)
                fresh = True
                jacobians += 1
                refreshed = True
                break
            else:
                damping *= growth
                growth *= 2.0
        if refreshed:
            continue
        if accepted:
            accepted_steps += 1
        if not quiet:
            report(problem, vector, cost, "step %d (%d Jacobians, %d residuals of all words)%s" % (
                accepted_steps, jacobians, residuals, "" if accepted else ": no damping lowered it"), current)
        if not accepted:
            break
        damping = max(damping, pyramid_word.DAMPING_START)
    return vector, heard


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
    if FAST:
        # One pool for every rung (the secant travels with each task): the
        # workers build the words once, not once a rung.
        workers = fields_word.POOL
        problem.workers = workers
        with concurrent.futures.ProcessPoolExecutor(workers, initializer=worker_start,
                                                    initargs=(word_specs, SECANTS[0], mode, cached, groups)) as pool:
            for secant in SECANTS:
                problem.secant = secant
                problem.steps = problem.steps_for(secant)
                print("%s secant %.3g s" % (label, secant), flush=True)
                vector, heard = levenberg_broyden(pool, problem, vector, frozen, fields_word.ITERATIONS)
        lag_report(problem, vector, heard)
        return problem, vector
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
