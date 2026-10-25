"""Outside-in recovery of a velar release on a clock local in time: the
front's clock (the cavity in front of the seal, the seal node capped at the
seal) is the hinted clock plus a cubic B-spline correction, one knot a tick.
The window [closure, end] is fitted to zero flow through the sealed link
(pyramid_word.link_flow) with the coefficients reaching from the window's
start to one front's length past its end free (the front peel reads that far
ahead); the window grows a tick at a time, warm-started; the release is the
first end the fit cannot hold at the sealed floor. (On the global pyramid,
fitting the hold rewrote the motion, freezing it needed a carry-on, and the
carry-on and the peel's reach coupled into the window's last steps.)"""
import os
import sys
from typing import List, Tuple

sys.path.insert(0, os.path.expanduser("~/ducklips/tract_fast/analysis"))
import numpy

import fields
import fields_word as f
import pyramid_word as p
import sequence_recover as s

FIT_STEPS = 20
HOLD_START = os.environ.get("HOLDSTART", "0") == "1"
HISTORY = 30       # ticks before the first window whose clock its fit also frees (the peel's memory, ~25 ms)
RECENT = int(os.environ.get("RECENT", "10"))   # ticks before a window's end whose coefficients a growing step frees
FIT_RATIO = 100.0   # a window is sealed if its fitted flow is within this of the established floor (insensitive 10-1e4)


def cubic_bspline_basis(ticks: int) -> numpy.ndarray:
    """Uniform cubic B-splines, one knot a tick (columns centred on each
    tick, support four ticks)."""
    basis = numpy.zeros((ticks, ticks))
    positions = numpy.arange(ticks, dtype=numpy.float64)
    for centre in range(ticks):
        distance = numpy.abs(positions - centre)
        inner = distance < 1.0
        outer = (distance >= 1.0) & (distance < 2.0)
        basis[inner, centre] = 2.0 / 3.0 - distance[inner] ** 2 + 0.5 * distance[inner] ** 3
        basis[outer, centre] = (2.0 - distance[outer]) ** 3 / 6.0
    return basis


class Front:
    def __init__(self, specs: List[str]) -> None:
        self.context = p.Context(specs, ["areas"])
        self.word = self.context.word
        sealing = [k for k in self.word.stops if k not in self.word.narrowings]
        # The stop to recover (STOP: its phone index; else the first that
        # seals and is not word-initial, else the first).
        chosen = int(os.environ["STOP"]) if os.environ.get("STOP") else None
        later = [k for k in sealing if k > 0]
        self.stop = chosen if chosen is not None else (later[0] if later else sealing[0])
        self.node = self.word.stops[self.stop]
        self.problem = f.TimeProblem(self.context, self.context.levels)
        self.hinted = fields.progress(self.word.hinted_arrivals, self.word.hinted_leavings, self.context.ticks, 0.001)
        # The front's base clock holds the stop's target from its hinted
        # arrival on (no motion of its own): the correction builds the move
        # toward the next phone when the sound calls for it. (With the hint's
        # own early move in the base, the fit had to cancel it and the added
        # ticks crept from -348 to -108 dB, read as Delta's release 35 ms
        # early.)
        held_from = int(self.word.hinted_arrivals[self.stop] * 1000.0)
        # The hint's own hold (before the base is held), for the anchors below.
        on_target = numpy.abs(self.hinted - numpy.round(self.hinted)) < 1e-12
        self.hinted = self.hinted.copy()
        self.hinted[held_from:] = numpy.minimum(self.hinted[held_from:], float(self.stop))
        self.basis = cubic_bspline_basis(self.context.ticks)
        self.reach = self.context.count - self.node
        # Anchored: the coefficients whose support lies wholly inside the
        # hinted hold (the hinted clock on its target there, and the true hold
        # contains the hinted one) stay zero, so the front holds its target.
        # (Free, the fit settled on a sound-alike front 7% toward the
        # fricative's shape, 1.932 against 2.000, and later ran backward; a
        # running maximum's kinks broke the fit.) From the hint's hold, not the
        # held base: the base stays on the stop to the word's end, and anchors
        # read from it froze every coefficient after the arrival (Delta's fit
        # then took no step as the true clock left the stop).
        self.anchored = numpy.array([bool(numpy.all(on_target[max(c - 2, 0):min(c + 3, self.context.ticks)]))
                                     for c in range(self.context.ticks)])

    def command(self, coefficients: numpy.ndarray, sealed_from: int = 0) -> numpy.ndarray:
        """The seal node capped at the seal only from `sealed_from` (a step)
        on: capped over the whole word, it changed the link before the
        closure and, under a fricative's jet, became the narrowest node there
        (Oscar's zero flow reached only -64.6 dB with the true clock)."""
        context = self.context
        clock = self.hinted + self.basis @ coefficients
        time = fields.same_time(self.hinted, context.count)
        time[:, self.node:] = clock[:, None]
        command = fields.read_postures(self.problem.rows, context.count, self.problem.space, time)
        tick = sealed_from // p.TICK
        command[tick:, self.node] = numpy.minimum(command[tick:, self.node], self.word.seal)
        return command

    def residual(self, coefficients: numpy.ndarray, window: Tuple[int, int]) -> numpy.ndarray:
        flow = p.link_flow(self.context, self.command(coefficients, window[0]), self.node)
        part = flow[window[0]:window[1]] / self.context.scale
        return part if numpy.all(numpy.isfinite(part)) else numpy.full(len(part), 1e3)

    def free_start(self, window: Tuple[int, int]) -> numpy.ndarray:
        """Growing backwards: the coefficients near the window's start (and
        one front's length past it: the peel at the start reads that far
        ahead)."""
        first = max(window[0] // p.TICK - 3, 0)
        last = min((window[0] + self.reach) // p.TICK + RECENT, self.context.ticks)
        return numpy.arange(first, last)

    def free(self, window: Tuple[int, int]) -> numpy.ndarray:
        """The coefficients near the window's end and one front's length past
        it: with a local basis the newly added stretch moves only those (the
        rest were fitted with the earlier windows). The whole window when it
        is short (the established one)."""
        first = max(window[0] // p.TICK - 2, window[1] // p.TICK - RECENT, 0)
        last = min((window[1] + self.reach) // p.TICK + 3, self.context.ticks)
        return numpy.arange(first, last)

    def fit(self, coefficients: numpy.ndarray, window: Tuple[int, int],
            backwards: bool = False, history: int = 0) -> Tuple[numpy.ndarray, float]:
        step_size = numpy.cbrt(numpy.finfo(numpy.float64).eps)
        current = coefficients.copy()
        free = self.free_start(window) if backwards else self.free(window)
        if history > 0:
            first = max(window[0] // p.TICK - history, 0)
            last = min((window[1] + self.reach) // p.TICK + 3, self.context.ticks)
            free = numpy.arange(first, last)
        free = free[~self.anchored[free]]
        if len(free) == 0:
            return current, p.level_db(self.residual(current, window))
        for iteration in range(FIT_STEPS):
            base = self.residual(current, window)
            matrix = numpy.empty((len(base), len(free)))
            for column, index in enumerate(free):
                ahead = current.copy()
                ahead[index] += step_size
                behind = current.copy()
                behind[index] -= step_size
                matrix[:, column] = (self.residual(ahead, window) - self.residual(behind, window)) / (2 * step_size)
            left, singular, rows = numpy.linalg.svd(matrix, full_matrices=False)
            if singular[0] == 0.0:
                # A flat window (deep in a seal every column is exactly zero):
                # no direction is heard, so there is no step to take.
                break
            projected = left.T @ (-base)
            damping = 1e-3 * float(singular[0] ** 2)
            accepted = False
            for _ in range(30):
                move = rows.T @ (singular / (singular ** 2 + damping) * projected)
                trial = current.copy()
                trial[free] += move
                if numpy.sum(self.residual(trial, window) ** 2) < numpy.sum(base ** 2):
                    current = trial
                    accepted = True
                    break
                damping *= 4.0
            if not accepted:
                break
        return current, p.level_db(self.residual(current, window))


def main():
    front = Front(sys.argv[1].split(","))
    context = front.context
    word = front.word
    hinted_command = front.problem.command(front.problem.vector_of(front.hinted))
    closure, _ = p.read_steps(context, hinted_command, word.hinted_arrivals, word.hinted_leavings)[front.stop]
    true_span = p.sealed_span(context, context.true_command, front.node, closure + 20 * p.TICK)
    begin = closure + 2 * p.TICK
    coefficients = numpy.zeros(context.ticks)
    # Check: the true clock through the basis.
    truth = fields.progress(word.true_arrivals, word.true_leavings, context.ticks, 0.001)
    true_coefficients = numpy.linalg.solve(front.basis, truth - front.hinted)
    check_end = true_span[1] - 2 * p.TICK if true_span[1] > begin else begin + 40 * p.TICK
    print("closure read %d; true sealed span %s; true clock through the basis: window to %d %.1f dB" % (
        closure, true_span, check_end, p.level_db(front.residual(true_coefficients, (begin, check_end)))), flush=True)
    # The first window in the middle of the stop's hinted hold, where the
    # front is surely still (from the read closure, Oscar's first window fell
    # in the /s/-to-/k/ transition and fitted only to -34.8 dB); grown
    # forward to the release and backward to the closure.
    if HOLD_START:
        # The first window in the middle of the hinted hold, where the front is
        # anchored on its target (the true hold contains the hinted one), the
        # past before it free: from there backward to the closure (the
        # transition fitted against the right boundary) and forward to the
        # release. (From the read closure, Oscar's transition settled on a
        # sound-alike front 6% short of its target.)
        middle = 0.5 * (word.hinted_arrivals[front.stop] + word.hinted_leavings[front.stop])
        begin = int((middle - 0.01) * s.STEP_RATE)
    end = begin + 20 * p.TICK
    # The first fit also frees the clock over the stretch before the window
    # (back to HISTORY ticks before it): the front peel's state entering the
    # window carries the front's past (Oscar's /s/ with its jet, at the
    # hinted timing), a transient lasting ~25 ms, so it is fitted, though
    # scored only inside the window. (Anchored in the hold with the past at
    # the hint, the first window fitted only to -56.4 dB.)
    # From a hold start the past to free reaches back past the read closure
    # (the front cavity, closed at the seal end, keeps the fricative's era in
    # the peel's state for longer than 30 ticks: anchored in the hold with 30
    # ticks free, -54.9 dB).
    history = HISTORY if not HOLD_START else (begin - (closure - HISTORY * p.TICK)) // p.TICK
    coefficients, floor = front.fit(coefficients, (begin, end), history=history)
    print("established window to %d: link flow %.1f dB" % (end, floor), flush=True)
    limit = int(word.hinted_arrivals[front.stop + 1] * s.STEP_RATE) + 40 * p.TICK
    # The sealed/open boundary: halfway, in dB, between the sealed floor and
    # the level of a window that surely reaches open link (to the next
    # phone's hinted arrival and 40 ticks on), measured unfitted. (A fixed
    # 100x above the floor misfired on Delta's -222 dB floor: a small
    # imperfection at -169 dB read as open, where open link is ~-30.)
    open_level = p.level_db(front.residual(coefficients, (begin, limit)))
    # The newly added tick is tested against this midpoint (below). (A jump
    # of the whole window by half the gap failed both ends on Delta: closure
    # 0, release 59336.)
    middle_level = 0.5 * (floor + open_level)
    print("sealed floor %.1f dB, open level %.1f dB: a new stretch above %.1f dB is open" % (
        floor, open_level, middle_level), flush=True)
    last = floor
    while end + p.TICK < limit:
        trial, level = front.fit(coefficients, (begin, end + p.TICK))
        # The newly added tick alone, after the fit: a sealed tick fits near
        # the floor, an open one cannot be pushed there whatever the fit does
        # elsewhere. (The whole window's level crept with the fit's
        # imperfection, -229.7 -> -111.7 dB, and absorbed a single open
        # tick, so neither a midpoint nor a jump of the whole window worked.)
        added = p.level_db(front.residual(trial, (end, end + p.TICK)))
        sealed = added <= middle_level
        clock = front.hinted + front.basis @ trial
        tick = (end + p.TICK) // p.TICK
        print("  window to %d: link flow %.1f dB, the added tick %.1f dB -> %s; front clock there %.4f (truth %.4f)" % (
            end + p.TICK, level, added, "sealed" if sealed else "open", clock[tick - 1], truth[tick - 1]), flush=True)
        if not sealed:
            break
        end, coefficients, last = end + p.TICK, trial, level
    release_level = last
    print("release found at step %d (true %d)" % (end, true_span[1]), flush=True)
    # The closure, backwards: the window's start moved earlier a tick at a
    # time while zero flow still fits (before the closure the link is open).
    while begin - p.TICK > 0:
        trial, level = front.fit(coefficients, (begin - p.TICK, end), backwards=True)
        added = p.level_db(front.residual(trial, (begin - p.TICK, end))[:p.TICK])
        sealed = added <= middle_level
        print("  window from %d: link flow %.1f dB -> %s" % (begin - p.TICK, level, "sealed" if sealed else "open"),
              flush=True)
        if not sealed:
            break
        begin, coefficients, last = begin - p.TICK, trial, level
    print("closure found at step %d (true %d)" % (begin, true_span[0]), flush=True)
    # To the step: bisection between the last sealed and the first open
    # window, at each end.
    sealed_at, open_at = begin, begin - p.TICK
    while sealed_at - open_at > 1:
        middle = (sealed_at + open_at) // 2
        trial, level = front.fit(coefficients, (middle, end), backwards=True)
        if p.level_db(front.residual(trial, (middle, end))[:max(sealed_at - middle, 1)]) <= middle_level:
            sealed_at, coefficients = middle, trial
        else:
            open_at = middle
    begin = sealed_at
    sealed_to, open_to = end, end + p.TICK
    while open_to - sealed_to > 1:
        middle = (sealed_to + open_to) // 2
        trial, level = front.fit(coefficients, (begin, middle))
        if p.level_db(front.residual(trial, (sealed_to, middle))) <= middle_level:
            sealed_to, coefficients = middle, trial
        else:
            open_to = middle
    end = sealed_to
    print("to the step: closure %d (true %d), release %d (true %d)" % (begin, true_span[0], end, true_span[1]),
          flush=True)


if __name__ == "__main__":
    main()
