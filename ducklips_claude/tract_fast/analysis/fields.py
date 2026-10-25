"""The motion as postures read through two coordinate fields:

    A(x, t) = P(xi(x, t), s(x, t))

P holds the word's postures: each phone's tube along the tube's node
coordinate u (read between nodes linearly), with s = k at phone k's target,
linear in s between targets. xi is the space field (where along the
postures node x sits at tick t: a constriction sliding, the place moving);
s is the time field (where in the word node x is at tick t: a transition
early or late, the lips leading the tongue). Both are smooth and carried by
pyramids (pyramid.py); timing and place are then coarse coordinates, not
fine area detail. The port and the folds have no place: they are read
through the time field alone (its first node's row).

A seal at node x at a moment read off the sound is a level set of the
fields there: s(x, t_event) at the value where P's area crosses the seal.
"""
from typing import Optional, Tuple

import numpy


def progress(arrivals: numpy.ndarray, leavings: numpy.ndarray, ticks: int, tick_seconds: float) -> numpy.ndarray:
    """A time field (one value per tick) for phone times: k held from
    arrival to leaving, each transition k + (1 - cos(pi f)) / 2 for f the
    share of the transition (smooth: its slope zero at both ends)."""
    field = numpy.zeros(ticks)
    for tick in range(ticks):
        now = tick * tick_seconds
        value = 0.0
        for k in range(len(arrivals)):
            if now >= arrivals[k]:
                value = float(k)
            if k + 1 < len(arrivals) and leavings[k] < now < arrivals[k + 1]:
                share = (now - leavings[k]) / (arrivals[k + 1] - leavings[k])
                value = k + 0.5 * (1.0 - numpy.cos(numpy.pi * share))
                break
        field[tick] = value
    return field


def read_postures(rows: numpy.ndarray, nodes: int, space: numpy.ndarray, time: numpy.ndarray,
                  side_time: Optional[numpy.ndarray] = None) -> numpy.ndarray:
    """The command (ticks x columns) from postures `rows` (phones x columns;
    the first `nodes` columns areas along the tube, the rest without place)
    read through the space field (ticks x nodes: the posture coordinate of
    each node) and the time field (ticks x nodes; the columns past the
    areas, the port and the folds, read through `side_time` (one value a
    tick: their own clock) or else node 0's)."""
    phones = rows.shape[0]
    ticks = time.shape[0]
    # The time field between targets: phone k and the share toward k + 1.
    clamped = numpy.clip(time, 0.0, phones - 1.0)
    lower = numpy.minimum(numpy.floor(clamped).astype(numpy.int64), phones - 2) if phones > 1 else \
        numpy.zeros_like(clamped, dtype=numpy.int64)
    share = clamped - lower if phones > 1 else numpy.zeros_like(clamped)
    # The space field between nodes: node j and the share toward j + 1.
    place = numpy.clip(space, 0.0, nodes - 1.0)
    left = numpy.minimum(numpy.floor(place).astype(numpy.int64), nodes - 2)
    across = place - left
    areas = rows[:, :nodes]
    upper = numpy.minimum(lower + 1, phones - 1)
    at_lower = (1.0 - across) * areas[lower, left] + across * areas[lower, left + 1]
    at_upper = (1.0 - across) * areas[upper, left] + across * areas[upper, left + 1]
    command = numpy.empty((ticks, rows.shape[1]))
    command[:, :nodes] = (1.0 - share) * at_lower + share * at_upper
    if rows.shape[1] > nodes:
        if side_time is None:
            first_lower = lower[:, 0]
            first_upper = upper[:, 0]
            first_share = share[:, 0:1]
        else:
            side = numpy.clip(side_time, 0.0, phones - 1.0)
            first_lower = numpy.minimum(numpy.floor(side).astype(numpy.int64), max(phones - 2, 0))
            first_upper = numpy.minimum(first_lower + 1, phones - 1)
            first_share = (side - first_lower)[:, None]
        command[:, nodes:] = (1.0 - first_share) * rows[first_lower, nodes:] + first_share * rows[first_upper, nodes:]
    return command


def identity_space(ticks: int, nodes: int) -> numpy.ndarray:
    """The space field that reads each node at its own place."""
    return numpy.tile(numpy.arange(nodes, dtype=numpy.float64), (ticks, 1))


def same_time(time: numpy.ndarray, nodes: int) -> numpy.ndarray:
    """A time field shared by every node (ticks x nodes)."""
    return numpy.tile(time[:, None], (1, nodes))


def level_crossing(rows: numpy.ndarray, node: int, phone: int, level: float) -> Tuple[float, float]:
    """Where between phone - 1, phone and phone + 1 the area of `node` (read
    at its own place) crosses `level`: the time field's values (closure,
    release) at the seal, P linear in s between targets. NaN where it does
    not cross."""
    found = []
    for first, second in ((phone - 1, phone), (phone, phone + 1)):
        if first < 0 or second >= rows.shape[0]:
            found.append(float("nan"))
            continue
        a, b = rows[first, node], rows[second, node]
        if (a - level) * (b - level) > 0.0 or a == b:
            found.append(float("nan"))
            continue
        found.append(first + (level - a) / (b - a))
    return found[0], found[1]
