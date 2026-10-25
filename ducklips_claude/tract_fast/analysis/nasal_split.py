"""The sound split between lips and nostrils, for the reversible core.

One microphone hears both openings. The nose's share s is the unknown: the
lips then carry y - s. The oral peel (from the lips) gives the pressure at
the vocal node the nose opens from; the nasal peel (from the nostrils) gives
the nose's first node and the wave the port delivered into it, and through
the port's series element (Tube's, one hop) the wave the vocal node must have
sent and the wave the port returned. Step by step the vocal node's pressure
must equal the two waves' sum there: one condition per step, met only by the
true split. Solved by Newton's method with GMRES on products of two peels
each; then the oral peel continues below the port with the returned wave
known, to the glottis."""
import os
import sys

import numpy
from scipy.sparse.linalg import LinearOperator, gmres

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reverse_core

RHO_C = reverse_core.RHO_C
NODE = reverse_core.NODE_LENGTH


class Branched:
    """A moving vocal tube with the nose: its geometry, the port's area per
    step, the nose (reverse_core.read_nose), the hop low-pass and the open
    end's extra length."""

    def __init__(self, geometry, port, nose, hop, end_extra):
        self.geometry = numpy.ascontiguousarray(geometry, float)
        self.port = numpy.asarray(port, float)
        self.port_node, self.nose_areas, self.nose_targets, self.nose_walls, self.nose_hop, scale, first = nose
        self.nose_walls = numpy.ascontiguousarray(self.nose_walls, float)
        self.glide_start = int(round(16.0 * float(numpy.mean(first / self.nose_targets)))) - 1
        self.hop = hop
        self.end_extra = end_extra
        self.port_y = self.port / RHO_C


def port_waves(model, incoming, nose_pressure):
    """From the nasal peel: per step the wave the vocal node sent toward the
    port the step before, and the wave the port returned into it."""
    count = len(incoming)
    leaving = nose_pressure - incoming
    sent_before = numpy.zeros(count)
    returned = numpy.zeros(count)
    for t in range(1, count):
        area = model.port[t]
        if area < 1e-6:
            continue
        z = RHO_C / area
        resistance = max(0.03 * max(area, 1e-4) ** -2.6, 0.0065 / max(area, 1e-4) * numpy.sqrt(2250.0)) * NODE
        flow = (incoming[t] - leaving[t - 1]) / z
        drive = (2.0 * z + resistance) * flow + 2.0 * reverse_core._kinetic_term(
            flow, model.geometry[t - 1, model.port_node, 0], model.nose_areas[0])
        sent_before[t] = leaving[t - 1] + 0.5 * drive
        returned[t] = sent_before[t] - z * flow
    return sent_before, returned


def mismatch(model, sound, share):
    """The junction's condition per step for the nose's share (zero for the
    true split)."""
    _, oral = reverse_core._peel_moving_side(model.geometry, sound - share, model.hop, True, model.end_extra, 0.0,
                                             reverse_core.JET_SEED, model.port_node, numpy.zeros(1), numpy.zeros(1),
                                             model.port_node)
    incoming, nose_pressure = reverse_core._peel_nose(model.nose_areas, model.nose_walls, model.nose_targets,
                                                      model.glide_start, model.nose_hop, True, share, model.port_y)
    sent_before, returned = port_waves(model, incoming, nose_pressure)
    condition = numpy.zeros(len(sound))
    for t in range(1, len(sound)):
        if model.port[t] >= 1e-6:
            condition[t - 1] = oral[t - 1] - returned[t - 1] - sent_before[t]
        else:
            condition[t - 1] = incoming[t] * RHO_C
    return condition


def split(model, sound, start=None, rounds=4, iterations=200, quiet=False):
    """The nose's share of the sound (Newton's method, GMRES inside)."""
    share = numpy.zeros(len(sound)) if start is None else numpy.array(start, float)
    for round_index in range(rounds):
        condition = mismatch(model, sound, share)
        size = float(numpy.sqrt(numpy.mean(condition ** 2)))
        if not quiet:
            print("  split round %d: condition rms %.3e" % (round_index, size))
            sys.stdout.flush()
        if size < 1e-9:
            break
        step = 1e-3 * max(float(numpy.sqrt(numpy.mean(share ** 2))), 1.0)
        operator = LinearOperator((len(sound), len(sound)), matvec=Product(model, sound, share, condition, step).apply,
                                  dtype=float)
        change, info = gmres(operator, -condition, restart=80, maxiter=iterations, rtol=1e-10)
        share = share + change
    return share


class Product:
    """The condition's derivative times a direction, by a difference."""

    def __init__(self, model, sound, share, condition, step):
        self.model = model
        self.sound = sound
        self.share = share
        self.condition = condition
        self.step = step

    def apply(self, direction):
        size = float(numpy.sqrt(numpy.mean(direction ** 2)))
        if size == 0.0:
            return numpy.zeros_like(direction)
        scale = self.step / size
        return (mismatch(self.model, self.sound, self.share + scale * direction) - self.condition) / scale


def glottal(model, sound, share):
    """The flow entering the vocal tube, the oral peel continued below the
    port with the wave it returned."""
    incoming, nose_pressure = reverse_core._peel_nose(model.nose_areas, model.nose_walls, model.nose_targets,
                                                      model.glide_start, model.nose_hop, True, share, model.port_y)
    sent_before, returned = port_waves(model, incoming, nose_pressure)
    arriving, pressure = reverse_core._peel_moving_side(model.geometry, sound - share, model.hop, True, model.end_extra,
                                                        0.0, reverse_core.JET_SEED, model.port_node, returned,
                                                        model.port_y, -1)
    return arriving, pressure
