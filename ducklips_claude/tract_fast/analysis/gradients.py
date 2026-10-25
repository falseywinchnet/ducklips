"""Exact Jacobians of a word's equation error by forward-mode derivatives
(the derivative toolchain: tract_fast/ad, built into build_ad by
`make -C ad`; findings in tract_fast/NOTES_gradients.md).

The residual of sequence_recover.sequence_residual (no seal, no nose, no
jet: a keyframed vowel sequence) is

    peel the sound through the hypothesis's moving tube -> implied flow U and
    the first node's pressure p; the returning wave w = (p - U rho c / A0) / 2;
    our folds forced by w (sustain FORCEWAVE) -> flow D; residual D - U over
    the target's window, mean removed.

OBSERVED (src/voice_body.cpp, VoiceBody::step): forced folds read only two
per-step quantities of the hypothesis, w and the first node's area A0; the
vocal tube, nose and branches step but never reach them. So the forced
render is replaced by a copy of the lungs, trachea and cord alone
(ad/forced_folds.hpp), in-process. The geometry and peel are copies of
reverse_core's (ad/geometry_peel.hpp). Both are templates over a scalar:
double for the residual, Dual<N> for N Jacobian columns per pass.

MEASURED (AA,IY,UW, 143 parameters, 39,360 steps; the Mini shared, load
~30): the native residual equals the reference's to 5e-12 relative (the
geometry and peel to the bit; the folds to 1 ulp at first) at 0.134 s CPU
against 0.303 s; the exact Jacobian costs 4.5 s CPU (33 native residuals,
15 reference ones) against 43.6 s CPU for the search's forward differences,
and agrees with central differences to their own error floor (1e-6).

The parameters reach the tube only through the AREASEQ schedule (rest areas
per 1 ms tick, linear between keyframes); its tangent is formed here in
closed form (keyframe values and times: keyframe_tangent), then carried
through geometry, peel and folds. A column moving one node's log area in one
phone leaves every other node's geometry unchanged, enters the peel only at
that node's layer and only from its phone's first keyframe interval; the
native code starts each block of columns there.

usage (in a search): SequenceJacobian(shape, target).jacobian(parameters) in
place of tube_recover.Problem.jacobian, or as least_squares(jac=...);
GRADIENT_THREADS=2 shares the blocks among native threads."""
import ctypes
import os
import subprocess
import sys
import tempfile

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import posture_recover
import reverse_core
import sequence_recover
import tube_recover

ROOT = posture_recover.ROOT
LIBRARY = os.path.join(ROOT, "build_ad", "libtract_ad.dylib")
SCRATCH = "/tmp/esf/grad"   # this toolchain's scratch (never /tmp/es, the reference's)
STEP_RATE = reverse_core.STEP_RATE
RHO_C = reverse_core.RHO_C
PITCH = 120.0               # sustain's arguments as posture_recover.render passes them
LUNG = 800.0
BLOCK = 16                  # tangent directions per native pass of the folds (CPU s for a 3-vowel word's
                            # 143 columns, the Mini at load ~30: folds 2.29 / 1.97 / 2.10 at 8 / 16 / 32)
PEEL_BLOCK = 8              # ... of the peel (2.37 / 2.70 / 3.59 s at 8 / 16 / 32: narrower blocks
                            # start nearer the glottis)
# 1: the link resistances glided as the engine's Tube does (one step ahead of
# reverse_core's _moving_geometry; an experiment: the residual is then no
# longer reverse_core's).
ENGINE_PHASE = int(os.environ.get("GRADIENT_ENGINE_PHASE", "0"))
THREADS = int(os.environ.get("GRADIENT_THREADS", "1"))   # native threads over the blocks (the Mini is shared:
                                                          # at most 2-3)

_library = None


def library():
    """The native library (ad/Makefile), loaded once."""
    global _library
    if _library is None:
        handle = ctypes.CDLL(LIBRARY)
        pointer = ctypes.POINTER(ctypes.c_double)
        integers = ctypes.POINTER(ctypes.c_int32)
        handle.sj_peel.restype = ctypes.c_int
        handle.sj_peel.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.c_int, pointer, pointer, ctypes.c_double,
                                   ctypes.c_double, ctypes.c_double, pointer, ctypes.c_int, pointer, ctypes.c_long,
                                   integers, ctypes.c_int, ctypes.c_int, ctypes.c_int, pointer, pointer, pointer, pointer,
                                   pointer, pointer]
        handle.sj_geometry.restype = ctypes.c_int
        handle.sj_geometry.argtypes = [ctypes.c_long, ctypes.c_int, pointer, pointer, ctypes.c_double, ctypes.c_int,
                                       pointer]
        handle.ff_primal.restype = ctypes.c_int
        handle.ff_primal.argtypes = [ctypes.c_long, ctypes.c_int, ctypes.c_double, ctypes.c_double, pointer, pointer,
                                     pointer, ctypes.c_long, pointer]
        handle.ff_tangent.restype = ctypes.c_int
        handle.ff_tangent.argtypes = [ctypes.c_long, ctypes.c_int, ctypes.c_double, ctypes.c_double, pointer, pointer,
                                      pointer, pointer, pointer, pointer, ctypes.c_long, ctypes.c_int, integers,
                                      ctypes.c_int, ctypes.c_int, pointer, pointer]
        _library = handle
    return _library


def integer_address(array):
    """An int32 C-contiguous array's pointer."""
    return array.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))


def address(array):
    """A float64 C-contiguous array's pointer (None passes NULL)."""
    if array is None:
        return None
    return array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def contiguous(array):
    return numpy.ascontiguousarray(array, dtype=numpy.float64)


# ---- The reference: sustain as posture_recover.render runs it ----

def sustain(tube_and_length, seconds, keyframes=None, forced=None, folds=None, geodump=None):
    """posture_recover.render's sustain run (the subset a vowel sequence
    uses), its scratch under SCRATCH: (flow entering or forced, radiated);
    with `geodump`, sustain's GEODUMP [step, node, 12] as well."""
    tube, length = tube_and_length
    os.makedirs(SCRATCH, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=SCRATCH) as folder:
        area_file = os.path.join(folder, "tube.txt")
        raw = os.path.join(folder, "tube.raw")
        with open(area_file, "w") as handle:
            handle.write("%.17g\n" % posture_recover.NODE)
            for area in tube:
                handle.write("%.17g\n" % area)
        environment = dict(os.environ, BRANCHES="", ASPIRATION="0", FRICATION="0", FLUX="0", FOLDTURB="0", CLOSING="0",
                           SKIN="0", AXIAL="0", WALLVAR="0", STATICAREA="1", VISCOUS="%g" % posture_recover.VISCOUS,
                           HOPLOSS="%g" % posture_recover.HOP, AREAFILE=area_file, RAW=raw, RAWDOUBLE="1",
                           ACTIVELENGTH="%.17g" % length)
        if keyframes is not None:
            with open(os.path.join(folder, "sequence.txt"), "w") as handle:
                for time, frame in zip(keyframes[0], keyframes[1]):
                    handle.write("%.17g " % time + " ".join("%.17g" % area for area in frame) + "\n")
            environment["AREASEQ"] = os.path.join(folder, "sequence.txt")
        if folds is not None:
            with open(os.path.join(folder, "folds.txt"), "w") as handle:
                for time, posture in zip(folds[0], folds[1]):
                    handle.write("%.17g " % time + " ".join("%.17g" % value for value in posture) + "\n")
            environment["FOLDKEYS"] = os.path.join(folder, "folds.txt")
        if forced is not None:
            numpy.asarray(forced, numpy.float64).tofile(os.path.join(folder, "forced.raw"))
            environment["FORCEWAVE"] = os.path.join(folder, "forced.raw")
        if geodump is not None:
            environment["GEODUMP"] = os.path.join(folder, "geometry.raw")
        subprocess.run([posture_recover.SUSTAIN, "aa", os.path.join(folder, "tube.wav"), "%g" % PITCH, "%g" % LUNG,
                        "%g" % seconds], env=environment, capture_output=True, cwd=ROOT)
        pairs = numpy.fromfile(raw, dtype=numpy.float64).reshape(-1, 2)
        geometry = None
        if geodump is not None:
            geometry = reverse_core.read_geometry(os.path.join(folder, "geometry.raw"), len(tube))
    return pairs[:, 0], pairs[:, 1], geometry


def render_steps(seconds):
    """The steps sustain renders for `seconds` as render passes it ("%g")."""
    return int(float("%g" % seconds) * STEP_RATE)


# ---- The native pieces ----

def even_blocks(directions, block=BLOCK):
    """Block starts cutting `directions` into blocks of `block`."""
    return numpy.array(list(range(0, directions, block)) + [directions], dtype=numpy.int32)


def peel(schedule, count, flow_out, end_extra, rest_dot=None, block_starts=None, threads=THREADS):
    """Geometry and peel (ad/peel_api.cpp): implied flow, first node's
    pressure (count steps), first node's area (the schedule's steps); with
    rest_dot [tick, node, direction], their tangents (blocks: block_starts,
    default even blocks of BLOCK)."""
    schedule = contiguous(schedule)
    steps_total, n = schedule.shape
    materials = contiguous(reverse_core.material_table(n))
    flow_out = contiguous(flow_out[:count])
    implied = numpy.zeros(count)
    pressure = numpy.zeros(count)
    area0 = numpy.zeros(steps_total)
    if rest_dot is None:
        code = library().sj_peel(count, steps_total, n, address(schedule), address(materials), posture_recover.VISCOUS,
                                 posture_recover.HOP, end_extra, address(flow_out), 0, None, 0, None, 0, 1,
                                 ENGINE_PHASE, address(implied), address(pressure), address(area0), None, None, None)
        assert code == 0
        return implied, pressure, area0
    rest_dot = contiguous(rest_dot)
    directions = rest_dot.shape[2]
    if block_starts is None:
        block_starts = even_blocks(directions)
    block_starts = numpy.ascontiguousarray(block_starts, dtype=numpy.int32)
    implied_dot = numpy.empty((count, directions))
    pressure_dot = numpy.empty((count, directions))
    area0_dot = numpy.empty((steps_total, directions))
    code = library().sj_peel(count, steps_total, n, address(schedule), address(materials), posture_recover.VISCOUS,
                             posture_recover.HOP, end_extra, address(flow_out), directions, address(rest_dot),
                             rest_dot.shape[0], integer_address(block_starts), len(block_starts) - 1, threads,
                             ENGINE_PHASE, address(implied), address(pressure), address(area0),
                             address(implied_dot), address(pressure_dot), address(area0_dot))
    assert code == 0
    return implied, pressure, area0, implied_dot, pressure_dot, area0_dot


def native_geometry(schedule):
    """GEODUMP's first three columns from the native copy of _moving_geometry."""
    schedule = contiguous(schedule)
    steps_total, n = schedule.shape
    materials = contiguous(reverse_core.material_table(n))
    out = numpy.zeros((steps_total, n, 3))
    library().sj_geometry(steps_total, n, address(schedule), address(materials), posture_recover.VISCOUS, ENGINE_PHASE,
                          address(out))
    return out


def forced_folds(waves, area0, vocal_count, waves_dot=None, area0_dot=None, block_starts=None, threads=THREADS):
    """The forced folds' flow (ad/forced_folds_api.cpp), and with tangents
    [step, direction] of the waves and the first node's area, its tangents
    (blocks: block_starts, default even blocks of BLOCK)."""
    waves = contiguous(waves)
    area0 = contiguous(area0)
    steps = len(waves)
    flow = numpy.zeros(steps)
    if waves_dot is None:
        code = library().ff_primal(steps, vocal_count, PITCH, LUNG, address(waves), address(area0), None, 0,
                                   address(flow))
        assert code == 0
        return flow
    waves_dot = contiguous(waves_dot)
    area0_dot = contiguous(area0_dot)
    directions = waves_dot.shape[1]
    if block_starts is None:
        block_starts = even_blocks(directions)
    block_starts = numpy.ascontiguousarray(block_starts, dtype=numpy.int32)
    flow_dot = numpy.empty((steps, directions))
    code = library().ff_tangent(steps, vocal_count, PITCH, LUNG, address(waves), address(waves_dot), address(area0),
                                address(area0_dot), None, None, 0, directions, integer_address(block_starts),
                                len(block_starts) - 1, threads, address(flow), address(flow_dot))
    assert code == 0
    return flow, flow_dot


# ---- The schedule's tangent ----

def time_jacobian(shape, parameters):
    """d(arrival_k, leaving_k)/d(time parameters) through SequenceShape.split's
    GAP clamps (max picks its first argument on a tie, as Python's does)."""
    phones = shape.phones
    count = shape.time_count
    times = parameters[phones * shape.count:phones * shape.count + count]
    arrivals = numpy.zeros(phones)
    leavings = numpy.zeros(phones)
    d_arrivals = numpy.zeros((phones, count))
    d_leavings = numpy.zeros((phones, count))
    for index in range(phones):
        if index > 0:
            floor = leavings[index - 1] + sequence_recover.GAP
            if times[2 * index - 1] >= floor:
                arrivals[index] = times[2 * index - 1]
                d_arrivals[index, 2 * index - 1] = 1.0
            else:
                arrivals[index] = floor
                d_arrivals[index] = d_leavings[index - 1]
        floor = arrivals[index] + sequence_recover.GAP
        if times[2 * index] >= floor:
            leavings[index] = times[2 * index]
            d_leavings[index, 2 * index] = 1.0
        else:
            leavings[index] = floor
            d_leavings[index] = d_arrivals[index]
    return d_arrivals, d_leavings


class ColumnOrder:
    """The native order of the Jacobian's columns and their blocks (each
    block one pass of the native code): the times, then each phone's log
    areas, the lips' node first, in blocks of at most `block`. A block of one
    phone's nodes starts the peel at its outermost node's layer and both the
    peel and the folds at that phone's first keyframe interval; mixing
    phones in a block would start it at the earliest phone's."""

    def __init__(self, shape, block=BLOCK):
        phones = shape.phones
        count = shape.count
        times = shape.time_count
        self.size = phones * count + times
        self.to_parameter = numpy.zeros(self.size, dtype=int)
        self.time_column = numpy.arange(times)
        self.tube_column = numpy.zeros((phones, count), dtype=int)
        for q in range(times):
            self.to_parameter[q] = phones * count + q
        for k in range(phones):
            base = times + k * count
            for node in range(count):
                column = base + (count - 1 - node)
                self.tube_column[k, node] = column
                self.to_parameter[column] = k * count + node
        self.phones = phones
        self.count = count
        self.times = times
        self.block_starts = self.starts(block)

    def starts(self, block):
        """Block starts in blocks of at most `block`, never mixing groups
        (the times; each phone's nodes)."""
        starts = [0]
        for first in range(0, self.times, block):
            starts.append(min(first + block, self.times))
        for k in range(self.phones):
            base = self.times + k * self.count
            for first in range(0, self.count, block):
                starts.append(base + min(first + block, self.count))
        return numpy.array(starts, dtype=numpy.int32)


def keyframe_tangent(times, frames, d_times, d_frames, steps_total):
    """The AREASEQ schedule's tangent per tick (sustain's 1 ms ticks, as
    reverse_core.keyframe_schedule samples them: linear between keyframes,
    held outside) for keyframe times' tangents d_times [key, column] and
    keyframe areas' d_frames [key, node, column]: [tick, node, column].
    Any parametrization that reaches the tube through keyframes (a word's
    phones, a stop's release frames, a held posture) differentiates by
    giving these two."""
    times = numpy.asarray(times, float)
    frames = numpy.asarray(frames, float)
    keys = len(times)
    ticks = (steps_total + 95) // 96
    moments = numpy.arange(ticks) * 96 / STEP_RATE
    out = numpy.zeros((ticks, frames.shape[1], d_frames.shape[2]))
    interval = numpy.searchsorted(times, moments, side="right") - 1
    for tick in range(ticks):
        j = interval[tick]
        if j < 0:
            out[tick] = d_frames[0]
        elif j >= keys - 1:
            out[tick] = d_frames[keys - 1]
        else:
            span = times[j + 1] - times[j]
            share = (moments[tick] - times[j]) / span
            out[tick] = (1.0 - share) * d_frames[j] + share * d_frames[j + 1]
            # The keyframes' times: d value / d t_j and d t_{j+1}.
            by_times = ((share - 1.0) / span) * d_times[j] - (share / span) * d_times[j + 1]
            if numpy.any(by_times != 0.0):
                out[tick] += numpy.outer(frames[j + 1] - frames[j], by_times)
    return out


def sequence_keyframe_jacobians(shape, parameters, order):
    """sequence_recover.keyframes' times and areas differentiated by a
    SequenceShape's parameters (columns in `order`): (times, frames, d_times,
    d_frames). Vowel sequences (each phone's tube at its arrival and its
    leaving; a stop's release frames are not yet differentiated)."""
    if shape.stops:
        raise NotImplementedError("stops add release keyframes: see NOTES_gradients.md")
    phones = shape.phones
    count = shape.count
    tubes, arrivals, leavings = shape.split(parameters)
    times, frames = sequence_recover.keyframes(tubes, arrivals, leavings)
    keys = len(times)
    assert keys == 2 * phones
    logs = parameters[:phones * count].reshape(phones, count)
    inside = (logs > numpy.log(0.02)) & (logs < numpy.log(30.0))
    d_tube = numpy.where(inside, tubes, 0.0)          # d area / d log area (the clip's)
    d_arrivals, d_leavings = time_jacobian(shape, parameters)
    d_times = numpy.zeros((keys, order.size))
    d_frames = numpy.zeros((keys, count, order.size))
    nodes = numpy.arange(count)
    for k in range(phones):
        d_times[2 * k, order.time_column] = d_arrivals[k]
        d_times[2 * k + 1, order.time_column] = d_leavings[k]
        for key in (2 * k, 2 * k + 1):
            d_frames[key, nodes, order.tube_column[k]] = d_tube[k]
    return times, frames, d_times, d_frames


def schedule_tangent(shape, parameters, order, steps_total):
    """The AREASEQ schedule's tangent for a SequenceShape's parameters:
    [tick, node, column] in `order`."""
    times, frames, d_times, d_frames = sequence_keyframe_jacobians(shape, parameters, order)
    return keyframe_tangent(times, frames, d_times, d_frames, steps_total)


# ---- The residual and its Jacobian ----

class SequenceJacobian:
    """sequence_recover.sequence_residual (target.equation, no seal, nose or
    jet; its bend rows too) in-process, and its exact Jacobian."""

    def __init__(self, shape, target, block=BLOCK, threads=THREADS, peel_block=PEEL_BLOCK):
        self.shape = shape
        self.target = target
        self.block = block
        self.threads = threads
        self.order = ColumnOrder(shape, block)
        self.peel_starts = self.order.starts(peel_block)

    def _layout(self, parameters):
        tubes, arrivals, leavings = self.shape.split(parameters)
        times, frames = sequence_recover.keyframes(tubes, arrivals, leavings)
        count = len(self.target.flow_out)
        steps = render_steps(self.target.seconds)
        steps_total = max(count, steps)
        schedule = reverse_core.keyframe_schedule(times, frames, steps_total)
        end_extra = self.shape.length - frames.shape[1] * posture_recover.NODE
        return count, steps, schedule, end_extra, frames.shape[1]

    def parts(self, parameters):
        """The implied flow, the forced folds' flow, the forced waves and the
        first node's area (primal)."""
        count, steps, schedule, end_extra, n = self._layout(parameters)
        implied, pressure, area0 = peel(schedule, count, self.target.flow_out, end_extra)
        forced = 0.5 * (pressure - implied * RHO_C / area0[:count])
        waves = numpy.zeros(steps)
        shared = min(count, steps)
        waves[:shared] = forced[:shared]
        driven = forced_folds(waves, area0[:steps], n)
        return implied, driven, forced, area0

    def residual(self, parameters):
        implied, driven, _, _ = self.parts(parameters)
        size = min(len(driven), len(implied))
        difference = driven[:size] - implied[:size]
        if self.target.filter is not None:
            difference = tube_recover.sosfiltfilt(self.target.filter, difference)
        residual = tube_recover.fixed_clip(difference, self.target)
        return self._with_bend(residual, parameters)

    def _with_bend(self, residual, parameters):
        bend = getattr(self.target, "bend", 0.0)
        if bend > 0.0:
            logs = parameters[:self.shape.phones * self.shape.count].reshape(self.shape.phones, self.shape.count)
            residual = numpy.concatenate([residual, bend * numpy.diff(logs, 2, axis=1).ravel()])
        return residual

    def jacobian(self, parameters, with_residual=False):
        """d residual / d parameters (rows: the residual's; columns: the
        parameters' order)."""
        count, steps, schedule, end_extra, n = self._layout(parameters)
        rest_dot = schedule_tangent(self.shape, parameters, self.order, schedule.shape[0])
        implied, pressure, area0, implied_dot, pressure_dot, area0_dot = peel(
            schedule, count, self.target.flow_out, end_extra, rest_dot, self.peel_starts, self.threads)
        forced = 0.5 * (pressure - implied * RHO_C / area0[:count])
        forced_dot = 0.5 * (pressure_dot - (RHO_C / area0[:count])[:, None] * implied_dot
                            + (implied * RHO_C / area0[:count] ** 2)[:, None] * area0_dot[:count])
        shared = min(count, steps)
        waves = numpy.zeros(steps)
        waves[:shared] = forced[:shared]
        waves_dot = numpy.zeros((steps, self.order.size))
        waves_dot[:shared] = forced_dot[:shared]
        driven, driven_dot = forced_folds(waves, area0[:steps], n, waves_dot, area0_dot[:steps],
                                          self.order.block_starts, self.threads)
        size = min(len(driven), len(implied))
        start = self.target.start
        length = self.target.length
        difference_dot = driven_dot[:size] - implied_dot[:size]
        if self.target.filter is not None:
            difference_dot = tube_recover.sosfiltfilt(self.target.filter, difference_dot, axis=0)
        clip = difference_dot[start:start + length]
        if len(clip) < length:
            clip = numpy.concatenate([clip, numpy.zeros((length - len(clip), clip.shape[1]))])
        clip = (clip - clip.mean(axis=0)) / self.target.scale
        jacobian = numpy.zeros((clip.shape[0], self.order.size))
        jacobian[:, self.order.to_parameter] = clip
        bend = getattr(self.target, "bend", 0.0)
        if bend > 0.0:
            phones = self.shape.phones
            nodes = self.shape.count
            rows = []
            for k in range(phones):
                for i in range(nodes - 2):
                    row = numpy.zeros(self.order.size)
                    row[k * nodes + i] += bend
                    row[k * nodes + i + 1] -= 2.0 * bend
                    row[k * nodes + i + 2] += bend
                    rows.append(row)
            jacobian = numpy.vstack([jacobian, numpy.array(rows)])
        if with_residual:
            difference = driven[:size] - implied[:size]
            if self.target.filter is not None:
                difference = tube_recover.sosfiltfilt(self.target.filter, difference)
            residual = self._with_bend(tube_recover.fixed_clip(difference, self.target), parameters)
            return jacobian, residual
        return jacobian


if __name__ == "__main__":
    print(__doc__)
