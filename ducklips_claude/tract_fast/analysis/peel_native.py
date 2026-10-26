"""The tube peeled from the sound in the native time-major kernel
(bindings/peel_core.cpp): observer.StatesPeel's and reverse_core.SidePeel's
values to the bit, with their interfaces (advance, reach, resumable,
resumed, boundary), so observer._implied_run and adjoint.ResidualCache use
either. Record and measurements: NOTES_residual_speed.md.

The kernel runs every layer (link j and node j) in one wavefront loop, layer
j's link d_j = 1 + 2 (n - 2 - j) wavefronts behind the open end; a
wavefront W reads inputs (sound, geometry rows, gate) of steps below W+1
only. `wave` is the next wavefront to run; `boundary` the steps whose
inputs the peel has read (min(wave, count)), the numba classes' meaning of
it. Node 0 is final below wave - 2 n + 3.

Snapshots (with snapshots=True / record=True) are taken at every wavefront
that is a multiple of `chunk`: a peel whose inputs first differ at step k
resumes from the last one at or before k."""
import ctypes
import os
from typing import Optional

import numpy

import kernels_jet
import reverse_core

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIBRARY = os.environ.get("PEEL_CORE_LIB", os.path.join(ROOT, "build_engine", "libpeel_core.dylib"))
VARIANT_STATES = 0
VARIANT_SIDE = 1
CHUNK = 512          # wavefronts between snapshots (the numba peels' chunk)
STATE_COLUMNS = 19   # observer.STATE_COLUMNS

_doubles = ctypes.POINTER(ctypes.c_double)
_longs = ctypes.POINTER(ctypes.c_int64)


class PeelInputs(ctypes.Structure):
    """bindings/peel_core.h's peel_inputs (field for field)."""
    _fields_ = [("geometry", _doubles), ("sound", _doubles), ("steps", ctypes.c_int64), ("nodes", ctypes.c_int64),
                ("t_start", ctypes.c_int64), ("record", ctypes.c_int64), ("init", _doubles),
                ("recorded", _doubles), ("flow_out", _doubles), ("pressure_out", _doubles),
                ("hop", ctypes.c_double), ("hop_zero", ctypes.c_double), ("end_extra", ctypes.c_double),
                ("jet_gain", ctypes.c_double), ("jet_closing", ctypes.c_double), ("jet_bound", ctypes.c_double),
                ("gate", _doubles), ("noise", _doubles),
                ("gate_length", ctypes.c_int64), ("noise_links", ctypes.c_int64), ("noise_steps", ctypes.c_int64),
                ("jet_seed", ctypes.c_int64), ("kinetic", ctypes.c_int32), ("has_gate", ctypes.c_int32),
                ("variant", ctypes.c_int32), ("reserved", ctypes.c_int32)]


class GeometryInputs(ctypes.Structure):
    """bindings/peel_core.h's geometry_inputs (field for field)."""
    _fields_ = [("rest", _doubles), ("materials", _doubles), ("out", _doubles), ("state", _doubles),
                ("counters", _doubles), ("snaps", _doubles), ("snap_counters", _doubles), ("least", _doubles),
                ("first_area", _doubles), ("jet_node", _longs), ("jet_area", _doubles), ("steps", ctypes.c_int64),
                ("nodes", ctypes.c_int64), ("begin", ctypes.c_int64), ("end", ctypes.c_int64),
                ("snap_every", ctypes.c_int64), ("rest_every", ctypes.c_int64), ("viscous", ctypes.c_double)]


class ForwardInputs(ctypes.Structure):
    """bindings/peel_core.h's forward_inputs (field for field)."""
    _fields_ = [("geometry", _doubles), ("flow_in", _doubles), ("state", _doubles), ("lips_out", _doubles),
                ("port_pressure", _doubles), ("gate", _doubles), ("noise", _doubles), ("steps", ctypes.c_int64),
                ("nodes", ctypes.c_int64), ("t0", ctypes.c_int64), ("t1", ctypes.c_int64),
                ("gate_length", ctypes.c_int64), ("noise_links", ctypes.c_int64), ("noise_steps", ctypes.c_int64),
                ("jet_seed", ctypes.c_int64), ("hop", ctypes.c_double), ("hop_zero", ctypes.c_double),
                ("end_extra", ctypes.c_double), ("seal", ctypes.c_double), ("jet_gain", ctypes.c_double),
                ("jet_closing", ctypes.c_double), ("jet_bound", ctypes.c_double), ("kinetic", ctypes.c_int32),
                ("seal_blend", ctypes.c_int32), ("has_gate", ctypes.c_int32), ("reserved", ctypes.c_int32)]


_lib = None


def library() -> ctypes.CDLL:
    """The kernel library (loaded once a process)."""
    global _lib
    if _lib is not None:
        return _lib
    lib = ctypes.CDLL(LIBRARY)
    lib.peel_state_size.restype = ctypes.c_int64
    lib.peel_state_size.argtypes = [ctypes.c_int64]
    lib.peel_complete.restype = ctypes.c_int64
    lib.peel_complete.argtypes = [ctypes.c_int64, ctypes.c_int64]
    lib.peel_begin.restype = None
    lib.peel_begin.argtypes = [ctypes.POINTER(PeelInputs), _doubles]
    lib.peel_advance.restype = None
    lib.peel_advance.argtypes = [ctypes.POINTER(PeelInputs), _doubles, ctypes.c_int64, ctypes.c_int64, _doubles,
                                 ctypes.c_int64, ctypes.c_int64]
    lib.geometry_steps.restype = None
    lib.geometry_steps.argtypes = [ctypes.POINTER(GeometryInputs)]
    lib.forward_steps.restype = None
    lib.forward_steps.argtypes = [ctypes.POINTER(ForwardInputs)]
    _lib = lib
    return lib


def forward_steps(geometry: numpy.ndarray, flow_in: numpy.ndarray, t0: int, t1: int, state: numpy.ndarray,
                  hop: float, hop_zero: float, kinetic: bool, end_extra: float, seal: float, lips_out: numpy.ndarray,
                  port_pressure: numpy.ndarray, jet_gain: float, jet_seed: int, jet_closing: float,
                  jet_gate: Optional[numpy.ndarray], noise: numpy.ndarray, seal_blend: bool,
                  jet_bound: float) -> None:
    """observer._forward_steps (its arguments, less the jet's coefficient
    cache, which holds a pure function of the corner); the same values to
    the bit. geometry, state, lips_out, port_pressure C contiguous float64."""
    flow = numpy.ascontiguousarray(flow_in, float)
    has_gate = jet_gate is not None
    gate = numpy.ascontiguousarray(jet_gate, float) if has_gate else numpy.zeros(1)
    table = numpy.ascontiguousarray(noise, float) if noise.size > 0 else numpy.zeros((1, 1))
    record = ForwardInputs()
    record.geometry = _pointer(geometry)
    record.flow_in = _pointer(flow)
    record.state = _pointer(state)
    record.lips_out = _pointer(lips_out)
    record.port_pressure = _pointer(port_pressure)
    record.gate = _pointer(gate)
    record.noise = _pointer(table)
    record.steps = int(geometry.shape[0])
    record.nodes = int(geometry.shape[1])
    record.t0 = int(t0)
    record.t1 = int(t1)
    record.gate_length = int(gate.shape[0])
    record.noise_links = int(table.shape[0]) if noise.size > 0 else 0
    record.noise_steps = int(table.shape[1]) if noise.size > 0 else 0
    record.jet_seed = int(jet_seed)
    record.hop = float(hop)
    record.hop_zero = float(hop_zero)
    record.end_extra = float(end_extra)
    record.seal = float(seal)
    record.jet_gain = float(jet_gain)
    record.jet_closing = float(jet_closing)
    record.jet_bound = float(jet_bound)
    record.kinetic = 1 if kinetic else 0
    record.seal_blend = 1 if seal_blend else 0
    record.has_gate = 1 if has_gate else 0
    record.reserved = 0
    library().forward_steps(ctypes.byref(record))


class GeometryNotes:
    """What the native geometry gathers per step besides its rows: the
    least area, the first node's area, the jet's narrowest node and area."""

    def __init__(self, steps: int) -> None:
        self.least = numpy.zeros(steps)
        self.first_area = numpy.zeros(steps)
        self.jet_node = numpy.full(steps, -1, dtype=numpy.int64)
        self.jet_area = numpy.ones(steps)


def geometry_steps(rest: numpy.ndarray, materials: numpy.ndarray, viscous: float, out: numpy.ndarray,
                   state: numpy.ndarray, counters: numpy.ndarray, begin: int, end: int, every: int,
                   snaps: numpy.ndarray, snap_counters: numpy.ndarray, notes: GeometryNotes,
                   rest_every: int = 1) -> None:
    """reverse_core._geometry_steps (its arguments; the same values to the
    bit), and the per-step notes of steps begin..end - 1. Every array C
    contiguous float64 (jet_node int64), as _geometry_steps takes them;
    `rest` may hold a row per rest_every steps (the command's ticks: the
    rows numpy.repeat would make, not made)."""
    record = GeometryInputs()
    record.rest = _pointer(rest)
    record.materials = _pointer(materials)
    record.out = _pointer(out)
    record.state = _pointer(state)
    record.counters = _pointer(counters)
    record.snaps = _pointer(snaps)
    record.snap_counters = _pointer(snap_counters)
    record.least = _pointer(notes.least)
    record.first_area = _pointer(notes.first_area)
    record.jet_node = notes.jet_node.ctypes.data_as(_longs)
    record.jet_area = _pointer(notes.jet_area)
    record.steps = int(out.shape[0])
    record.nodes = int(rest.shape[1])
    record.begin = int(begin)
    record.end = int(end)
    record.snap_every = int(every)
    record.rest_every = int(rest_every)
    record.viscous = float(viscous)
    library().geometry_steps(ctypes.byref(record))


def _pointer(array: Optional[numpy.ndarray]) -> Optional[ctypes.POINTER(ctypes.c_double)]:
    if array is None:
        return None
    return array.ctypes.data_as(_doubles)


class NativePeel:
    """The state, inputs and snapshots of one native peel. The arrays it
    borrows (geometry, sound, flow/pressure outputs, recorded, init, gate,
    noise) are held as fields, so they live as long as it."""

    def __init__(self, variant: int, geometry: numpy.ndarray, sound: numpy.ndarray, t_start: int,
                 init: Optional[numpy.ndarray], record: int, hop: float, hop_zero: float, kinetic: bool,
                 end_extra: float, flow_out: numpy.ndarray, pressure_out: numpy.ndarray,
                 recorded: Optional[numpy.ndarray], jet_gain: float, jet_seed: int, jet_closing: float,
                 jet_gate: Optional[numpy.ndarray], snapshots: bool, chunk: int) -> None:
        self.variant = int(variant)
        self.geometry = numpy.ascontiguousarray(geometry, float)
        self.sound = numpy.ascontiguousarray(sound, float)
        count = self.sound.shape[0]
        n = self.geometry.shape[1]
        self.count = count
        self.nodes = n
        self.t_start = int(t_start)
        self.init = None if init is None else numpy.ascontiguousarray(init, float)
        self.record = int(record)
        self.hop = float(hop)
        self.hop_zero = float(hop_zero)
        self.kinetic = bool(kinetic)
        self.end_extra = float(end_extra)
        self.flow_out = flow_out
        self.pressure_out = pressure_out
        self.recorded = recorded if recorded is not None else numpy.zeros((n, STATE_COLUMNS))
        self.jet_gain = float(jet_gain)
        self.jet_seed = int(jet_seed)
        self.jet_closing = float(jet_closing)
        self.has_gate = jet_gate is not None
        self.gate = numpy.ascontiguousarray(jet_gate, float) if self.has_gate else numpy.zeros(1)
        self.chunk = int(chunk)
        if self.jet_gain > 0.0:
            self.noise = kernels_jet.noise_table(self.jet_seed, n - 1, count)
        else:
            self.noise = numpy.zeros((1, 1))
        lib = library()
        self.size = int(lib.peel_state_size(n))
        self.complete = int(lib.peel_complete(count, n))
        self.state = numpy.zeros(self.size)
        self.wave = self.t_start
        self.snap_every = self.chunk if snapshots else 0
        slots = self.complete // self.chunk + 1 if snapshots else 1
        self.snaps = numpy.zeros((slots, self.size))
        lib.peel_begin(ctypes.byref(self.inputs()), _pointer(self.state))
        if snapshots and self.wave % self.chunk == 0:
            self.snaps[self.wave // self.chunk] = self.state

    @property
    def boundary(self) -> int:
        """The steps whose inputs the peel has read."""
        return min(self.wave, self.count)

    def inputs(self) -> PeelInputs:
        """The C record over this peel's fields (valid while they are)."""
        record = PeelInputs()
        record.geometry = _pointer(self.geometry)
        record.sound = _pointer(self.sound)
        record.steps = self.count
        record.nodes = self.nodes
        record.t_start = self.t_start
        record.record = self.record
        record.init = _pointer(self.init)
        record.recorded = _pointer(self.recorded)
        record.flow_out = _pointer(self.flow_out)
        record.pressure_out = _pointer(self.pressure_out)
        record.hop = self.hop
        record.hop_zero = self.hop_zero
        record.end_extra = self.end_extra
        record.jet_gain = self.jet_gain
        record.jet_closing = self.jet_closing
        record.jet_bound = float(reverse_core.JET_BOUND)
        record.gate = _pointer(self.gate)
        record.noise = _pointer(self.noise)
        record.gate_length = int(self.gate.shape[0])
        record.noise_links = int(self.noise.shape[0])
        record.noise_steps = int(self.noise.shape[1])
        record.jet_seed = self.jet_seed
        record.kinetic = 1 if self.kinetic else 0
        record.has_gate = 1 if self.has_gate else 0
        record.variant = self.variant
        record.reserved = 0
        return record

    def reach(self, until: int) -> int:
        """The wavefront to which a peel wanting node 0's values below
        `until` and its record must run (node 0 runs step t at wavefront
        t + 2 n - 2; the record's last pieces at record + 2 n - 2)."""
        n = self.nodes
        target = max(int(until) + 2 * n - 2, self.record + 2 * n - 1, self.t_start)
        return min(self.complete, target)

    def advance(self, target: int) -> None:
        """Run to wavefront `target`; `count` or more: complete."""
        target = int(target)
        if target >= self.count:
            target = self.complete
        if target <= self.wave:
            return
        library().peel_advance(ctypes.byref(self.inputs()), _pointer(self.state), self.wave, target,
                               _pointer(self.snaps), self.snap_every, self.snaps.shape[0])
        self.wave = target

    def resume_wave(self, first_changed: int) -> int:
        """The snapshot wavefront a resumed peel would start from."""
        if self.snap_every <= 0:
            return self.t_start
        wave = min(int(first_changed), self.wave, self.count - 1) // self.chunk * self.chunk
        return wave

    def copy_resumed(self, wave: int, geometry: numpy.ndarray, sound: numpy.ndarray,
                     jet_gate: Optional[numpy.ndarray], flow_out: numpy.ndarray, pressure_out: numpy.ndarray,
                     recorded: Optional[numpy.ndarray]) -> "NativePeel":
        """A peel at snapshot `wave` for inputs equal to this one's before
        it (the caller's outputs hold this run's values)."""
        other = NativePeel.__new__(type(self))
        other.__dict__.update(self.__dict__)
        other.geometry = numpy.ascontiguousarray(geometry, float)
        other.sound = numpy.ascontiguousarray(sound, float)
        other.flow_out = flow_out
        other.pressure_out = pressure_out
        if recorded is not None:
            other.recorded = recorded
        if self.has_gate:
            other.gate = numpy.ascontiguousarray(jet_gate, float)
        other.state = self.snaps[wave // self.chunk].copy()
        other.snaps = self.snaps.copy()
        other.wave = wave
        return other


class NativeStatesPeel(NativePeel):
    """observer.StatesPeel's interface and values (the native kernel)."""

    def __init__(self, geometry, sound, t_start, init, record, hop, hop_zero, kinetic, end_extra, flow_in_out,
                 pressure_out, recorded, jet_gain=0.0, jet_seed=0, jet_closing=0.0, jet_gate=None, snapshots=False,
                 chunk=CHUNK):
        start = int(t_start)
        NativePeel.__init__(self, VARIANT_STATES, geometry, sound, start, init if start > 0 else None, record,
                            hop, hop_zero, kinetic, end_extra, flow_in_out, pressure_out, recorded, jet_gain,
                            jet_seed, jet_closing, jet_gate, snapshots, chunk)
        if self.init is None:
            self.init = numpy.zeros((self.nodes, STATE_COLUMNS))

    def resumable(self, first_changed: int) -> int:
        """The boundary a resumed peel would start from (t_start: none)."""
        wave = self.resume_wave(first_changed)
        return wave if wave > self.t_start else self.t_start

    def resumed(self, geometry, sound, first_changed, jet_gate, flow_in_out, pressure_out, recorded):
        """StatesPeel.resumed: None when no snapshot past t_start applies."""
        wave = self.resume_wave(first_changed)
        if self.snap_every <= 0 or wave <= self.t_start:
            return None
        return self.copy_resumed(wave, geometry, sound, jet_gate, flow_in_out, pressure_out, recorded)


class NativeSidePeel(NativePeel):
    """reverse_core.SidePeel's values for a tube without a side port or a
    stop node (peel_moving_record's peel): `arriving` (node 0's flow) and
    `pressure` (its pressure) are complete after advance(count)."""

    def __init__(self, geometry, flow_out, hop, kinetic, end_extra, jet_gain, jet_seed, jet_gate=None,
                 record=False, chunk=CHUNK):
        count = numpy.asarray(flow_out).shape[0]
        self.arriving = numpy.zeros(count)
        self.pressure = numpy.zeros(count)
        NativePeel.__init__(self, VARIANT_SIDE, geometry, flow_out, 0, None, -1, hop, 0.0, kinetic, end_extra,
                            self.arriving, self.pressure, None, jet_gain, jet_seed, 0.0, jet_gate, record, chunk)

    def resumed(self, geometry, flow_out, first_changed, jet_gate=None):
        """SidePeel.resumed (from wavefront 0 when no later snapshot applies)."""
        wave = self.resume_wave(first_changed)
        arriving = self.arriving.copy()
        pressure = self.pressure.copy()
        other = self.copy_resumed(wave, geometry, flow_out, jet_gate, arriving, pressure, None)
        other.arriving = arriving
        other.pressure = pressure
        return other


def peel_moving_record(geometry, flow_out, hop=0.0, kinetic=True, end_extra=0.0, jet=0.0,
                       seed=reverse_core.JET_SEED, gate=None):
    """reverse_core.peel_moving_record's values and a resumable record."""
    run = NativeSidePeel(geometry, flow_out, hop, kinetic, end_extra, jet, seed, gate, record=True)
    run.advance(run.count)
    return run.arriving.copy(), run.pressure.copy(), run


def peel_moving_resume(record, geometry, flow_out, first_changed, gate=None):
    """reverse_core.peel_moving_resume's values."""
    run = record.resumed(geometry, flow_out, first_changed, gate)
    run.advance(run.count)
    return run.arriving.copy(), run.pressure.copy(), run
