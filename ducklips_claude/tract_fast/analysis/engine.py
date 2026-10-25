"""The engine in-process: posture_recover.render without a subprocess.

render(...) takes exactly posture_recover.render's arguments and returns
exactly its values (the reference runs `sustain` with keyframe text files
and reads its RAW file; here the same run, tract::SustainRun, is called
through bindings/tract_engine.h with the numbers themselves). Where the
reference writes text with %.17g, the doubles are passed as they are (the
text reads back to the same double); where it writes %g (the pitch and the
seconds on the command line), the same rounding is applied here first.

Equivalence with the reference tool: bindings/equivalence.py (byte-identical
outputs over its matrix; NOTES_engine.md).

Lower level, for searches that re-run only the tail of a word:
    settings = prepare(tube_and_length, seconds, ...)   # render's arguments
    run = Run(settings)
    run.advance(k)                  # steps 0..k-1; returns their flows
    snap = run.snapshot()           # the whole state at step k
    tail = snap.advance(n, other)   # steps k..n-1 with changed inputs
"""
import ctypes
import os

import numpy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIBRARY = os.environ.get("TRACT_ENGINE_LIB", os.path.join(ROOT, "build_engine", "libtract_engine.dylib"))

# posture_recover's constants (the reference's values; bindings/equivalence.py
# checks them against the reference module).
NODE = 35000.0 / 96000.0
VISCOUS = 1.0
HOP = 0.02
JET_CLOSING = 0.04
NOSE_POLE = 0.8
NOSE_ZERO = 0.79

# The sustain tool's variables. A caller's environment reaches the reference
# tool through posture_recover.render's subprocess; those render does not set
# are passed through here the same way (parsed by the tool's own parser).
# TRACE only prints (the reference discards its output) and is dropped.
SUSTAIN_VARIABLES = (
    "AREAFILE", "AREAFILE2", "GLIDE", "AREASEQ", "VELUMSEQ", "FOLDSEQ", "JETSEQ", "FOLDKEYS", "AREATICKS", "VELUMTICKS", "FOLDTICKS", "PITCHEND", "RELEASE",
    "BRANCHES", "ASPIRATION", "FRICATION", "FLUX", "CLOSING", "TREMOR", "ADDTREMOR", "BREATHTREMOR", "PITCHFB",
    "PORTRESIST", "REVJET", "REVJETCLOSE", "SEALBLEND", "JETBOUND", "FOLDTURB", "SKIN", "WALLS", "WALLMASS", "AXIAL", "WALLVAR", "VISCOUS",
    "ACTIVELENGTH", "STATICAREA", "REVNOSE", "SEAL", "NOSERAD", "NOSEVISCOUS", "NOSESHELF", "NOSEHOP", "HOPLOSS",
    "REST", "CHINK", "FOLDDAMP", "VELUM", "SOURCE", "RAW", "FORCEWAVE", "OBSERVE", "FORCEFLOW", "RAWDOUBLE",
    "RAWCORD", "PORTWAVES", "PORTDUMP", "GEODUMP", "AREAS", "WALLDUMP", "NOSEDUMP")

# Variables that name files: the reference tool runs in ROOT, so a relative
# path is ROOT's.
FILE_VARIABLES = ("AREAFILE", "AREAFILE2", "AREASEQ", "VELUMSEQ", "FOLDSEQ", "JETSEQ", "FOLDKEYS", "AREATICKS",
                  "VELUMTICKS", "FOLDTICKS", "SOURCE", "RAW",
                  "FORCEWAVE", "FORCEFLOW", "PORTWAVES", "PORTDUMP", "GEODUMP", "AREAS", "WALLDUMP", "NOSEDUMP")

_lib = None


def _beside_root(path):
    return os.path.join(ROOT, path)


def _library():
    global _lib
    if _lib is not None:
        return _lib
    lib = ctypes.CDLL(LIBRARY)
    pointer = ctypes.c_void_p
    doubles = ctypes.POINTER(ctypes.c_double)
    lib.tract_error.restype = ctypes.c_char_p
    lib.tract_settings_new.restype = pointer
    lib.tract_settings_copy.restype = pointer
    lib.tract_settings_copy.argtypes = [pointer]
    lib.tract_settings_free.argtypes = [pointer]
    lib.tract_settings_number.argtypes = [pointer, ctypes.c_char_p, ctypes.c_double]
    lib.tract_settings_pair.argtypes = [pointer, ctypes.c_char_p, ctypes.c_double, ctypes.c_double]
    lib.tract_settings_text.argtypes = [pointer, ctypes.c_char_p, ctypes.c_char_p]
    lib.tract_settings_array.argtypes = [pointer, ctypes.c_char_p, doubles, ctypes.c_int64, ctypes.c_int64]
    lib.tract_settings_environment.argtypes = [pointer, ctypes.c_int, ctypes.POINTER(ctypes.c_char_p),
                                               ctypes.POINTER(ctypes.c_char_p)]
    lib.tract_settings_steps.restype = ctypes.c_int64
    lib.tract_settings_steps.argtypes = [pointer]
    lib.tract_settings_vocal_nodes.restype = ctypes.c_int64
    lib.tract_settings_vocal_nodes.argtypes = [pointer]
    if hasattr(lib, "tract_settings_divergence"):  # (absent from builds before it; benchmarks load those)
        lib.tract_settings_divergence.restype = ctypes.c_int64
        lib.tract_settings_divergence.argtypes = [pointer, pointer]
    lib.tract_run_new.restype = pointer
    lib.tract_run_new.argtypes = [pointer]
    lib.tract_run_copy.restype = pointer
    lib.tract_run_copy.argtypes = [pointer]
    lib.tract_run_free.argtypes = [pointer]
    lib.tract_run_step.restype = ctypes.c_int64
    lib.tract_run_step.argtypes = [pointer]
    lib.tract_run_total.restype = ctypes.c_int64
    lib.tract_run_total.argtypes = [pointer]
    lib.tract_run_advance.restype = ctypes.c_int64
    lib.tract_run_advance.argtypes = [pointer, pointer, ctypes.c_int64, doubles, doubles, doubles,
                                      ctypes.POINTER(ctypes.c_int32)]
    lib.tract_run_finish.argtypes = [pointer, pointer]
    lib.tract_run_nose.restype = ctypes.c_int64
    lib.tract_run_nose.argtypes = [pointer, doubles, doubles, ctypes.c_int64]
    lib.tract_render.restype = ctypes.c_int64
    lib.tract_render.argtypes = [pointer, doubles, doubles, doubles, ctypes.POINTER(ctypes.c_int32)]
    _lib = lib
    return lib


def _check(code):
    if code != 0:
        raise ValueError(_library().tract_error().decode())


def _pointer(array):
    if array is None:
        return None
    return array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


class Settings:
    """A render's inputs (tract_settings), owned."""

    def __init__(self, handle=None, lib=None):
        self.lib = lib if lib is not None else _library()
        self.handle = handle if handle is not None else self.lib.tract_settings_new()

    def __del__(self):
        if getattr(self, "handle", None):
            self.lib.tract_settings_free(self.handle)
            self.handle = None

    def copy(self):
        return Settings(self.lib.tract_settings_copy(self.handle), self.lib)

    def number(self, name, value):
        _check(self.lib.tract_settings_number(self.handle, name.encode(), float(value)))
        return self

    def pair(self, name, first, second):
        _check(self.lib.tract_settings_pair(self.handle, name.encode(), float(first), float(second)))
        return self

    def text(self, name, value):
        _check(self.lib.tract_settings_text(self.handle, name.encode(), value.encode()))
        return self

    def array(self, name, values, columns=None):
        values = numpy.ascontiguousarray(values, dtype=numpy.float64)
        if values.ndim == 1:
            values = values.reshape(1, -1) if columns is None else values.reshape(-1, columns)
        rows, width = values.shape if values.size else (0, columns if columns is not None else 0)
        _check(self.lib.tract_settings_array(self.handle, name.encode(), _pointer(values), rows, width))
        return self

    def environment(self, variables):
        names = (ctypes.c_char_p * len(variables))(*[name.encode() for name in variables])
        texts = (ctypes.c_char_p * len(variables))(*[variables[name].encode() for name in variables])
        _check(self.lib.tract_settings_environment(self.handle, len(variables), names, texts))
        return self

    def steps(self):
        return int(self.lib.tract_settings_steps(self.handle))

    def divergence(self, other):
        """The first step whose inputs differ from `other`'s (-1: they build
        different bodies): a run of `other` may resume from this one's
        state at any step up to it."""
        return int(self.lib.tract_settings_divergence(self.handle, other.handle))

    def vocal_nodes(self):
        return int(self.lib.tract_settings_vocal_nodes(self.handle))


class Run:
    """A render in progress (tract_run): advance it, snapshot it (a copy of
    the whole state), resume a snapshot with changed inputs."""

    def __init__(self, settings, handle=None):
        self.settings = settings
        self.lib = settings.lib
        self.handle = handle if handle is not None else self.lib.tract_run_new(settings.handle)
        if not self.handle:
            raise ValueError(self.lib.tract_error().decode())

    def __del__(self):
        if getattr(self, "handle", None):
            self.lib.tract_run_free(self.handle)
            self.handle = None

    @property
    def step(self):
        return int(self.lib.tract_run_step(self.handle))

    @property
    def total(self):
        # (The steps of the settings it is advanced with: a resumed run takes
        # its length from them.)
        return self.settings.steps()

    def snapshot(self):
        return Run(self.settings, self.lib.tract_run_copy(self.handle))

    def advance(self, end=None, settings=None, port_waves=False, geometry=False):
        """Steps to `end` (default: the render's end) with `settings` (default:
        the run's own; changed ones act from this step on). Returns the RAW
        pairs (steps, 2), and the PORTWAVES rows (steps, 6) and the GEODUMP
        rows (steps, nodes, 12) when asked, and whether every output was
        finite."""
        if settings is not None:
            self.settings = settings
        end = self.total if end is None else min(int(end), self.total)
        count = max(end - self.step, 0)
        raw = numpy.empty((count, 2))
        waves = numpy.empty((count, 6)) if port_waves else None
        nodes = self.settings.vocal_nodes()
        rows = numpy.empty((count, nodes, 12)) if geometry else None
        status = ctypes.c_int32(0)
        written = self.lib.tract_run_advance(self.handle, self.settings.handle, end, _pointer(raw), _pointer(waves),
                                               _pointer(rows), ctypes.byref(status))
        out = [raw[:written]]
        if port_waves:
            out.append(waves[:written])
        if geometry:
            out.append(rows[:written])
        out.append(status.value == 0)
        return tuple(out)

    def finish(self):
        self.lib.tract_run_finish(self.handle, self.settings.handle)

    def advance_into(self, end, raw, waves=None):
        """Steps to `end`, writing the RAW pairs into raw[step] and the
        PORTWAVES rows into waves[step] (arrays covering the whole render).
        Returns the step reached and whether every output was finite."""
        start = self.step
        status = ctypes.c_int32(0)
        doubles = ctypes.POINTER(ctypes.c_double)
        raw_at = ctypes.cast(raw.ctypes.data + 16 * start, doubles)
        waves_at = None if waves is None else ctypes.cast(waves.ctypes.data + 48 * start, doubles)
        written = self.lib.tract_run_advance(self.handle, self.settings.handle, int(end), raw_at, waves_at, None,
                                             ctypes.byref(status))
        return start + int(written), status.value == 0

    def nose(self):
        """NOSEDUMP's numbers for the run as it stands: the header (port node,
        node count, hop low-pass, viscous scale) and the rows (nodes, 12)."""
        header = numpy.zeros(4)
        count = int(self.lib.tract_run_nose(self.handle, _pointer(header), None, 0))
        rows = numpy.zeros((count, 12))
        self.lib.tract_run_nose(self.handle, _pointer(header), _pointer(rows), count)
        return header, rows


def _finite_text(values, what):
    # The reference writes these as text (%.17g); a non-finite value there
    # stops the tool's reader, which this path does not imitate.
    if not numpy.all(numpy.isfinite(numpy.asarray(values, dtype=numpy.float64))):
        raise ValueError("%s: non-finite values are written as text by the reference; not supported" % what)


def _stream(rows, width):
    """Keyframe rows as the tool reads their text: one stream of numbers
    taken `width` at a time (an incomplete tail dropped)."""
    tokens = []
    for row in rows:
        tokens.extend(float(value) for value in row)
    usable = len(tokens) // width * width
    return numpy.array(tokens[:usable], dtype=numpy.float64).reshape(-1, width)


def field_parts(field):
    """A control field's parts as float64 arrays (None where absent): the
    areas (ticks x nodes), the velum (ticks x 1), the fold postures (ticks x
    3: lower rest, upper rest, chink), one row a 1 ms control tick."""
    if field is None:
        return None, None, None
    areas = field.get("areas")
    velum = field.get("velum")
    folds = field.get("folds")
    areas = None if areas is None else numpy.ascontiguousarray(areas, dtype=numpy.float64).reshape(len(areas), -1)
    velum = None if velum is None else numpy.ascontiguousarray(velum, dtype=numpy.float64).reshape(-1, 1)
    folds = None if folds is None else numpy.ascontiguousarray(folds, dtype=numpy.float64).reshape(-1, 3)
    return areas, velum, folds


def prepare(tube_and_length, seconds=0.3, pitch=120.0, second=None, glide=None, forced=None, keyframes=None,
            wave=None, jet=0.0, seal=0.0, velum=None, nose_dump=None, port_waves=None, folds=None, observed=None,
            jet_gate=None, field=None):
    """posture_recover.render's inputs as engine Settings (port_waves is
    taken by render, which writes it; here it is ignored)."""
    tube, length = tube_and_length
    settings = Settings()
    tube = numpy.asarray(tube, dtype=numpy.float64)
    _finite_text(tube, "the tube")
    settings.array("AREAFILE", numpy.concatenate([[NODE], tube]))
    own = {"BRANCHES", "ASPIRATION", "FRICATION", "FLUX", "FOLDTURB", "CLOSING", "SKIN", "AXIAL", "WALLVAR",
           "STATICAREA", "VISCOUS", "HOPLOSS", "AREAFILE", "RAW", "RAWDOUBLE", "ACTIVELENGTH", "TRACE"}
    if keyframes is not None:
        own.add("AREASEQ")
    if jet > 0.0:
        own.add("REVJET")
        if os.environ.get("JETBOUND", "0") == "1":
            own.add("JETBOUND")
        if jet_gate is not None:
            own.add("JETSEQ")
        if seal > 0.0:
            own.add("REVJETCLOSE")
    if seal > 0.0:
        own.add("SEAL")
        if os.environ.get("SEALBLEND", "0") == "1":
            own.add("SEALBLEND")
    if velum is not None:
        own.update(("VELUMSEQ", "REVNOSE", "NOSEHOP", "NOSESHELF"))
        if nose_dump is not None:
            own.add("NOSEDUMP")
        if port_waves is not None:
            own.add("PORTWAVES")
    if folds is not None:
        own.add("FOLDKEYS" if len(folds[1]) > 0 and len(folds[1][0]) == 3 else "FOLDSEQ")
    field_areas, field_velum, field_folds = field_parts(field)
    if field_areas is not None:
        own.add("AREATICKS")
    if field_velum is not None:
        own.update(("VELUMTICKS", "REVNOSE", "NOSEHOP", "NOSESHELF"))
    if field_folds is not None:
        own.add("FOLDTICKS")
    if forced is not None:
        own.add("FORCEWAVE")
        if observed is not None and observed[1] > 0.0:
            own.update(("FORCEFLOW", "OBSERVE"))
    if second is not None:
        own.update(("AREAFILE2", "GLIDE"))
    inherited = {name: os.environ[name] for name in SUSTAIN_VARIABLES if name in os.environ and name not in own}
    for name in FILE_VARIABLES:
        if name in inherited:
            inherited[name] = _beside_root(inherited[name])
    if inherited:
        settings.environment(inherited)
    # The command line: "aa" WAVE "%g" pitch "800" "%g" seconds.
    settings.number("PITCH", float("%g" % pitch))
    settings.number("LUNG", 800.0)
    settings.number("SECONDS", float("%g" % seconds))
    if wave is not None:
        settings.text("WAVE", _beside_root(wave))
    settings.text("BRANCHES", "")
    for name in ("ASPIRATION", "FRICATION", "FLUX", "FOLDTURB", "CLOSING", "SKIN", "AXIAL", "WALLVAR"):
        settings.number(name, 0.0)
    settings.number("STATICAREA", 1.0)
    settings.number("VISCOUS", float("%g" % VISCOUS))
    settings.number("HOPLOSS", float("%g" % HOP))
    settings.number("RAWDOUBLE", 1.0)
    settings.number("ACTIVELENGTH", float(length))
    if keyframes is not None:
        frames = [numpy.asarray(frame, dtype=numpy.float64) for frame in keyframes[1]]
        rows = [numpy.concatenate([[float(time)], frame]) for time, frame in zip(keyframes[0], frames)]
        for row in rows:
            _finite_text(row, "AREASEQ")
            if len(row) - 1 != len(tube):
                raise ValueError("AREASEQ keyframe at %g s has %d areas, not %d" % (row[0], len(row) - 1, len(tube)))
        settings.array("AREASEQ", numpy.array(rows).reshape(len(rows), len(tube) + 1), len(tube) + 1)
    if jet > 0.0:
        settings.number("REVJET", float(jet))
        if os.environ.get("JETBOUND", "0") == "1":
            settings.number("JETBOUND", 1.0)
        if jet_gate is not None:
            rows = _stream(zip(jet_gate[0], jet_gate[1]), 2)
            _finite_text(rows, "JETSEQ")
            settings.array("JETSEQ", rows, 2)
        if seal > 0.0:
            settings.number("REVJETCLOSE", JET_CLOSING)
    if seal > 0.0:
        settings.number("SEAL", float(seal))
        if os.environ.get("SEALBLEND", "0") == "1":
            # The engine's sub-step seal (Tube::set_seal_blend), as the
            # replica's observer.SEAL_BLEND.
            settings.number("SEALBLEND", 1.0)
    if velum is not None:
        rows = _stream(zip(velum[0], velum[1]), 2)
        _finite_text(rows, "VELUMSEQ")
        settings.array("VELUMSEQ", rows, 2)
        settings.number("REVNOSE", 1.0)
        settings.number("NOSEHOP", NOSE_POLE)
        settings.number("NOSESHELF", NOSE_ZERO)
        if nose_dump is not None:
            settings.text("NOSEDUMP", _beside_root(nose_dump))
    if folds is not None:
        smooth = len(folds[1]) > 0 and len(folds[1][0]) == 3
        rows = _stream([[time] + list(posture) for time, posture in zip(folds[0], folds[1])], 4 if smooth else 5)
        _finite_text(rows, "FOLDKEYS" if smooth else "FOLDSEQ")
        settings.array("FOLDKEYS" if smooth else "FOLDSEQ", rows, 4 if smooth else 5)
    if forced is not None:
        settings.array("FORCEWAVE", numpy.asarray(forced, numpy.float64).ravel())
        if observed is not None and observed[1] > 0.0:
            settings.array("FORCEFLOW", numpy.asarray(observed[0], numpy.float64).ravel())
            settings.number("OBSERVE", float(observed[1]))
    if field_areas is not None:
        _finite_text(field_areas, "AREATICKS")
        if field_areas.shape[1] != len(tube):
            raise ValueError("AREATICKS rows have %d areas, not %d" % (field_areas.shape[1], len(tube)))
        settings.array("AREATICKS", field_areas, len(tube))
    if field_velum is not None:
        _finite_text(field_velum, "VELUMTICKS")
        settings.array("VELUMTICKS", field_velum, 1)
        settings.number("REVNOSE", 1.0)
        settings.number("NOSEHOP", NOSE_POLE)
        settings.number("NOSESHELF", NOSE_ZERO)
    if field_folds is not None:
        _finite_text(field_folds, "FOLDTICKS")
        settings.array("FOLDTICKS", field_folds, 3)
    if second is not None:
        second = numpy.asarray(second, dtype=numpy.float64)
        _finite_text(second, "AREAFILE2")
        if len(second) != len(tube) or glide[1] <= 0.0:
            raise ValueError("AREAFILE2 needs %d areas and GLIDE=start,duration" % len(tube))
        settings.array("AREAFILE2", numpy.concatenate([[NODE], second]))
        settings.pair("GLIDE", float(glide[0]), float(glide[1]))
    return settings


def render(tube_and_length, seconds=0.3, pitch=120.0, second=None, glide=None, forced=None, keyframes=None,
           wave=None, jet=0.0, seal=0.0, velum=None, nose_dump=None, port_waves=None, folds=None, observed=None,
           jet_gate=None, field=None):
    """posture_recover.render, in-process: the flow entering and the flow
    radiated, per step (see posture_recover.render for every argument)."""
    settings = prepare(tube_and_length, seconds, pitch, second, glide, forced, keyframes, wave, jet, seal, velum,
                       nose_dump, port_waves, folds, observed, jet_gate, field)
    steps = settings.steps()
    pairs = numpy.empty((steps, 2))
    waves = numpy.empty((steps, 6)) if (velum is not None and port_waves is not None) else None
    status = ctypes.c_int32(0)
    written = settings.lib.tract_render(settings.handle, _pointer(pairs), _pointer(waves), None, ctypes.byref(status))
    if written < 0:
        raise ValueError(settings.lib.tract_error().decode())
    if waves is not None:
        waves[:written].tofile(_beside_root(port_waves))
    pairs = pairs[:written]
    return pairs[:, 0], pairs[:, 1]


class Incremental:
    """posture_recover.render, re-running only what changed since the last
    call. Each render keeps snapshots of its run every `every` steps (10 ms
    by default) and its outputs; the next render finds the first step whose
    inputs differ (Settings.divergence: the body's construction compared
    whole, then every control tick, fold trigger, forced wave and observed
    flow), resumes from the latest snapshot at or before it, and takes the
    earlier outputs as they were. Exact: the state at a step depends only on
    the construction and the inputs before it (bindings/equivalence.py
    checks the results byte for byte). A search that moves a later phone
    re-runs only from that phone's first tick. One per process (each
    worker its own)."""

    def __init__(self, every=960):
        self.every = int(every)
        self.settings = None
        self.snapshots = []   # (step, Run), ascending
        self.raw = None
        self.waves = None
        self.reached = 0      # the steps the last render wrote
        self.reused = 0       # steps taken from the last render, this call

    def render(self, tube_and_length, seconds=0.3, pitch=120.0, second=None, glide=None, forced=None, keyframes=None,
               wave=None, jet=0.0, seal=0.0, velum=None, nose_dump=None, port_waves=None, folds=None, observed=None,
               jet_gate=None):
        settings = prepare(tube_and_length, seconds, pitch, second, glide, forced, keyframes, wave, jet, seal, velum,
                           nose_dump, port_waves, folds, observed, jet_gate)
        steps = settings.steps()
        nasal = velum is not None
        resume = -1
        if self.settings is not None and self.settings.lib is settings.lib:
            resume = min(self.settings.divergence(settings), self.reached)
        snapshots = [(step, run) for step, run in self.snapshots if step <= resume]
        raw = numpy.empty((steps, 2))
        waves = numpy.empty((steps, 6)) if nasal else None
        if snapshots:
            start, saved = snapshots[-1]
            run = saved.snapshot()
            run.settings = settings
            raw[:start] = self.raw[:start]
            if nasal:
                waves[:start] = self.waves[:start]
        else:
            start = 0
            run = Run(settings)
            snapshots = [(0, run.snapshot())]
        self.reused = start
        reached = start
        finite = True
        while finite and reached < steps:
            end = min((reached // self.every + 1) * self.every, steps)
            before = reached
            reached, finite = run.advance_into(end, raw, waves)
            if reached == before:
                raise RuntimeError("the engine did not advance at step %d" % reached)
            if finite and reached == end and reached < steps:
                snapshots.append((reached, run.snapshot()))
        if finite:
            run.finish()
        self.settings = settings
        self.snapshots = snapshots
        self.raw = raw
        self.waves = waves
        self.reached = reached
        if nasal and port_waves is not None:
            waves[:reached].tofile(_beside_root(port_waves))
        pairs = raw[:reached].copy()
        return pairs[:, 0], pairs[:, 1]
