"""The inverse kernels' benchmark and equivalence harness.

Reproduces the kernel calls of sequence_recover's residuals for five words
(sealed_residual's observer path, march_residual's split march, the plain
peel with the gated jet, a nasal word with the jet, and a sealed word with
the jet), times each kernel, and compares every kernel output of the fast
tree with the frozen reference's, bit for bit (numpy.array_equal, NaN equal
to NaN).

usage (on the Mini):
  bench_kernels.py ab OUT_PREFIX [--repeat R] [--cases a,b,c,d,e] [--variants truth,time,node]
      both trees in one process (the reference's modules under ref_* names),
      evaluated alternately; every output compared in process (with --save,
      also written to OUT_PREFIX_{ref,fast}.npz)
  bench_kernels.py run OUT.npz [--repeat R] [--cases ...] [--variants ...]
      the fast tree alone (profile)
  bench_kernels.py resume - [--cases ...] [--variants time,node]
      the resume gate: each perturbed hypothesis evaluated whole and resumed
      from the truth's record (geometry, peel / observer / march), compared
  bench_kernels.py compare FIRST.npz SECOND.npz

The engine's renders (our folds) are memoized on disk under /tmp/esf (keyed
by every input's bytes), so both trees see the same engine output for the
same input and the kernels' time is measured apart from the engine's. The
helpers that build a residual's inputs (keyframes, plan, velum_keyframes,
jet and fold keyframes) are sequence_recover's own, from the fast tree; the
kernels are each tree's. numba's caches go to /tmp/esf/numba_cache (never
into the reference tree)."""
import os

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/esf/numba_cache")

import hashlib
import importlib.util
import sys
import tempfile
import time

import numpy

FAST_DIR = os.path.dirname(os.path.abspath(__file__))
REFERENCE_DIR = os.path.expanduser("~/ducklips/tract/analysis")
SCRATCH = "/tmp/esf"
CACHE = os.path.join(SCRATCH, "render_cache")
KERNEL_MODULES = ("reverse_core", "kernels_jet", "marching", "observer")
SAVE = [False]   # --save: write OUT_PREFIX_{ref,fast}.npz

WORDS = {
    "a": "AA,B@lips,AA",
    "b": "AA,M@lips^0.47,AA",
    "c": "AA,S@alveolar~0.08!,AA",
    "d": "IY,UW,N@alveolar^0.5,IH,F@lips~0.15!,AO,ER,M@lips^0.5",
    "e": "G@velar,AA,L@alveolar=0.25,F@lips~0.15!",
}


class Timer:
    """Accumulated wall and CPU seconds per label."""

    def __init__(self) -> None:
        self.wall: dict = {}
        self.cpu: dict = {}
        self.calls: dict = {}

    def add(self, label: str, wall: float, cpu: float) -> None:
        self.wall[label] = self.wall.get(label, 0.0) + wall
        self.cpu[label] = self.cpu.get(label, 0.0) + cpu
        self.calls[label] = self.calls.get(label, 0) + 1

    def reset(self) -> None:
        self.wall.clear()
        self.cpu.clear()
        self.calls.clear()


TIMER = Timer()


class Timed:
    """A callable that times another (a named class, not a closure)."""

    def __init__(self, label: str, function) -> None:
        self.label = label
        self.function = function

    def __call__(self, *arguments, **keywords):
        wall = time.perf_counter()
        cpu = time.process_time()
        try:
            return self.function(*arguments, **keywords)
        finally:
            TIMER.add(self.label, time.perf_counter() - wall, time.process_time() - cpu)


class ScratchTemp:
    """tempfile's TemporaryDirectory, always under /tmp/esf (the reference's
    /tmp/es is never touched)."""

    @staticmethod
    def TemporaryDirectory(dir: str = None):
        return tempfile.TemporaryDirectory(dir=SCRATCH)


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def digest(parts: list) -> str:
    hasher = hashlib.sha1()
    for part in parts:
        if isinstance(part, numpy.ndarray):
            hasher.update(str(part.dtype).encode())
            hasher.update(str(part.shape).encode())
            hasher.update(numpy.ascontiguousarray(part).tobytes())
        elif isinstance(part, (list, tuple)):
            hasher.update(b"(")
            hasher.update(digest(list(part)).encode())
            hasher.update(b")")
        else:
            hasher.update(repr(part).encode())
    return hasher.hexdigest()


class Kernels:
    """One tree's kernel modules (reverse_core, marching, observer)."""

    def __init__(self, name: str, modules: dict) -> None:
        self.name = name
        self.reverse_core = modules["reverse_core"]
        self.marching = modules["marching"]
        self.observer = modules["observer"]


def load_fast() -> Kernels:
    """The fast tree's kernels under their own names (what posture_recover,
    tube_recover and sequence_recover import)."""
    if FAST_DIR not in sys.path:
        sys.path.insert(0, FAST_DIR)
    import reverse_core
    import marching
    import observer
    return Kernels("fast", {"reverse_core": reverse_core, "marching": marching, "observer": observer})


def load_reference(directory: str) -> Kernels:
    """The reference tree's kernels as ref_reverse_core, ref_marching,
    ref_observer: while each loads, its siblings' plain names point at the
    reference's modules (their `import marching` finds them), then the fast
    tree's are put back."""
    saved = {name: sys.modules.get(name) for name in KERNEL_MODULES}
    path_before = list(sys.path)
    loaded = {}
    try:
        for name in ("reverse_core", "marching", "observer"):
            spec = importlib.util.spec_from_file_location("ref_" + name, os.path.join(directory, name + ".py"))
            module = importlib.util.module_from_spec(spec)
            sys.modules["ref_" + name] = module
            sys.modules[name] = module
            spec.loader.exec_module(module)
            loaded[name] = module
    finally:
        for name, module in saved.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module
        sys.path[:] = path_before
    return Kernels("ref", loaded)


class Harness:
    """The fast tree's helpers with a memoized render, and the kernel sets."""

    def __init__(self, with_reference: bool) -> None:
        self.fast = load_fast()
        self.posture_recover = load_module("posture_recover", os.path.join(FAST_DIR, "posture_recover.py"))
        self.posture_recover.tempfile = ScratchTemp
        self.tube_recover = load_module("tube_recover", os.path.join(FAST_DIR, "tube_recover.py"))
        self.sequence_recover = load_module("sequence_recover", os.path.join(FAST_DIR, "sequence_recover.py"))
        self.reference = load_reference(REFERENCE_DIR) if with_reference else None
        self.original_render = self.posture_recover.render
        self.posture_recover.render = self.render
        self.renders = 0
        os.makedirs(CACHE, exist_ok=True)

    @staticmethod
    def first_difference(first, second) -> int:
        import reverse_core
        return reverse_core.first_difference(first, second)

    def sets(self) -> list:
        return [self.fast] if self.reference is None else [self.reference, self.fast]

    def render(self, tube_and_length, seconds=0.3, pitch=120.0, second=None, glide=None, forced=None, keyframes=None,
               wave=None, jet=0.0, seal=0.0, velum=None, nose_dump=None, port_waves=None, folds=None, observed=None,
               jet_gate=None):
        """posture_recover.render memoized (no wave or nose dumps)."""
        key = digest([numpy.asarray(tube_and_length[0], float), float(tube_and_length[1]), seconds, pitch, second,
                      glide, None if forced is None else numpy.asarray(forced, float), keyframes, jet, seal, velum,
                      port_waves is not None, folds, None if observed is None else
                      (numpy.asarray(observed[0], float), observed[1]), jet_gate])
        path = os.path.join(CACHE, key + ".npz")
        if os.path.exists(path) and wave is None and nose_dump is None:
            stored = numpy.load(path)
            first, second_column = stored["first"], stored["second"]
            if port_waves is not None:
                stored["waves"].tofile(port_waves)
            return first, second_column
        first, second_column = self.original_render(tube_and_length, seconds, pitch, second, glide, forced, keyframes,
                                                    wave, jet, seal, velum, nose_dump, port_waves, folds, observed,
                                                    jet_gate)
        self.renders += 1
        if wave is None and nose_dump is None:
            waves = numpy.fromfile(port_waves, dtype=numpy.float64) if port_waves is not None else numpy.zeros(0)
            numpy.savez(path, first=first, second=second_column, waves=waves)
        return first, second_column


INSTRUMENTED = {
    "reverse_core": ("keyframe_schedule", "moving_geometry", "velum_port", "jet_gate_steps", "peel_moving_glottis"),
    "observer": ("implied", "closures", "_peel_states", "_forward_steps", "_forward_branched_steps",
                 "_nose_steps_isolated"),
    "marching": ("march_windows", "_march_run", "_peel_below", "jet_table"),
}


def instrument(kernels: Kernels) -> None:
    """Time every kernel the residuals reach (module attributes replaced;
    the kernels' own module-global lookups find the timed ones)."""
    for module_name, labels in INSTRUMENTED.items():
        module = getattr(kernels, module_name)
        for label in labels:
            if hasattr(module, label):
                setattr(module, label, Timed("%s:%s.%s" % (kernels.name, module_name, label), getattr(module, label)))


def word_setup(harness: Harness, spec: str) -> dict:
    """sequence_recover.main's parsing and known answer, without its search."""
    sr = harness.sequence_recover
    pr = harness.posture_recover
    profile = sr.table(os.path.join(sr.DATA, "profile_bdl5.organs"))
    specs = spec.split(",")
    names = [item.split("@")[0] for item in specs]
    phones = len(names)
    first, length = pr.nodes(profile[[item for item in specs if "@" not in item][0]])
    count = len(first)
    places = {"lips": count - 2, "alveolar": count - 6, "velar": count - 18}
    stops = {}
    narrowings = {}
    fricatives = []
    voiceless = {}
    true_velum = numpy.zeros(phones)
    for index, item in enumerate(specs):
        if "@" in item:
            place = item.split("@")[1]
            unvoiced = place.endswith("!")
            place = place.rstrip("!")
            if "^" in place:
                place, area = place.split("^")
                true_velum[index] = float(area)
            if "~" in place:
                place, area = place.split("~")
                narrowings[index] = float(area)
                fricatives.append(index)
            if "=" in place:
                place, area = place.split("=")
                narrowings[index] = float(area)
            stops[index] = places[place]
            if unvoiced:
                voiceless[index] = sr.FRICATIVE_FOLDS if index in narrowings else sr.PLOSIVE_FOLDS
    nasals = [k for k in range(phones) if true_velum[k] > 0.0]
    seal = sr.SEAL if len(stops) > len(narrowings) else 0.0
    jet = sr.JET_GAIN if fricatives else 0.0
    fricatives = tuple(fricatives)
    true_tubes = numpy.zeros((phones, count))
    for index, name in enumerate(names):
        if index not in stops:
            true_tubes[index] = pr.nodes(profile[name], count)[0]
    for index in sorted(stops):
        true_tubes[index] = true_tubes[sr.vowel_before(index, stops)]
        true_tubes[index, stops[index]] = narrowings.get(index, sr.STOP_AREA)
    arrivals, leavings = sr.plan(phones)
    seconds = leavings[-1] + 0.05
    times, frames = sr.keyframes(true_tubes, arrivals, leavings)
    port = sr.velum_keyframes(arrivals, leavings, true_velum) if nasals else None
    flow_in, flow_out = harness.render((true_tubes[0], length), seconds, keyframes=(times, frames), seal=seal,
                                       velum=port, jet=jet,
                                       jet_gate=sr.jet_keyframes(arrivals, leavings, fricatives),
                                       folds=sr.fold_keyframes(arrivals, leavings, voiceless))
    return dict(spec=spec, phones=phones, count=count, length=length, stops=stops, narrowings=narrowings,
                fricatives=fricatives, voiceless=voiceless, velum=true_velum, nasals=nasals, seal=seal, jet=jet,
                tubes=true_tubes, arrivals=arrivals, leavings=leavings, seconds=seconds, flow_out=flow_out)


def make_target(harness: Harness, word: dict):
    tr = harness.tube_recover
    start, size = harness.sequence_recover.window(0.02, word["seconds"] - 0.01)
    target = tr.Target(word["flow_out"], start, size, 0.0, 1.0, True)
    target.seal = word["seal"]
    target.jet = word["jet"]
    target.voiceless = word["voiceless"]
    target.fricatives = word["fricatives"]
    return target


def nose_constants(harness: Harness):
    """sequence_recover.nose_constants with its dump under /tmp/esf."""
    sr = harness.sequence_recover
    if sr.NOSE is None:
        path = os.path.join(SCRATCH, "nose_constants_%d.txt" % os.getpid())
        harness.original_render((numpy.full(46, 2.0), 46.5 * harness.posture_recover.NODE), 0.01,
                                velum=([0.0, 1.0], [0.0, 0.0]), nose_dump=path)
        sr.NOSE = harness.fast.reverse_core.read_nose(path)
        os.remove(path)
    return sr.NOSE


def gate_steps(kernels: Kernels, harness: Harness, arrivals, leavings, fricatives, count):
    """sequence_recover.jet_gate_for with this tree's jet_gate_steps."""
    frames = harness.sequence_recover.jet_keyframes(arrivals, leavings, fricatives)
    if frames is None:
        return None
    return kernels.reverse_core.jet_gate_steps(frames[0], frames[1], count)


def geometry_digest_untimed(geometry) -> numpy.ndarray:
    """The geometry as sha1 bytes (its 170+ MB not stored), plus its shape."""
    array = numpy.ascontiguousarray(numpy.asarray(geometry))
    value = hashlib.sha1(array.tobytes()).digest()
    return numpy.frombuffer(value + str(array.shape).encode(), numpy.uint8).copy()


geometry_digest = Timed("bench.digest (not a kernel)", geometry_digest_untimed)


def sealed_case(kernels: Kernels, harness: Harness, tubes, arrivals, leavings, length, target,
                outputs: dict) -> numpy.ndarray:
    """sequence_recover.sealed_residual, its kernel outputs kept."""
    sr = harness.sequence_recover
    rc = kernels.reverse_core
    pr = harness.posture_recover
    times, frames = sr.keyframes(tubes, arrivals, leavings)
    schedule = rc.keyframe_schedule(times, frames, len(target.flow_out))
    geometry = rc.moving_geometry(schedule, pr.VISCOUS)
    free = Timed("engine.free_flow", sr.FreeFlow(frames[0], length, (len(target.flow_out) + 1) / sr.STEP_RATE, times,
                                                 frames, target.seal,
                                                 folds=sr.fold_keyframes(arrivals, leavings, target.voiceless),
                                                 jet=target.jet,
                                                 jet_gate=sr.jet_keyframes(arrivals, leavings, target.fricatives)))
    gate = gate_steps(kernels, harness, arrivals, leavings, target.fricatives, len(target.flow_out))
    implied, pressure, revealed, modelled = kernels.observer.implied(
        geometry, target.flow_out, target.seal, pr.HOP, length - frames.shape[1] * pr.NODE, free, hide=sr.HIDE_AREA,
        jet=target.jet, jet_closing=pr.JET_CLOSING, jet_gate=gate)
    forced = 0.5 * (pressure - implied * rc.RHO_C / geometry[:, 0, 0])
    forced[~revealed] = numpy.nan
    driven = free(forced, None)
    size = min(len(driven), len(implied))
    difference = driven[:size] - implied[:size]
    hidden = ~revealed[:size]
    difference[hidden] = modelled[:size][hidden] - target.flow_out[:size][hidden]
    outputs.update(schedule=schedule, geometry_digest=geometry_digest(geometry), implied=implied, pressure=pressure,
                   revealed=revealed, modelled=modelled)
    return harness.tube_recover.fixed_clip(difference, target)


def plain_case(kernels: Kernels, harness: Harness, tubes, arrivals, leavings, length, target,
               outputs: dict) -> numpy.ndarray:
    """sequence_recover.sequence_residual without a seal (the gated jet)."""
    sr = harness.sequence_recover
    rc = kernels.reverse_core
    pr = harness.posture_recover
    times, frames = sr.keyframes(tubes, arrivals, leavings)
    schedule = rc.keyframe_schedule(times, frames, len(target.flow_out))
    geometry = rc.moving_geometry(schedule, pr.VISCOUS)
    gate = gate_steps(kernels, harness, arrivals, leavings, target.fricatives, len(target.flow_out))
    implied, pressure = rc.peel_moving_glottis(geometry, target.flow_out, hop=pr.HOP,
                                               end_extra=length - frames.shape[1] * pr.NODE, jet=target.jet, gate=gate)
    forced = 0.5 * (pressure - implied * rc.RHO_C / geometry[:, 0, 0])
    driven, _ = Timed("engine.render", pr.render)((frames[0], length), target.seconds, keyframes=(times, frames),
                                                   forced=forced, jet=target.jet,
                                                   folds=sr.fold_keyframes(arrivals, leavings, target.voiceless),
                                                   jet_gate=sr.jet_keyframes(arrivals, leavings, target.fricatives))
    size = min(len(driven), len(implied))
    difference = driven[:size] - implied[:size]
    outputs.update(schedule=schedule, geometry_digest=geometry_digest(geometry), implied=implied, pressure=pressure)
    return harness.tube_recover.fixed_clip(difference, target)


def march_windows_of(harness: Harness, geometry, port, frames, port_node: int, m: int, count: int, seal: float):
    """sequence_recover.march_residual's folds-decided windows."""
    sr = harness.sequence_recover
    windows = []
    oral_lag = frames.shape[1] - 1 - port_node
    front = numpy.zeros_like(geometry[:, :, 0])
    front[:, port_node + 1:] = geometry[:, port_node + 1:, 0]
    front[:, :port_node + 1] = numpy.inf
    narrow = (front.min(axis=1) <= sr.HIDE_AREA) | ((front.min(axis=1) <= sr.WEAK_ORAL) & (port >= 1e-6))
    shut = front.min(axis=1) <= seal
    t = 0
    while t < count:
        if not narrow[t]:
            t += 1
            continue
        start = t
        while t < count and narrow[t]:
            t += 1
        stretch = numpy.nonzero(shut[start:t])[0]
        begin = max(start - oral_lag + 1 - sr.RELEASE_MARGIN, 1)
        if len(stretch) == 0 or numpy.any(port[start + stretch] < 1e-6):
            windows.append((begin, min(t + 2 + sr.RELEASE_MARGIN, count)))
            continue
        first_shut = start + int(stretch[0])
        last_shut = start + int(stretch[-1])
        closing_end = min(first_shut - oral_lag + 1 + sr.RELEASE_MARGIN, count)
        if closing_end > begin:
            windows.append((begin, closing_end))
        windows.append((max(last_shut - m + 2 - sr.RELEASE_MARGIN, 1), min(t + 2 + sr.RELEASE_MARGIN, count)))
    merged = []
    for begin, end in sorted(windows):
        if merged and begin <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((begin, end))
    return merged


class MarchPieces:
    """march_residual's callbacks as named objects: the peel below the port
    (returning) and the engine's port waves (render_waves, render_window)."""

    def __init__(self, kernels: Kernels, harness: Harness, geometry, port_node, port_y, area0, jet_steps, frames,
                 length, seconds, times, target, velum, arrivals, leavings) -> None:
        self.kernels = kernels
        self.harness = harness
        self.geometry = geometry
        self.port_node = port_node
        self.port_y = port_y
        self.area0 = area0
        self.jet_steps = jet_steps
        self.frames = frames
        self.length = length
        self.seconds = seconds
        self.times = times
        self.target = target
        self.velum = velum
        self.arrivals = arrivals
        self.leavings = leavings
        self.window_forced = []

    def returning(self, arrays):
        pr = self.harness.posture_recover
        flow, pressure = self.kernels.marching._peel_below(self.geometry, self.port_node, arrays["port_pressure"],
                                                           arrays["from_front"], arrays["returned"], self.port_y,
                                                           pr.HOP, True, self.jet_steps)
        return flow, pressure, 0.5 * (pressure - flow * self.kernels.reverse_core.RHO_C / self.area0)

    def render_waves(self, forced):
        sr = self.harness.sequence_recover
        pr = self.harness.posture_recover
        path = os.path.join(SCRATCH, "waves_%d_%d.raw" % (os.getpid(), numpy.random.randint(1 << 30)))
        flow, _ = Timed("engine.render_waves", pr.render)(
            (self.frames[0], self.length), self.seconds, keyframes=(self.times, self.frames), forced=forced,
            seal=self.target.seal, velum=self.velum, port_waves=path, jet=self.target.jet,
            jet_gate=sr.jet_keyframes(self.arrivals, self.leavings, self.target.fricatives),
            folds=sr.fold_keyframes(self.arrivals, self.leavings, self.target.voiceless))
        waves = numpy.fromfile(path, dtype=numpy.float64).reshape(-1, 6)
        os.remove(path)
        return flow, waves

    def render_window(self, arrays, begin):
        flow, pressure, wave = self.returning(arrays)
        forced = wave.copy()
        forced[max(begin - self.port_node - 2, 0):] = numpy.nan
        self.window_forced.append(forced)
        return self.render_waves(forced)[1]


def march_case(kernels: Kernels, harness: Harness, tubes, arrivals, leavings, velum_areas, length, target,
               outputs: dict) -> numpy.ndarray:
    """sequence_recover.march_residual (its temporary files under /tmp/esf)."""
    sr = harness.sequence_recover
    rc = kernels.reverse_core
    pr = harness.posture_recover
    marching = kernels.marching
    times, frames = sr.keyframes(tubes, arrivals, leavings)
    count = len(target.flow_out)
    schedule = rc.keyframe_schedule(times, frames, count)
    geometry = numpy.ascontiguousarray(rc.moving_geometry(schedule, pr.VISCOUS))
    velum = sr.velum_keyframes(arrivals, leavings, velum_areas)
    port = rc.velum_port(velum[0], velum[1], count)
    nose = nose_constants(harness)
    port_node = nose[0]
    m = len(nose[1])
    end_extra = length - frames.shape[1] * pr.NODE
    if sr.HIDE_AREA > target.seal:
        windows = march_windows_of(harness, geometry, port, frames, port_node, m, count, target.seal)
    else:
        windows = []
        for first, last, node in kernels.observer.closures(geometry, target.seal):
            if last + 1 < count and port[last + 1] >= 1e-6:
                windows.append((max(last - m + 2 - sr.RELEASE_MARGIN, 1), min(last + 2 + sr.RELEASE_MARGIN, count)))
    seconds = (count + 1) / sr.STEP_RATE
    area0 = geometry[:, 0, 0]
    port_y = port / rc.RHO_C
    gate = gate_steps(kernels, harness, arrivals, leavings, target.fricatives, count)
    jet = (target.jet, rc.JET_SEED, pr.JET_CLOSING if target.seal > 0.0 else 0.0, gate)
    jet_steps = marching.jet_table(numpy.ascontiguousarray(geometry[:, :, 0]), jet, gate)
    pieces = MarchPieces(kernels, harness, geometry, port_node, port_y, area0, jet_steps, frames, length, seconds,
                         times, target, velum, arrivals, leavings)
    free = numpy.zeros(count, bool)
    for begin, end in windows:
        free[max(begin - port_node - 2, 0):end] = True
    share, returned, port_pressure, from_front, mode, waves, in_window = marching.march_windows(
        target.flow_out, geometry, port, nose, pr.HOP, end_extra, target.seal, windows, pieces.render_window,
        pr.NOSE_POLE, pr.NOSE_ZERO, jet)
    arrays = {"port_pressure": port_pressure, "from_front": from_front, "returned": returned}
    implied, pressure, wave = pieces.returning(arrays)
    forced = wave.copy()
    forced[free] = numpy.nan
    driven, waves_b = pieces.render_waves(forced)
    size = min(len(driven), len(implied), len(waves_b))
    difference = driven[:size] - implied[:size]
    hidden = free[:size]
    difference[hidden] = waves_b[:size, 4][hidden] + waves_b[:size, 5][hidden] - target.flow_out[:size][hidden]
    outputs.update(schedule=schedule, geometry_digest=geometry_digest(geometry), port=port, share=share,
                   returned=returned, port_pressure=port_pressure, from_front=from_front, mode=mode,
                   implied=implied, pressure=pressure, windows=numpy.array(windows, numpy.int64).reshape(-1, 2))
    for index, value in enumerate(pieces.window_forced):
        outputs["window_forced_%d" % index] = value
    return harness.tube_recover.fixed_clip(difference, target)


def hypothesis(word: dict, variant: str):
    """The truth, or it with one time +0.3 ms (the consonant's leaving), or
    one node of the last phone's tube +1e-3 in log area."""
    tubes = word["tubes"].copy()
    arrivals = word["arrivals"].copy()
    leavings = word["leavings"].copy()
    if variant == "time":
        index = min(word["stops"]) if word["stops"] else 0
        leavings[index] += 3e-4
    elif variant == "node":
        tubes[-1, word["count"] // 2] *= numpy.exp(1e-3)
    return tubes, arrivals, leavings


def evaluate(kernels: Kernels, harness: Harness, word: dict, variant: str, outputs: dict) -> numpy.ndarray:
    target = make_target(harness, word)
    tubes, arrivals, leavings = hypothesis(word, variant)
    if word["nasals"]:
        return march_case(kernels, harness, tubes, arrivals, leavings, word["velum"], word["length"], target, outputs)
    if word["seal"] > 0.0:
        return sealed_case(kernels, harness, tubes, arrivals, leavings, word["length"], target, outputs)
    return plain_case(kernels, harness, tubes, arrivals, leavings, word["length"], target, outputs)


def kernel_seconds(labels: dict, prefix: str) -> float:
    """The kernels' own seconds: the outermost timed kernels of a set (the
    residual's geometry, peels, observer, march; not the engine)."""
    outer = ("reverse_core.keyframe_schedule", "reverse_core.moving_geometry", "reverse_core.velum_port",
             "reverse_core.jet_gate_steps", "reverse_core.peel_moving_glottis", "observer.implied",
             "marching.march_windows", "marching.jet_table")
    total = 0.0
    for name in outer:
        total += labels.get(prefix + ":" + name, 0.0)
    # The peel below the port outside march_windows (the final returning).
    return total


def run(out_prefix: str, repeat: int, cases: list, variants: list, with_reference: bool) -> int:
    harness = Harness(with_reference)
    sets = harness.sets()
    for kernels in sets:
        instrument(kernels)
    saved = {kernels.name: {} for kernels in sets}
    failures = 0
    for case in cases:
        word = word_setup(harness, WORDS[case])
        for variant in variants:
            key = "%s_%s" % (case, variant)
            outputs = {}
            residuals = {}
            for kernels in sets:
                # Warm-up (numba's cache, the engine's memo), kept as the output.
                outputs[kernels.name] = {}
                residuals[kernels.name] = evaluate(kernels, harness, word, variant, outputs[kernels.name])
            TIMER.reset()
            totals = {kernels.name: [] for kernels in sets}
            for index in range(repeat):
                order = sets if index % 2 == 0 else list(reversed(sets))
                for kernels in order:
                    wall = time.perf_counter()
                    cpu = time.process_time()
                    evaluate(kernels, harness, word, variant, {})
                    totals[kernels.name].append((time.perf_counter() - wall, time.process_time() - cpu))
            lines = ["%s %s (%s, %d steps):" % (case, variant, word["spec"], len(word["flow_out"]))]
            for kernels in sets:
                walls = [pair[0] for pair in totals[kernels.name]]
                cpus = [pair[1] for pair in totals[kernels.name]]
                lines.append("  %-4s residual wall median %.3f min %.3f, cpu median %.3f; kernels (outer) cpu %.4f"
                             " wall %.4f per eval" % (
                                 kernels.name, float(numpy.median(walls)), min(walls), float(numpy.median(cpus)),
                                 kernel_seconds(TIMER.cpu, kernels.name) / repeat,
                                 kernel_seconds(TIMER.wall, kernels.name) / repeat))
            for label in sorted(TIMER.wall):
                lines.append("    %-48s wall %.4f s  cpu %.4f s  (%d calls/eval)" % (
                    label, TIMER.wall[label] / repeat, TIMER.cpu[label] / repeat, TIMER.calls[label] // repeat))
            for kernels in sets:
                store = saved[kernels.name]
                store[key + "__residual"] = residuals[kernels.name]
                for name, value in outputs[kernels.name].items():
                    store[key + "__" + name] = numpy.asarray(value)
            if with_reference:
                first = saved["ref"]
                second = saved["fast"]
                differing = [name for name in first if name.startswith(key + "__") and
                             not same_array(first[name], second.get(name))]
                failures += len(differing)
                lines.append("  equivalence: %d arrays, %s" % (
                    len([name for name in first if name.startswith(key + "__")]),
                    "all bit-identical" if not differing else "DIFFER: " + " ".join(differing)))
            print("\n".join(lines))
            sys.stdout.flush()
    if SAVE[0]:
        # (Each schedule is 14-37 MB: kept only when asked, the Mini's disk is shared.)
        for name, store in saved.items():
            numpy.savez("%s_%s.npz" % (out_prefix, name), **store)
    print("engine renders actually run: %d; equivalence failures: %d" % (harness.renders, failures))
    return failures


class CountingCall:
    """A callback that counts its calls (the engine's, inside the kernels)."""

    def __init__(self, function) -> None:
        self.function = function
        self.calls = 0

    def __call__(self, *arguments):
        self.calls += 1
        return self.function(*arguments)


def keyframe_change(first, second) -> float:
    """The time from which two keyframe lists' interpolation can differ:
    the keyframe before their first difference (0 if the first)."""
    if first is None and second is None:
        return numpy.inf
    if first is None or second is None:
        return 0.0
    times_a, values_a = first
    times_b, values_b = second
    size = min(len(times_a), len(times_b))
    for index in range(size):
        if times_a[index] != times_b[index] or not numpy.array_equal(numpy.asarray(values_a[index], float),
                                                                     numpy.asarray(values_b[index], float)):
            return float(times_a[index - 1]) if index > 0 else 0.0
    if len(times_a) != len(times_b):
        return float(times_a[size - 1]) if size > 0 else 0.0
    return numpy.inf


class Inputs:
    """A hypothesis's per-step inputs (the residual's construction)."""

    def __init__(self, harness: Harness, word: dict, variant: str) -> None:
        sr = harness.sequence_recover
        rc = harness.fast.reverse_core
        pr = harness.posture_recover
        self.tubes, self.arrivals, self.leavings = hypothesis(word, variant)
        self.times, self.frames = sr.keyframes(self.tubes, self.arrivals, self.leavings)
        count = len(word["flow_out"])
        self.count = count
        self.schedule = rc.keyframe_schedule(self.times, self.frames, count)
        self.gate = gate_steps(harness.fast, harness, self.arrivals, self.leavings, word["fricatives"], count)
        self.folds = sr.fold_keyframes(self.arrivals, self.leavings, word["voiceless"])
        self.velum = sr.velum_keyframes(self.arrivals, self.leavings, word["velum"]) if word["nasals"] else None
        self.port = rc.velum_port(self.velum[0], self.velum[1], count) if self.velum is not None else None
        self.end_extra = word["length"] - self.frames.shape[1] * pr.NODE


def first_changed(base: Inputs, other: Inputs) -> int:
    """The first step at which anything the kernels or the engine read can
    differ: the schedule's rows, the gate, the port, the folds' keyframes
    (from the keyframe before their first difference, a tick early)."""
    rc_first = Harness.first_difference
    steps = [rc_first(base.schedule, other.schedule)]
    if base.gate is not None:
        steps.append(rc_first(base.gate, other.gate))
    if base.port is not None:
        steps.append(rc_first(base.port, other.port))
    moment = keyframe_change(base.folds, other.folds)
    if numpy.isfinite(moment):
        steps.append(max(int(moment * 96000.0) // 96 * 96 - 96, 0))
    return int(min(steps))


def timed(function, *arguments):
    wall = time.perf_counter()
    cpu = time.process_time()
    result = function(*arguments)
    return result, time.perf_counter() - wall, time.process_time() - cpu


class ResumeGate:
    """For one word: the truth recorded, each perturbed hypothesis evaluated
    whole and resumed from the truth's record; every output compared."""

    def __init__(self, harness: Harness, word: dict) -> None:
        self.harness = harness
        self.word = word
        self.kernels = harness.fast
        self.rc = harness.fast.reverse_core
        self.pr = harness.posture_recover
        self.sr = harness.sequence_recover

    def free_flow(self, inputs: Inputs):
        sr = self.sr
        word = self.word
        target = make_target(self.harness, word)
        return sr.FreeFlow(inputs.frames[0], word["length"], (len(word["flow_out"]) + 1) / sr.STEP_RATE, inputs.times,
                           inputs.frames, word["seal"], folds=inputs.folds, jet=word["jet"],
                           jet_gate=sr.jet_keyframes(inputs.arrivals, inputs.leavings, word["fricatives"]))

    def plain(self, inputs: Inputs, geometry, record=None, changed=0):
        word = self.word
        if record is None:
            return self.rc.peel_moving_record(geometry, word["flow_out"], self.pr.HOP, True, inputs.end_extra,
                                              word["jet"], self.rc.JET_SEED, inputs.gate)
        return self.rc.peel_moving_resume(record, geometry, word["flow_out"], changed, inputs.gate)

    def sealed(self, inputs: Inputs, geometry, free, record=None, changed=0):
        word = self.word
        observer = self.kernels.observer
        if record is None:
            return observer.implied_record(geometry, word["flow_out"], word["seal"], self.pr.HOP, inputs.end_extra,
                                           free, hide=self.sr.HIDE_AREA, jet=word["jet"],
                                           jet_closing=self.pr.JET_CLOSING, jet_gate=inputs.gate)
        return observer.implied_resume(record, geometry, word["flow_out"], free, changed, inputs.gate)

    def march(self, inputs: Inputs, geometry, render, record=None, changed=0, below=None):
        word = self.word
        sr = self.sr
        marching = self.kernels.marching
        nose = nose_constants(self.harness)
        port_node = nose[0]
        m = len(nose[1])
        count = inputs.count
        windows = march_windows_of(self.harness, geometry, inputs.port, inputs.frames, port_node, m, count,
                                   word["seal"])
        jet = (word["jet"], self.rc.JET_SEED, self.pr.JET_CLOSING if word["seal"] > 0.0 else 0.0, inputs.gate)
        arguments = (word["flow_out"], geometry, inputs.port, nose, self.pr.HOP, inputs.end_extra, word["seal"],
                     windows, render, self.pr.NOSE_POLE, self.pr.NOSE_ZERO, jet)
        if record is None:
            outputs = marching.march_windows_record(*arguments)
        else:
            outputs = marching.march_windows_resume(record, *arguments, first_changed=changed)
        port_y = inputs.port / self.rc.RHO_C
        jet_steps = marching.jet_table(numpy.ascontiguousarray(geometry[:, :, 0]), jet, inputs.gate)
        share, returned, port_pressure, from_front = outputs[:4]
        if below is None:
            run = marching.BelowPeel(geometry, port_node, port_pressure, from_front, returned, port_y, self.pr.HOP,
                                     True, jet_steps, record=True)
        else:
            changed_below = min(changed, Harness.first_difference(port_pressure, below.pressure_input),
                                Harness.first_difference(from_front, below.from_front),
                                Harness.first_difference(returned, below.returned))
            run = below.resumed(port_pressure, from_front, returned, changed_below, geometry, port_y, jet_steps)
        run.pressure_input = numpy.array(port_pressure)
        run.advance(run.count)
        return outputs, run

    def window_renderer(self, inputs: Inputs, geometry):
        """march_residual's render_window for this hypothesis (counted)."""
        word = self.word
        nose = nose_constants(self.harness)
        port_node = nose[0]
        count = inputs.count
        port_y = inputs.port / self.rc.RHO_C
        jet = (word["jet"], self.rc.JET_SEED, self.pr.JET_CLOSING if word["seal"] > 0.0 else 0.0, inputs.gate)
        jet_steps = self.kernels.marching.jet_table(numpy.ascontiguousarray(geometry[:, :, 0]), jet, inputs.gate)
        target = make_target(self.harness, word)
        pieces = MarchPieces(self.kernels, self.harness, geometry, port_node, port_y, geometry[:, 0, 0], jet_steps,
                             inputs.frames, word["length"], (count + 1) / self.sr.STEP_RATE, inputs.times, target,
                             inputs.velum, inputs.arrivals, inputs.leavings)
        return CountingCall(pieces.render_window)


def resume_gate(cases: list, variants: list) -> int:
    harness = Harness(False)
    rc = harness.fast.reverse_core
    pr = harness.posture_recover
    failures = 0
    for case in cases:
        word = word_setup(harness, WORDS[case])
        gate = ResumeGate(harness, word)
        base = Inputs(harness, word, "truth")
        base_geometry, geometry_record = rc.moving_geometry_record(base.schedule, pr.VISCOUS)
        if word["nasals"]:
            render = gate.window_renderer(base, base_geometry)
            base_outputs, base_below = gate.march(base, base_geometry, render)
            base_record = base_outputs[7]
        elif word["seal"] > 0.0:
            free = CountingCall(gate.free_flow(base))
            base_outputs = gate.sealed(base, base_geometry, free)
            base_record = base_outputs[4]
        else:
            base_outputs = gate.plain(base, base_geometry)
            base_record = base_outputs[2]
        for variant in variants:
            other = Inputs(harness, word, variant)
            changed = first_changed(base, other)
            full_geometry, full_geometry_seconds, _ = timed(rc.moving_geometry, other.schedule, pr.VISCOUS)
            (resumed_geometry, _), resumed_geometry_seconds, _ = timed(rc.moving_geometry_resume, other.schedule,
                                                                      geometry_record, base_geometry, changed)
            if variant == variants[0]:
                workspace = rc.GeometryWorkspace(base_geometry, geometry_record)
            workspace_geometry, workspace_seconds, _ = timed(workspace.resume, other.schedule, changed)
            same = [("geometry", same_array(full_geometry, resumed_geometry)),
                    ("geometry(workspace)", same_array(full_geometry, workspace_geometry))]
            calls = ""
            if word["nasals"]:
                render_full = gate.window_renderer(other, full_geometry)
                (full, full_below), full_seconds, full_cpu = timed(gate.march, other, full_geometry, render_full)
                render_resumed = gate.window_renderer(other, resumed_geometry)
                (resumed, resumed_below), resumed_seconds, resumed_cpu = timed(
                    gate.march, other, resumed_geometry, render_resumed, base_record, changed, base_below)
                for index, name in enumerate(("share", "returned", "port_pressure", "from_front", "mode",
                                              "in_window")):
                    position = index if index < 5 else 6
                    same.append((name, same_array(full[position], resumed[position])))
                # The window waves on the rows the march reads (each window's
                # from two steps before it to a nose length past it); the rows
                # between windows are a render's free run, never read (and
                # march_residual does not use them): a reused window keeps
                # the record's there.
                read = numpy.zeros(len(full[5]), bool)
                nose_length = len(nose_constants(harness)[1])
                for begin, end in resumed[7].windows:
                    read[max(begin - 2, 0):min(end + nose_length + 1, len(read))] = True
                same.append(("waves(read rows)", same_array(full[5][read], resumed[5][read])))
                if not same_array(full[5], resumed[5]):
                    calls += "; waves differ on %d unread rows" % int(numpy.sum(numpy.any(full[5] != resumed[5],
                                                                                         axis=1)))
                same.append(("below_flow", same_array(full_below.arriving, resumed_below.arriving)))
                same.append(("below_pressure", same_array(full_below.pressure, resumed_below.pressure)))
                calls = "window renders %d -> %d" % (render_full.calls, render_resumed.calls) + calls
            elif word["seal"] > 0.0:
                free_full = CountingCall(gate.free_flow(other))
                full, full_seconds, full_cpu = timed(gate.sealed, other, full_geometry, free_full)
                free_resumed = CountingCall(gate.free_flow(other))
                resumed, resumed_seconds, resumed_cpu = timed(gate.sealed, other, resumed_geometry, free_resumed,
                                                              base_record, changed)
                for index, name in enumerate(("flow", "pressure", "revealed", "modelled")):
                    same.append((name, same_array(full[index], resumed[index])))
                calls = "free_flow calls %d -> %d" % (free_full.calls, free_resumed.calls)
            else:
                full, full_seconds, full_cpu = timed(gate.plain, other, full_geometry)
                resumed, resumed_seconds, resumed_cpu = timed(gate.plain, other, resumed_geometry, base_record,
                                                              changed)
                same.append(("flow", same_array(full[0], resumed[0])))
                same.append(("pressure", same_array(full[1], resumed[1])))
            bad = [name for name, equal in same if not equal]
            failures += len(bad)
            print("%s %s: first change at step %d of %d (%.0f%% left); geometry %.4f -> %.4f s (workspace %.4f); "
                  "kernel %.4f -> %.4f s wall (cpu %.4f -> %.4f); %s; %s" % (
                      case, variant, changed, base.count, 100.0 * (base.count - changed) / base.count,
                      full_geometry_seconds, resumed_geometry_seconds, workspace_seconds, full_seconds,
                      resumed_seconds, full_cpu,
                      resumed_cpu, calls, "bit-identical (%d arrays)" % len(same) if not bad else "DIFFER: " +
                      " ".join(bad)))
            sys.stdout.flush()
    print("resume gate failures: %d" % failures)
    return failures


def same_array(first, second) -> bool:
    if second is None:
        return False
    a = numpy.asarray(first)
    b = numpy.asarray(second)
    if a.shape != b.shape or a.dtype != b.dtype:
        return False
    return bool(numpy.array_equal(a, b, equal_nan=a.dtype.kind == "f"))


def compare(first_path: str, second_path: str) -> int:
    first = numpy.load(first_path)
    second = numpy.load(second_path)
    failures = 0
    for key in sorted(set(first.files) | set(second.files)):
        if key not in first.files or key not in second.files:
            print("MISSING %s" % key)
            failures += 1
            continue
        a = first[key]
        b = second[key]
        if a.shape != b.shape:
            print("SHAPE %s %s %s" % (key, a.shape, b.shape))
            failures += 1
            continue
        if not same_array(a, b):
            failures += 1
            scale = float(numpy.nanmax(numpy.abs(a))) if a.size and a.dtype.kind == "f" else 1.0
            deviation = float(numpy.nanmax(numpy.abs(a.astype(float) - b.astype(float)))) if a.size else 0.0
            print("DIFF %s max |d| %.3e (relative %.3e)" % (key, deviation, deviation / max(scale, 1e-300)))
    print("%d arrays, %d differ" % (len(first.files), failures))
    return failures


def main() -> None:
    mode = sys.argv[1]
    if mode == "compare":
        sys.exit(1 if compare(sys.argv[2], sys.argv[3]) else 0)
    out = sys.argv[2]
    repeat = 3
    cases = ["a", "b", "c", "d", "e"]
    variants = ["truth", "time", "node"]
    arguments = sys.argv[3:]
    if "--save" in arguments:
        SAVE[0] = True
        arguments = [argument for argument in arguments if argument != "--save"]
    for index in range(0, len(arguments), 2):
        if arguments[index] == "--repeat":
            repeat = int(arguments[index + 1])
        elif arguments[index] == "--cases":
            cases = arguments[index + 1].split(",")
        elif arguments[index] == "--variants":
            variants = arguments[index + 1].split(",")
    if mode == "resume":
        failures = resume_gate(cases, [variant for variant in variants if variant != "truth"])
    else:
        failures = run(out, repeat, cases, variants, mode == "ab")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
