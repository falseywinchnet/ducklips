"""Where a word's search spends its time: evaluations, wall and CPU per
stage of the driver (sequence_recover, tube_recover), and the evaluations
the exact shortcuts saved.

A stage is a path ("phone 2 AA/recover/basin 3/node rung 1"); every
counter is inclusive (a stage's numbers contain its sub-stages'). The
pool's work is timed in the workers (their own CPU plus the engine
subprocesses they waited for: os.times' children); the main process's by
its own os.times over the stage. Wall time is the main process's.

The pool's size comes from the environment (POOL, default 10: the
reference's) and its tasks are handed out one at a time (CHUNK, default
1): with pool.map's default chunks (ceil(tasks / 4 / workers)) a
145-column Jacobian on 4 workers ran as 15 chunks of 10, 40 tasks deep on
the busiest worker against 37 with single tasks."""
import json
import os
import sys
import time
from typing import Any, Callable, Dict, Iterator, List, Optional, Sequence, Tuple
from contextlib import contextmanager

POOL = int(os.environ.get("POOL", "10"))
CHUNK = int(os.environ.get("CHUNK", "1"))


def process_cpu() -> float:
    """This process's CPU seconds, its waited-for children's included (the
    engine runs as a subprocess)."""
    times = os.times()
    return times.user + times.system + times.children_user + times.children_system


class StageRecord:
    """One stage path's inclusive counters."""

    def __init__(self) -> None:
        self.calls: int = 0               # times the stage was entered
        self.evaluations: int = 0         # residuals computed (pool and main)
        self.pool_tasks: int = 0          # of them, in the pool
        self.saved: int = 0               # evaluations an exact shortcut made unnecessary
        self.worker_cpu: float = 0.0      # s, in the pool's workers
        self.main_cpu: float = 0.0        # s, in the main process
        self.wall: float = 0.0            # s

    def as_dict(self) -> Dict[str, float]:
        return {"calls": self.calls, "evaluations": self.evaluations, "pool_tasks": self.pool_tasks,
                "saved": self.saved, "worker_cpu": self.worker_cpu, "main_cpu": self.main_cpu, "wall": self.wall}


class Ledger:
    """The stage tree's counters (one per process; only the main process's
    is reported)."""

    def __init__(self) -> None:
        self.records: Dict[str, StageRecord] = {}
        self.stack: List[str] = []
        self.order: List[str] = []

    def path(self) -> str:
        return "/".join(self.stack) if self.stack else "(top)"

    def record(self, path: str) -> StageRecord:
        if path not in self.records:
            self.records[path] = StageRecord()
            self.order.append(path)
        return self.records[path]

    def open_paths(self) -> List[str]:
        """Every path on the stack (the counters are inclusive)."""
        paths = ["(top)"]
        for depth in range(1, len(self.stack) + 1):
            paths.append("/".join(self.stack[:depth]))
        return paths

    def add(self, evaluations: int = 0, pool_tasks: int = 0, saved: int = 0, worker_cpu: float = 0.0) -> None:
        for path in self.open_paths():
            entry = self.record(path)
            entry.evaluations += evaluations
            entry.pool_tasks += pool_tasks
            entry.saved += saved
            entry.worker_cpu += worker_cpu

    @contextmanager
    def stage(self, name: str) -> Iterator[None]:
        self.stack.append(name)
        path = self.path()
        self.record(path).calls += 1
        wall_start = time.time()
        cpu_start = process_cpu()
        try:
            yield
        finally:
            entry = self.record(path)
            entry.wall += time.time() - wall_start
            entry.main_cpu += process_cpu() - cpu_start
            self.stack.pop()

    def report(self, stream: Any = None, depth: int = 3) -> None:
        """The stage tree to `depth` levels: evaluations, saved, pool tasks,
        wall and CPU (minutes)."""
        stream = stream if stream is not None else sys.stdout
        stream.write("%-58s %6s %8s %7s %7s %8s %8s\n" % ("stage", "calls", "evals", "saved", "wall_m", "cpu_m",
                                                          "s/eval"))
        for path in self.order:
            if path == "(top)" or path.count("/") >= depth:
                continue
            entry = self.records[path]
            cpu = entry.worker_cpu + entry.main_cpu
            stream.write("%-58s %6d %8d %7d %7.2f %8.2f %8.3f\n" % (
                path[-58:], entry.calls, entry.evaluations, entry.saved, entry.wall / 60.0, cpu / 60.0,
                cpu / max(entry.evaluations, 1)))
        stream.flush()

    def dump(self, file_path: str) -> None:
        with open(file_path, "w") as handle:
            json.dump({"order": self.order, "records": {k: v.as_dict() for k, v in self.records.items()}}, handle,
                      indent=1)


LEDGER = Ledger()


def timed_task(task: Tuple[Callable[[Any], Any], Any]) -> Tuple[Any, float]:
    """A pool task (function, argument) and the CPU it took in its worker."""
    function, argument = task
    start = process_cpu()
    value = function(argument)
    return value, process_cpu() - start


def pool_map(pool: Any, function: Callable[[Any], Any], arguments: Sequence[Any],
             evaluations: Optional[int] = None) -> List[Any]:
    """pool.map(function, arguments), timed into the ledger's open stage,
    one task at a time (CHUNK). `evaluations`: the residuals the tasks
    compute in all (default one a task)."""
    if len(arguments) == 0:
        return []
    results = pool.map(timed_task, [(function, argument) for argument in arguments], chunksize=CHUNK)
    values = [value for value, _ in results]
    LEDGER.add(evaluations=len(arguments) if evaluations is None else evaluations, pool_tasks=len(arguments),
               worker_cpu=sum(cpu for _, cpu in results))
    return values


def count_main(evaluations: int = 1) -> None:
    """Residuals computed in the main process (their CPU is in the stage's
    main_cpu)."""
    LEDGER.add(evaluations=evaluations)


def count_saved(evaluations: int) -> None:
    LEDGER.add(saved=evaluations)


def stage(name: str) -> Any:
    return LEDGER.stage(name)


def summary_by_kind(ledger: Optional[Ledger] = None) -> Dict[str, StageRecord]:
    """Exclusive totals by stage kind (a path's last component, digits and
    phone names dropped), from the leaves up: each evaluation counted once,
    in the deepest stage that holds it."""
    ledger = ledger if ledger is not None else LEDGER
    children: Dict[str, List[str]] = {}
    for path in ledger.order:
        if path == "(top)":
            continue
        parent = path.rsplit("/", 1)[0] if "/" in path else "(top)"
        children.setdefault(parent, []).append(path)
    kinds: Dict[str, StageRecord] = {}
    for path in ledger.order:
        if path == "(top)":
            continue
        entry = ledger.records[path]
        own = StageRecord()
        own.evaluations = entry.evaluations
        own.saved = entry.saved
        own.worker_cpu = entry.worker_cpu
        own.main_cpu = entry.main_cpu
        own.wall = entry.wall
        own.calls = entry.calls
        for child in children.get(path, []):
            sub = ledger.records[child]
            own.evaluations -= sub.evaluations
            own.saved -= sub.saved
            own.worker_cpu -= sub.worker_cpu
            own.main_cpu -= sub.main_cpu
            own.wall -= sub.wall
        kind = kind_of(path)
        total = kinds.setdefault(kind, StageRecord())
        total.calls += own.calls
        total.evaluations += own.evaluations
        total.saved += own.saved
        total.worker_cpu += own.worker_cpu
        total.main_cpu += own.main_cpu
        total.wall += own.wall
    return kinds


def kind_of(path: str) -> str:
    """A path's kind: its last component without numbers ("node rung 3" ->
    "node rung"), prefixed by its phone-level parent kind."""
    last = path.rsplit("/", 1)[-1]
    words = [word for word in last.split() if not word.replace(".", "").replace("-", "").isdigit()]
    return " ".join(words) if words else last


def report_kinds(stream: Any = None, ledger: Optional[Ledger] = None) -> None:
    stream = stream if stream is not None else sys.stdout
    kinds = summary_by_kind(ledger)
    stream.write("%-34s %8s %8s %8s %8s\n" % ("kind (exclusive)", "evals", "saved", "wall_m", "cpu_m"))
    for kind in sorted(kinds, key=wall_of(kinds), reverse=True):
        entry = kinds[kind]
        stream.write("%-34s %8d %8d %8.2f %8.2f\n" % (kind[:34], entry.evaluations, entry.saved, entry.wall / 60.0,
                                                      (entry.worker_cpu + entry.main_cpu) / 60.0))
    stream.flush()


def wall_of(kinds: Dict[str, StageRecord]) -> Callable[[str], float]:
    def key(kind: str) -> float:
        return kinds[kind].wall
    return key
