"""Policy optimization of the engine network (engine_net.py) against the
NATO reference recordings: the neurons' latents are the policy, the engine
is the environment, -cost is the reward. Nothing is differentiated through
the engine; gradients are estimated from renders.

adamw (or wavlm, the same optimizer on wavlm_objective.py's cost): every
neuron at once. Each step draws PAIRS directions in latent
space and renders the policy nudged both ways along each (antithetic, the
same seed within a pair, a different seed for each pair, so the engine's
chance averages out). The costs are rank-shaped; their gradient feeds
AdamW. Dampers:
- each latent moves at most STEP_CAP per step (Adam's step is clipped), so a
  latent unit takes at least 1 / STEP_CAP steps;
- decoupled weight decay pulls every latent back toward 0 (the start);
- the neurons' activations saturate toward their bounds;
- every EVAL_EVERY steps the policy itself is rendered at two seeds; if its
  cost rose more than ROLLBACK over the best, it returns to the best, the
  learning rate halves and the moments clear.

als: block by block (each block of the parameter tree, then each vowel's
posture; ORDER=alternate puts a timing block before each shape block), the
rest held. The Jacobian of the residuals by forward differences in latent
space (DIFFERENCE; every render re-aligned to the reference, so the timing
residuals follow durations; the pooled residual keeps its length); a
Levenberg-Marquardt step on
mean(r^2) + DECAY * sum(z^2), each latent's step capped by the block's
radius; kept only if the rendered cost falls (else lambda x4, radius /2,
up to TRIES), as in fit_blocks.py.

Files in WORKDIR: reference.pkl, policy.json (latents, moments, rates;
resume), best.json, policy.log (one JSON line per step or block turn),
best.params, best.organs and best.plan (the best policy's engine files).
GROUPS=organs,plan (environment) fits those groups only (default body,organs).
usage: policy_fit.py WORKDIR adamw [STEPS]   (NATO reference, frame features)
       policy_fit.py WORKDIR wavlm [STEPS]   (bdl's ARCTIC, WavLM distance; wavlm_objective.py)
       policy_fit.py WORKDIR als [SWEEPS]
       policy_fit.py WORKDIR sentence OUT.wav   (the owner's sentence at the best policy)"""
import json
import multiprocessing
import os
import subprocess
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import engine_net  # noqa: E402
import fit_blocks  # noqa: E402
import wavlm_objective  # noqa: E402

FIT_SEED = fit_blocks.FIT_SEED
CHECK_SEED = fit_blocks.CHECK_SEED
SEEDS = ("73021", "11939", "40487", "68023", "25147", "90001", "33331", "57179")
WORKERS = 10
FAILED_COST = fit_blocks.FAILED_COST

# adamw
SIGMA = 0.2
PAIRS = 6
LEARNING_RATE = 0.03
STEP_CAP = 0.05
DECAY = 0.01
FIRST_MOMENT = 0.9
SECOND_MOMENT = 0.99
EVAL_EVERY = 5
ROLLBACK = 0.03

# als
DIFFERENCE = 0.3
LAMBDA_START = 1.0
LAMBDA_FLOOR = 0.05
TRIES = 3
RADIUS_START = 0.2
RADIUS_MAX = 0.4
RADIUS_MIN = 0.02
ALS_DECAY = 0.01

SENTENCE = "Hello, my name is Papa, and I speak with a body. The quick brown fox jumps over the lazy dog."


def candidate_cost(references, rendered):
    vector, whole, _ = engine_net.score(references, rendered)
    return engine_net.behaviour_cost(vector, whole), engine_net.groups_of(vector)


def load_policy(workdir, net):
    path = os.path.join(workdir, "policy.json")
    if os.path.exists(path):
        with open(path) as handle:
            state = json.load(handle)
        for key in ("latents", "first", "second"):
            state[key] = numpy.array(state[key])
        return state
    zeros = numpy.zeros(net.size())
    latents = zeros.copy()
    if os.environ.get("START_FROM"):
        # START_FROM=best.json of an earlier run over the same neurons.
        with open(os.environ["START_FROM"]) as handle:
            latents = numpy.array(json.load(handle)["latents"])
    return {"step": 0, "sweep": 0, "turn": 0, "latents": latents, "first": zeros.copy(), "second": zeros.copy(),
            "rate": LEARNING_RATE, "best": None, "damping": {}, "radius": {}}


def save_policy(workdir, net, state):
    path = os.path.join(workdir, "policy.json")
    plain = dict(state)
    for key in ("latents", "first", "second"):
        plain[key] = state[key].tolist()
    with open(path + ".tmp", "w") as handle:
        json.dump(plain, handle)
    os.replace(path + ".tmp", path)


def keep_best(workdir, net, state, cost):
    state["best"] = cost
    with open(os.path.join(workdir, "best.json"), "w") as handle:
        json.dump({"cost": cost, "latents": state["latents"].tolist(),
                   "values": net.values(state["latents"])}, handle, indent=1)
    net.write(state["latents"], os.path.join(workdir, "best"))


def policy_cost(pool, workdir, net, references, latents):
    """The policy itself at the fit seed and the check seed: (mean cost,
    [cost per seed], groups at the fit seed)."""
    rendered = engine_net.render_latents(pool, workdir, net, [("fit", latents, FIT_SEED), ("check", latents, CHECK_SEED)])
    fit, groups = candidate_cost(references, rendered["fit"])
    check, _ = candidate_cost(references, rendered["check"])
    return 0.5 * (fit + check), [fit, check], groups


def rank_utilities(costs):
    """Centered ranks in [-0.5, 0.5]: the highest cost gets +0.5."""
    order = numpy.argsort(numpy.argsort(costs))
    return order / (len(costs) - 1.0) - 0.5


class NatoObjective:
    """The NATO reference recordings, frame features (fit_blocks.py)."""

    def __init__(self, pool, workdir, net, references):
        self.pool = pool
        self.workdir = workdir
        self.net = net
        self.references = references

    def candidates(self, candidates):
        rendered = engine_net.render_latents(self.pool, self.workdir, self.net, candidates)
        return numpy.array([candidate_cost(self.references, rendered[tag])[0] for tag, _, _ in candidates])

    def policy(self, latents):
        return policy_cost(self.pool, self.workdir, self.net, self.references, latents)


class WavlmObjective:
    """bdl's ARCTIC recordings, WavLM distance (wavlm_objective.py): the
    policy is judged on the fitted sentences; the held-out ones are
    reported beside them, never used."""

    def __init__(self, pool, workdir, net):
        self.pool = pool
        self.workdir = workdir
        self.net = net

    def candidates(self, candidates):
        found = wavlm_objective.costs(self.pool, self.workdir, self.net, candidates, wavlm_objective.FIT)
        return numpy.array([found[tag] for tag, _, _ in candidates])

    def policy(self, latents):
        fit = wavlm_objective.costs(self.pool, self.workdir, self.net, [("fit", latents, FIT_SEED)],
                                    wavlm_objective.FIT)["fit"]
        held = wavlm_objective.costs(self.pool, self.workdir, self.net, [("held", latents, FIT_SEED)],
                                     wavlm_objective.HELD)["held"]
        return fit, [fit, held], {"fit": fit, "held": held}


def adamw(objective, workdir, net, steps):
    state = load_policy(workdir, net)
    log = open(os.path.join(workdir, "policy.log"), "a")
    generator = numpy.random.default_rng(int(FIT_SEED) + state["step"])
    if state["best"] is None:
        cost, per_seed, groups = objective.policy(state["latents"])
        keep_best(workdir, net, state, cost)
        log.write(json.dumps({"step": 0, "policy": cost, "seeds": per_seed, "groups": groups}) + "\n")
        print("start: policy cost %.4f (%s)" % (cost, " ".join("%.4f" % value for value in per_seed)), flush=True)
    for _ in range(steps):
        state["step"] += 1
        directions = generator.standard_normal((PAIRS, net.size()))
        candidates = []
        for pair in range(PAIRS):
            seed = SEEDS[(state["step"] * PAIRS + pair) % len(SEEDS)]
            candidates.append(("p%d" % pair, state["latents"] + SIGMA * directions[pair], seed))
            candidates.append(("m%d" % pair, state["latents"] - SIGMA * directions[pair], seed))
        costs = objective.candidates(candidates)
        utilities = rank_utilities(costs)
        gradient = numpy.zeros(net.size())
        for pair in range(PAIRS):
            gradient += (utilities[2 * pair] - utilities[2 * pair + 1]) * directions[pair]
        gradient /= 2.0 * PAIRS * SIGMA
        state["first"] = FIRST_MOMENT * state["first"] + (1.0 - FIRST_MOMENT) * gradient
        state["second"] = SECOND_MOMENT * state["second"] + (1.0 - SECOND_MOMENT) * gradient ** 2
        first = state["first"] / (1.0 - FIRST_MOMENT ** state["step"])
        second = state["second"] / (1.0 - SECOND_MOMENT ** state["step"])
        step = numpy.clip(state["rate"] * first / (numpy.sqrt(second) + 1e-8), -STEP_CAP, STEP_CAP)
        state["latents"] = state["latents"] - step - state["rate"] * DECAY * state["latents"]
        record = {"step": state["step"], "population": float(costs.mean()), "lowest": float(costs.min()),
                  "rate": state["rate"], "moved": float(numpy.abs(step).max()),
                  "distance": float(numpy.sqrt(numpy.mean(state["latents"] ** 2)))}
        line = "%4d population %.4f (lowest %.4f) rms latent %.3f" % (
            state["step"], record["population"], record["lowest"], record["distance"])
        if state["step"] % EVAL_EVERY == 0:
            cost, per_seed, groups = objective.policy(state["latents"])
            record.update({"policy": cost, "seeds": per_seed, "groups": groups})
            line += "  policy %.4f (%s)" % (cost, " ".join("%.4f" % value for value in per_seed))
            if cost < state["best"]:
                keep_best(workdir, net, state, cost)
                line += " best"
            elif cost > state["best"] * (1.0 + ROLLBACK):
                with open(os.path.join(workdir, "best.json")) as handle:
                    state["latents"] = numpy.array(json.load(handle)["latents"])
                state["rate"] *= 0.5
                state["first"][:] = 0.0
                state["second"][:] = 0.0
                record["rolled_back"] = True
                line += " rolled back, rate %.4f" % state["rate"]
        log.write(json.dumps(record) + "\n")
        log.flush()
        save_policy(workdir, net, state)
        print(line, flush=True)


FULL_COST_SLACK = 0.01  # a view's step may raise the whole cost by at most this share


def view_of(block):
    """Which residuals a block is solved against: with ORDER=alternate a
    timing block sees the timing residuals (every render re-aligned) and a
    shape block the frame residuals and pitch (the alignment held at the
    block's start); otherwise every residual (re-aligned)."""
    if os.environ.get("ORDER") != "alternate":
        return "all"
    return "timing" if block in TIMING_BLOCKS else "shape"


def view_score(references, rendered, view, alignments):
    if view == "shape":
        return engine_net.score(references, rendered, alignments)
    return engine_net.score(references, rendered)


def als_turn(pool, workdir, net, references, state, block, log):
    indices = net.block_indices(block)
    latents = state["latents"]
    view = view_of(block)
    candidates = [("base", latents, FIT_SEED)]
    for index in indices:
        nudged = latents.copy()
        nudged[index] += DIFFERENCE
        candidates.append(("d%d" % index, nudged, FIT_SEED))
    rendered = engine_net.render_latents(pool, workdir, net, candidates)
    full, whole, alignments = engine_net.score(references, rendered["base"])
    mask = engine_net.residual_mask(view, len(full))
    base = full[mask]
    decay = ALS_DECAY * float(numpy.sum(latents ** 2))
    base_full = engine_net.behaviour_cost(full, whole) + decay
    base_view = engine_net.behaviour_cost(base, whole) + decay
    jacobian = numpy.zeros((len(base), len(indices)))
    for column, index in enumerate(indices):
        vector, nudged_whole, _ = view_score(references, rendered["d%d" % index], view, alignments)
        if nudged_whole:
            jacobian[:, column] = (vector[mask] - base) / DIFFERENCE
    count = len(base)
    normal = jacobian.T @ jacobian / count
    gradient = jacobian.T @ base / count + ALS_DECAY * latents[indices]
    damping = state["damping"].get(block, LAMBDA_START)
    radius = state["radius"].get(block, RADIUS_START)
    record = {"sweep": state["sweep"], "block": block, "view": view, "base": base_full, "base_view": base_view,
              "tries": []}
    after = base_full
    for _ in range(TRIES):
        system = normal + damping * numpy.diag(numpy.diag(normal)) + ALS_DECAY * numpy.eye(len(indices))
        step = numpy.clip(numpy.linalg.solve(system, -gradient), -radius, radius)
        candidate = latents.copy()
        candidate[indices] += step
        tried = engine_net.render_latents(pool, workdir, net, [("step", candidate, FIT_SEED)])["step"]
        candidate_decay = ALS_DECAY * float(numpy.sum(candidate ** 2))
        vector, tried_whole, _ = view_score(references, tried, view, alignments)
        cost_view = engine_net.behaviour_cost(vector[mask], tried_whole) + candidate_decay
        whole_vector, whole_ok, _ = engine_net.score(references, tried)
        cost_full = engine_net.behaviour_cost(whole_vector, whole_ok) + candidate_decay
        predicted = float(numpy.mean((base + jacobian @ step) ** 2)) + candidate_decay
        record["tries"].append({"damping": damping, "radius": radius, "view_cost": cost_view, "predicted": predicted,
                                "cost": cost_full})
        if cost_view < base_view and cost_full <= base_full * (1.0 + FULL_COST_SLACK):
            state["latents"] = candidate
            damping = max(damping * 0.5, LAMBDA_FLOOR)
            radius = min(radius * 1.25, RADIUS_MAX)
            after = cost_full
            record["kept"] = cost_full
            break
        damping *= 4.0
        radius = max(radius * 0.5, RADIUS_MIN)
    state["damping"][block] = damping
    state["radius"][block] = radius
    log.write(json.dumps(record) + "\n")
    log.flush()
    return base_full, after


TIMING_BLOCKS = ("plan.duration", "plan.timing")


def als_schedule(net):
    """The blocks of one sweep. ORDER=alternate (environment): a timing
    block (TIMING_BLOCKS in turn) before each shape block (the consonant
    shapes, each vowel's posture), so timing and shape converge together;
    otherwise every block once, in the network's order."""
    if os.environ.get("ORDER") != "alternate":
        return list(net.blocks)
    timing = [block for block in TIMING_BLOCKS if block in net.blocks]
    shape = [block for block in net.blocks if block == "plan.articulation" or block.startswith("organs.")]
    schedule = []
    for index, block in enumerate(shape):
        if timing:
            schedule.append(timing[index % len(timing)])
        schedule.append(block)
    return schedule


def als(pool, workdir, net, references, sweeps):
    state = load_policy(workdir, net)
    schedule = als_schedule(net)
    log = open(os.path.join(workdir, "policy.log"), "a")
    if state["best"] is None:
        cost, per_seed, groups = policy_cost(pool, workdir, net, references, state["latents"])
        keep_best(workdir, net, state, cost)
        print("start: policy cost %.4f (%s)" % (cost, " ".join("%.4f" % value for value in per_seed)), flush=True)
    for _ in range(sweeps):
        while state["turn"] < len(schedule):
            block = schedule[state["turn"]]
            before, after = als_turn(pool, workdir, net, references, state, block, log)
            state["turn"] += 1
            save_policy(workdir, net, state)
            print("sweep %d %-16s %.4f -> %.4f" % (state["sweep"], block, before, after), flush=True)
        cost, per_seed, groups = policy_cost(pool, workdir, net, references, state["latents"])
        log.write(json.dumps({"sweep": state["sweep"], "policy": cost, "seeds": per_seed, "groups": groups}) + "\n")
        log.flush()
        print("sweep %d done: policy cost %.4f (%s)" % (state["sweep"], cost, " ".join("%.4f" % v for v in per_seed)),
              flush=True)
        if cost < state["best"]:
            keep_best(workdir, net, state, cost)
        state["sweep"] += 1
        state["turn"] = 0
        save_policy(workdir, net, state)


def sentence(workdir, out):
    score = os.path.join(workdir, "sentence.score")
    subprocess.run([os.path.join(engine_net.ROOT, "voice", "build", "say"), SENTENCE, "/dev/null", "--score", score,
                    "--plan-only", "--pressure", "750", "--phones",
                    os.path.join(engine_net.ROOT, "tract", "data", "phones_tract.txt"), "--organs",
                    os.path.join(workdir, "best.organs"), "--plan-params", os.path.join(workdir, "best.plan")],
                   cwd=os.path.join(engine_net.ROOT, "voice"), check=True, stdout=subprocess.DEVNULL)
    environment = dict(os.environ)
    # ENGINE_SOURCE=kinematic: the kinematic folds (VocalTractLab 2025 glottis, 2026-10-18) instead of
    # the fold tube.
    if os.environ.get("ENGINE_SOURCE") == "kinematic":
        environment.pop("FOLDTUBE", None)
        environment["KINEMATIC"] = "1"
    else:
        environment["FOLDTUBE"] = "1"
    environment["PARAMS"] = os.path.join(workdir, "best.params")
    subprocess.run([os.path.join(engine_net.ROOT, "tract", "build", "speak"), score, os.path.abspath(out), FIT_SEED],
                   env=environment, check=True)


def main():
    workdir = os.path.abspath(sys.argv[1])
    mode = sys.argv[2]
    os.makedirs(workdir, exist_ok=True)
    if mode == "sentence":
        sentence(workdir, sys.argv[3])
        return
    # GROUPS=body,organs,plan: the neuron groups being fitted (engine_net.py).
    net = engine_net.EngineNet(os.environ.get("GROUPS", "body,organs").split(","))
    print("%d neurons in %d blocks" % (net.size(), len(net.blocks)), flush=True)
    if mode == "wavlm":
        with multiprocessing.Pool(WORKERS, initializer=wavlm_objective.start_worker) as pool:
            adamw(WavlmObjective(pool, workdir, net), workdir, net, int(sys.argv[3]) if len(sys.argv) > 3 else 50)
        return
    references = fit_blocks.reference_frames(workdir)
    with multiprocessing.Pool(WORKERS) as pool:
        if mode == "adamw":
            adamw(NatoObjective(pool, workdir, net, references), workdir, net,
                  int(sys.argv[3]) if len(sys.argv) > 3 else 50)
        elif mode == "als":
            als(pool, workdir, net, references, int(sys.argv[3]) if len(sys.argv) > 3 else 1)


if __name__ == "__main__":
    main()
