"""A perceptual cost for the engine policy: how far our reading of an ARCTIC
sentence lies from bdl's real recording of it, in WavLM's representation
of speech (microsoft/wavlm-base-plus, stock transformers, weights in
~/models).

Both recordings are trimmed (35 dB), normalized, and run through WavLM;
the frames of layers 6 and 9 (unit length) are aligned by dynamic time
warping on cosine distance, and the distance is the mean cosine distance
along the path, averaged over the two layers.

MEASURED (a0001-a0010): rms, another person reading the same sentences,
lies 0.283 / 0.156 (layers 6 / 9) from bdl; m64 0.588 / 0.489, the
kinematic source 0.564 / 0.469, the m66 policy 0.571 / 0.473. rms is
nearer bdl than each of ours in every sentence.

Sentences: FIT (a0001-a0010) steer; HELD (a0011-a0020) are never fitted and
check that the policy learns speech, not ten recordings."""
import os
import re
import subprocess

os.environ.setdefault("USE_FLAX", "0")
os.environ.setdefault("USE_TF", "0")

import librosa  # noqa: E402
import numpy  # noqa: E402
import torch  # noqa: E402
from transformers import WavLMModel  # noqa: E402

import engine_net  # noqa: E402

ROOT = engine_net.ROOT
RATE = 16000
LAYERS = (6, 9)
MODEL = os.path.expanduser("~/models/wavlm-base-plus")
BDL = os.path.join(ROOT, "voice", "ref", "arctic_bdl_full", "wav")
TEXTS = os.path.join(ROOT, "voice", "ref", "arctic_bdl", "arctic.data")
FIT = ["arctic_a%04d" % number for number in range(1, 11)]
HELD = ["arctic_a%04d" % number for number in range(11, 21)]

# Per worker process: the model and bdl's frames, loaded once.
WORKER = {}


def sentence_texts():
    texts = {}
    with open(TEXTS) as handle:
        for line in handle:
            found = re.match(r'\( (\S+) "(.*)" \)', line.strip())
            if found:
                texts[found.group(1)] = found.group(2)
    return texts


def prepared(signal):
    trimmed, _ = librosa.effects.trim(signal, top_db=35)
    return (trimmed - trimmed.mean()) / (trimmed.std() + 1e-9)


def frames(signal):
    """{layer: unit-length frames} of a prepared clip."""
    with torch.no_grad():
        hidden = WORKER["model"](torch.tensor(signal, dtype=torch.float32)[None], output_hidden_states=True).hidden_states
    layers = {}
    for layer in LAYERS:
        values = hidden[layer][0].numpy()
        layers[layer] = values / (numpy.linalg.norm(values, axis=1, keepdims=True) + 1e-9)
    return layers


def distance(first, second):
    total = 0.0
    for layer in LAYERS:
        cost = 1.0 - first[layer] @ second[layer].T
        accumulated, path = librosa.sequence.dtw(C=cost)
        total += float(accumulated[-1, -1]) / len(path)
    return total / len(LAYERS)


def start_worker():
    """Pool initializer: one thread each, the model and bdl's frames."""
    torch.set_num_threads(1)
    WORKER["model"] = WavLMModel.from_pretrained(MODEL)
    WORKER["model"].eval()
    WORKER["texts"] = sentence_texts()
    WORKER["bdl"] = {}
    for name in FIT + HELD:
        signal, _ = librosa.load(os.path.join(BDL, name + ".wav"), sr=RATE)
        WORKER["bdl"][name] = frames(prepared(signal))


def render_sentence(task):
    """task = (workdir, sentence, (body, organs, plan) paths, seed, tag).
    Returns (tag, sentence, distance to bdl or None)."""
    workdir, name, (tree_path, organs_path, plan_path), seed, tag = task
    folder = os.path.join(workdir, "renders", tag)
    os.makedirs(folder, exist_ok=True)
    score = os.path.join(folder, name + ".score")
    wave = os.path.join(folder, name + ".wav")
    result = (tag, name, None)
    planned = subprocess.run([os.path.join(ROOT, "voice", "build", "say"), WORKER["texts"][name], "/dev/null",
                              "--score", score, "--plan-only", "--pressure", "750", "--phones",
                              os.path.join(ROOT, "tract", "data", "phones_tract.txt"), "--organs", organs_path,
                              "--plan-params", plan_path],
                             cwd=os.path.join(ROOT, "voice"), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if planned.returncode == 0:
        environment = dict(os.environ)
        environment["FOLDTUBE"] = "1"
        environment["PARAMS"] = tree_path
        finished = subprocess.run([os.path.join(ROOT, "tract", "build", "speak"), score, wave, seed],
                                  env=environment, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if finished.returncode == 0:
            signal, _ = librosa.load(wave, sr=RATE)
            if numpy.all(numpy.isfinite(signal)) and numpy.abs(signal).max() > 0.0:
                result = (tag, name, distance(frames(prepared(signal)), WORKER["bdl"][name]))
    for path in (score, wave):
        if os.path.exists(path):
            os.remove(path)
    return result


def costs(pool, workdir, net, candidates, sentences):
    """candidates: [(tag, latents, seed)]. Returns {tag: mean distance over
    the sentences (FAILED when any failed)}."""
    tasks = []
    for tag, latents, seed in candidates:
        paths = net.write(latents, os.path.join(workdir, "net_" + tag))
        for name in sentences:
            tasks.append((workdir, name, paths, seed, tag))
    distances = {tag: [] for tag, _, _ in candidates}
    for tag, _, value in pool.imap_unordered(render_sentence, tasks):
        distances[tag].append(value)
    results = {}
    for tag, values in distances.items():
        results[tag] = engine_net.fit_blocks.FAILED_COST if any(value is None for value in values) else float(
            numpy.mean(values))
    return results
