"""Pretrained listeners that tell synthetic speech from people, as costs for
the engine policy (policy_fit.py). Each is a stock transformers audio
classifier (no code of its own; only config and safetensors weights,
downloaded to ~/models): a clip at 16 kHz goes in, the probability of its
"synthetic" label comes out.

- ast5:    MattyB95/AST-ASVspoof5-Synthetic-Voice-Detection (Audio
           Spectrogram Transformer fine-tuned on ASVspoof 5)
- ast2019: MattyB95/AST-ASVspoof2019-Synthetic-Voice-Detection (on
           ASVspoof 2019 LA)
- w2v_mo:  mo-thecreator/Deepfake-audio-detection (wav2vec2-base)
- w2v_mm:  MelodyMachine/Deepfake-audio-detection-V2 (wav2vec2-base)

usage: listeners.py WAV...   (prints each listener's synthetic probability)"""
import os
import sys

# (transformers would import the Mini's jax, whose native library is broken.)
os.environ.setdefault("USE_FLAX", "0")
os.environ.setdefault("USE_TF", "0")

import librosa  # noqa: E402
import numpy  # noqa: E402
import torch  # noqa: E402
from transformers import AutoConfig, AutoFeatureExtractor, AutoModelForAudioClassification  # noqa: E402

RATE = 16000
MODELS = os.path.expanduser("~/models")
LISTENERS = {
    "ast5": "AST-ASVspoof5-Synthetic-Voice-Detection",
    "ast2019": "AST-ASVspoof2019-Synthetic-Voice-Detection",
    "w2v_mo": "Deepfake-audio-detection",
    "w2v_mm": "Deepfake-audio-detection-V2",
}
SYNTHETIC_WORDS = ("spoof", "fake", "synthetic", "deepfake")


def synthetic_label(config):
    """The index of the label that means synthetic."""
    for index, name in config.id2label.items():
        if any(word in name.lower() for word in SYNTHETIC_WORDS):
            return int(index)
    raise ValueError("no synthetic label in %s" % config.id2label)


class Listener:
    def __init__(self, name):
        folder = os.path.join(MODELS, LISTENERS[name])
        self.name = name
        self.config = AutoConfig.from_pretrained(folder)
        self.extractor = AutoFeatureExtractor.from_pretrained(folder)
        self.model = AutoModelForAudioClassification.from_pretrained(folder, use_safetensors=True)
        self.model.eval()
        self.synthetic = synthetic_label(self.config)

    def synthetic_probability(self, signals):
        """P(synthetic) for each 16 kHz clip in a list."""
        probabilities = []
        with torch.no_grad():
            for signal in signals:
                inputs = self.extractor(signal, sampling_rate=RATE, return_tensors="pt")
                logits = self.model(**inputs).logits[0]
                probabilities.append(float(torch.softmax(logits, dim=-1)[self.synthetic]))
        return numpy.array(probabilities)


def load(path):
    signal, _ = librosa.load(path, sr=RATE)
    return signal


def main():
    torch.set_num_threads(4)
    signals = [load(path) for path in sys.argv[1:]]
    for name in LISTENERS:
        listener = Listener(name)
        probabilities = listener.synthetic_probability(signals)
        print("%-8s labels %s" % (name, listener.config.id2label))
        for path, probability in zip(sys.argv[1:], probabilities):
            print("  %.3f  %s" % (probability, os.path.basename(path)))


if __name__ == "__main__":
    main()
