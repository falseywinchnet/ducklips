"""Before a listener steers the engine: does it hear people as people and
our voices as synthetic? Each listener (listeners.py) scores:
- words: the NATO reference reader and the owner's first pass (people);
  m64, the kinematic source (Perfect Paul) and the m66 policy (ours);
- sentences: bdl's ARCTIC a0001-a0010 (people); the same sentences by
  m64, kinematic and m66 (ours).
Prints each set's median P(synthetic) and the area under the ROC curve of
people against each of ours (1: every clip of ours scored more synthetic
than every person's; 0.5: chance).
usage: listener_check.py VALDIR [NAME,NAME...]   (VALDIR holds m64/ kin/ m66/
arc_m64/ arc_kin/ arc_m66/ of WAVs; default: every listener)"""
import glob
import os
import sys

import numpy
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import listeners  # noqa: E402
import nato_snr  # noqa: E402

WORDS = nato_snr.nato_compare.WORDS
BDL = os.path.expanduser("~/ducklips/voice/ref/arctic_bdl_full/wav")


def area_under_curve(people, ours):
    """P(a clip of ours scores above a person's), ties counting half."""
    above = 0.0
    for score in ours:
        above += numpy.sum(score > people) + 0.5 * numpy.sum(score == people)
    return above / (len(people) * len(ours))


def folder_clips(folder):
    return [listeners.load(path) for path in sorted(glob.glob(os.path.join(folder, "*.wav")))]


def main():
    torch.set_num_threads(8)
    folder = sys.argv[1]
    sets = {
        "nato reader": [nato_snr.reference_nato(index) for index in range(len(WORDS))],
        "owner": [nato_snr.reference_owner(index, 0) for index in range(len(WORDS))],
        "m64 words": folder_clips(os.path.join(folder, "m64")),
        "kinematic words": folder_clips(os.path.join(folder, "kin")),
        "m66 words": folder_clips(os.path.join(folder, "m66")),
        "bdl": [listeners.load(os.path.join(BDL, "arctic_a%04d.wav" % number)) for number in range(1, 11)],
        "m64 sentences": folder_clips(os.path.join(folder, "arc_m64")),
        "kinematic sentences": folder_clips(os.path.join(folder, "arc_kin")),
        "m66 sentences": folder_clips(os.path.join(folder, "arc_m66")),
    }
    comparisons = ((("nato reader", "owner"), ("m64 words", "kinematic words", "m66 words")),
                   (("bdl",), ("m64 sentences", "kinematic sentences", "m66 sentences")))
    names = sys.argv[2].split(",") if len(sys.argv) > 2 else list(listeners.LISTENERS)
    for name in names:
        listener = listeners.Listener(name)
        scores = {key: listener.synthetic_probability(clips) for key, clips in sets.items()}
        print("%s: median P(synthetic)" % name)
        for key in sets:
            print("  %-20s %.3f  (%d clips)" % (key, numpy.median(scores[key]), len(scores[key])))
        for people_keys, our_keys in comparisons:
            people = numpy.concatenate([scores[key] for key in people_keys])
            print("  AUC vs %s: %s" % ("+".join(people_keys), "  ".join(
                "%s %.2f" % (key, area_under_curve(people, scores[key])) for key in our_keys)), flush=True)


if __name__ == "__main__":
    main()
