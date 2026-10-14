"""Vowels in Hillenbrand's own context: each hVd word is spoken alone and its
vowel measured (LPC F1-F3, middle 40%) against the 1995 adult male means.
Separates the planner's running-speech effects from the engine's vowel shapes.
usage: hvd_test.py [say arguments...]"""
import os
import subprocess
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca
from invert import lpc_formants

WORDS = [("heed", "iy", 342, 2322, 3000), ("hid", "ih", 427, 2034, 2684), ("hayed", "ey", 476, 2089, 2691),
         ("head", "eh", 580, 1799, 2605), ("had", "ae", 588, 1952, 2601), ("hod", "aa", 768, 1333, 2522),
         ("hawed", "ao", 652, 997, 2538), ("hoed", "ow", 497, 910, 2459), ("hood", "uh", 469, 1122, 2434),
         ("who'd", "uw", 378, 997, 2343), ("hud", "ah", 623, 1200, 2550), ("heard", "er", 474, 1379, 1710)]


def measure(word, vowel):
    folder = os.path.join(ca.ROOT, "out", "hvd")
    os.makedirs(folder, exist_ok=True)
    wave = os.path.join(folder, word.replace("'", "") + ".wav")
    score = os.path.join(folder, word.replace("'", "") + ".score")
    subprocess.run([os.path.join(ca.ROOT, "build", "say"), word.capitalize() + ".", wave, "--score", score, "--dry"] +
                   sys.argv[1:], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    data, rate = ca.read_wave(wave)
    data = ca.resample(data, rate, 16000)
    phones = ca.our_phones(score)
    found = [(s, e) for p, s, e in phones if p == vowel]
    if not found:
        return " ".join(p for p, s, e in phones), None
    start, end = found[0]
    values = []
    for t in numpy.linspace(start + 0.3 * (end - start), end - 0.3 * (end - start), 5):
        segment = data[int((t - 0.0125) * 16000):int((t + 0.0125) * 16000)]
        f = lpc_formants(segment, 16000)
        if f:
            values.append(f)
    return " ".join(p for p, s, e in phones), numpy.median(values, axis=0) if values else None


if __name__ == "__main__":
    total = []
    for word, vowel, f1, f2, f3 in WORDS:
        phones, f = measure(word, vowel)
        if f is None:
            print("%-6s %-14s no measurement" % (word, phones))
            continue
        error = numpy.abs(numpy.log(numpy.array(f[:2]) / numpy.array([f1, f2])))
        total.append(error.mean())
        print("%-6s %-14s target %4d %4d %4d   ours %4.0f %4.0f %4.0f" % (word, phones, f1, f2, f3, *f))
    print("HVD mean |log error| F1,F2: %.3f" % numpy.mean(total))
