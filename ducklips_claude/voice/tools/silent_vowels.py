"""Vowels that fail to sound: level more than 26 dB under the sentence's 30 ms peak, over our ARCTIC renders a0001-a0060. A guard for fold changes (stressed vowels once went silent at phrase starts)."""
import os, sys, numpy
sys.path.insert(0, "tools"); import compare_arctic as ca
silent = []; total = 0
for k in range(1, 61):
    name = "arctic_a%04d" % k
    d, r = ca.read_wave(os.path.join(ca.OUT, name + "_ours48.wav")); ph = ca.our_phones(os.path.join(ca.OUT, name + ".score"))
    peak = numpy.sqrt(numpy.max(numpy.convolve(d**2, numpy.ones(int(0.03*r))/int(0.03*r), "same")))
    for p, s, e in ph:
        if p not in ca.VOWELS: continue
        total += 1
        v = d[int(s*r):int(e*r)]
        if numpy.sqrt(numpy.mean(v**2)) < 0.05 * peak: silent.append((name, p, round(s, 2)))
print("vowels %d, under -26 dB re sentence peak: %d" % (total, len(silent)), silent[:12])
