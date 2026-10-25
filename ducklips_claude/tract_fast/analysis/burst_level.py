"""Burst level of stops: the loudest 10 ms frame inside each instance of a
stop (from its start to 20 ms past its end), re the mean level of the
following vowel; bdl against ours over ARCTIC a0001-a0060.
usage: burst_level.py PHONE..."""
import os, sys, numpy
sys.path.insert(0, os.path.expanduser("~/ducklips/voice/tools"))
import compare_arctic as ca
VOWELS = {"aa","ae","ah","ao","aw","ay","eh","er","ey","ih","iy","ow","oy","uh","uw","ax"}
for phone in sys.argv[1:]:
    for who in ("bdl", "ours"):
        diffs = []
        for number in range(1, 61):
            name = "arctic_a%04d" % number
            if who == "ours":
                path = os.path.join(ca.OUT, name + "_ours.wav"); phones = ca.our_phones(os.path.join(ca.OUT, name + ".score"))
            else:
                path = os.path.join(ca.REF, "wav", name + ".wav"); phones = ca.reference_phones(name)
            if not os.path.exists(path): continue
            x, r = ca.read_wave(path); x = numpy.asarray(x, float)
            for i, (p, s, e) in enumerate(phones[:-1]):
                n, ns, ne = phones[i + 1]
                if p != phone or n not in VOWELS or ne - ns < 0.04: continue
                f = int(0.01 * r)
                frames = [numpy.mean(x[j:j + f] ** 2) for j in range(int(s * r), int((e + 0.02) * r) - f, f // 2)]
                vowel = numpy.mean(x[int(ns * r):int(ne * r)] ** 2)
                if frames and vowel > 0: diffs.append(10 * numpy.log10(max(frames) / vowel + 1e-20))
        print("%-3s %-4s n=%3d burst re vowel: median %+.1f dB, quartiles %+.1f %+.1f" % (phone, who, len(diffs), numpy.median(diffs), numpy.percentile(diffs, 25), numpy.percentile(diffs, 75)))
