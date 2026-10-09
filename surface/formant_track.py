"""Tracks F1-F3 (LPC, 25 ms frames every 5 ms) through voiced stretches.

usage: python3 formant_track.py FILE.wav [from_s to_s]
"""
import sys

import numpy

exec(open(__file__.replace("formant_track.py", "voice_compare.py")).read().split("for argument in sys.argv")[0])


def main():
    signal, rate = read_wav(sys.argv[1])
    begin = float(sys.argv[2]) if len(sys.argv) > 3 else 0.0
    end = float(sys.argv[3]) if len(sys.argv) > 3 else len(signal) / rate
    frame = int(0.025 * rate)
    for start in range(int(begin * rate), int(end * rate) - frame, int(0.005 * rate)):
        piece = signal[start:start + frame]
        f0, strength = periodicity(piece, rate)
        if strength < 0.5:
            continue
        resonances = formant_bandwidths(piece, rate)
        text = " ".join("%5.0f" % frequency for frequency, bandwidth in resonances)
        print("%.3f  %s" % (start / rate, text))


main()
