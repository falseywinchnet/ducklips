"""Which phones the local clicks (clicks2.py, K) fall in: counts per phone
and per second of that phone. usage: clickphones.py WAV SCORE K"""
import sys, re, subprocess, collections
wav, score, K = sys.argv[1], sys.argv[2], sys.argv[3]
out = subprocess.run(["/Users/joshuahkuttenkuler/miniforge3/bin/python3", "/tmp/clicks2.py", wav, K], capture_output=True, text=True).stdout.split("\n")
times = [float(t) for t in out[1].split()] if len(out) > 1 and out[1].strip() else []
ph = [(m.group(1).rstrip("-0123456789"), float(m.group(2)), float(m.group(3))) for m in (re.match(r"# (\S+) ([\d.]+) ([\d.]+)", l) for l in open(score)) if m]
count = collections.Counter(); dur = collections.Counter()
for n, s, d in ph: dur[n] += d
for t in times:
    for n, s, d in ph:
        if s <= t < s + d: count[n] += 1
for n, c in count.most_common(12):
    print("%-4s clicks %3d  %.0f/s" % (n, c, c / dur[n]))
