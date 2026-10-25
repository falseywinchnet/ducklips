"""Sheet cross-sections at the middle of instances of a phone: lip area and
the narrowest tract node, per instance (first N in vt_cur).
usage: midareas.py PHONE N"""
import glob, re, subprocess, sys
phone, count = sys.argv[1].upper(), int(sys.argv[2])
done = 0
for score in sorted(glob.glob("/tmp/es/vt_cur/*.score")):
    for line in open(score):
        m = re.match(r"# %s-?\d? ([\d.]+) ([\d.]+)" % phone, line)
        if not m:
            continue
        mid = float(m.group(1)) + float(m.group(2)) / 2
        rows = [l.split() for l in subprocess.run(["/tmp/tb/sheet_areas", score, "%.3f" % mid], capture_output=True, text=True).stdout.splitlines()]
        tract = [r for r in rows if int(r[0]) >= 35]
        narrow = min(tract, key=lambda r: float(r[4]))
        print("%s %s mid %.3f dur %s: lips %s cm^2, narrowest node %s (%s cm) %s cm^2" % (score.split("/")[-1], phone, mid, m.group(2), tract[-1][4], narrow[0], narrow[1], narrow[4]))
        done += 1
        if done >= count:
            sys.exit()
