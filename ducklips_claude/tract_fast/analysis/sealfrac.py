"""Share of nasal instances whose oral tract is sealed (narrowest tract
area under 0.01 cm^2) at 25, 50 and 75% of the phone (sheet_areas).
usage: sealfrac.py PHONE N"""
import glob, re, subprocess, sys
phone, count = sys.argv[1].upper(), int(sys.argv[2])
hits = {0.25: [], 0.5: [], 0.75: []}
done = 0
for score in sorted(glob.glob("/tmp/es/vt_cur/*.score")):
    for line in open(score):
        m = re.match(r"# %s-?\d? ([\d.]+) ([\d.]+)" % phone, line)
        if not m or done >= count:
            continue
        s, d = float(m.group(1)), float(m.group(2))
        for q in hits:
            rows = [l.split() for l in subprocess.run(["/tmp/tb/sheet_areas", score, "%.3f" % (s + q * d)], capture_output=True, text=True).stdout.splitlines()]
            tract = [float(r[4]) for r in rows if int(r[0]) >= 35]
            hits[q].append(min(tract) < 0.01)
        done += 1
print(phone, "n", done, " ".join("%d%%: sealed %.0f%%" % (100 * q, 100 * sum(v) / len(v)) for q, v in hits.items()))
