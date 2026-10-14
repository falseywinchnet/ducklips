"""Voice onset time of every voiceless stop in an utterance, from the body's
own trace (say with CONSTRICTION=1): release is the first step after a
sealed closure where the narrowest oral area exceeds 0.02 cm^2; voicing
begins at the first glottal closure (area below 0.03 cm^2) after it.
usage: vot.py "text" [say arguments...]"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

text = sys.argv[1]
environment = dict(os.environ, CONSTRICTION="1")
result = subprocess.run([os.path.join(ROOT, "build", "say"), text, "/tmp/vot.wav", "--dry"] + sys.argv[2:],
                        capture_output=True, text=True, env=environment, check=True)
rows = []
for line in result.stderr.splitlines():
    parts = line.split()
    if len(parts) >= 9 and parts[1] == "narrowest":
        rows.append((float(parts[0]), float(parts[2]), float(parts[8])))
closed_since = None
index = 0
while index < len(rows):
    time, narrow, glottal = rows[index]
    if narrow < 0.002:
        if closed_since is None:
            closed_since = time
    elif closed_since is not None:
        if narrow > 0.02 and time - closed_since > 0.02:
            release = time
            voiced = None
            for later_time, later_narrow, later_glottal in rows[index:]:
                if later_glottal < 0.03:
                    voiced = later_time
                    break
                if later_narrow < 0.002:
                    break
            gap_at_release = glottal
            if voiced is not None:
                print("release %.3f s  glottis %.3f cm2 at release  voice onset +%.0f ms" %
                      (release, gap_at_release, 1000 * (voiced - release)))
            closed_since = None
        elif narrow >= 0.02:
            closed_since = None
    index += 1
