"""A speaker's vowel formants from ARCTIC labels: LPC F1-F3 over the middle
40% of each vowel by f123.py's method (validated on bdl): 10 kHz, pre-
emphasis 0.63, 25 ms frames every 10 ms, order 12, roots under 600 Hz wide;
medians per vowel. (At 16 kHz with order 12 the poles were too few: bdl's
/aa/ came out with F2 2825 Hz against f123's 1221.) Writes JSON {vowel: [F1, F2, F3]} for
organ_vopt3.py's TARGETS (its vowels: iy ih eh ae aa ah ao ow uh uw ey ax).
usage: vowel_targets.py REFDIR IDS OUT.json"""
import json
import os
import sys

import numpy
from scipy.linalg import solve_toeplitz
from scipy.signal import lfilter, resample_poly

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fidelity

WANTED = {"iy", "ih", "eh", "ae", "aa", "ah", "ao", "ow", "uh", "uw", "ey", "ax"}


def formants(piece):
    y = lfilter([1, -0.63], [1], resample_poly(piece, 10000, fidelity.RATE))
    out = []
    for i in range(0, len(y) - 250, 100):
        f = y[i:i + 250] * numpy.hamming(250)
        c = numpy.correlate(f, f, "full")[249:262]
        if c[0] <= 0:
            continue
        a = solve_toeplitz(c[:12], -c[1:13])
        roots = numpy.roots(numpy.r_[1, a])
        roots = roots[numpy.imag(roots) > 0]
        fr = numpy.angle(roots) * 10000 / (2 * numpy.pi)
        bw = -numpy.log(numpy.abs(roots)) * 10000 / numpy.pi
        q = sorted(a_ for a_, b_ in zip(fr, bw) if a_ > 150 and b_ < 600)[:3]
        if len(q) == 3:
            out.append(q)
    return out


def main():
    ref_dir, ids, out_path = sys.argv[1:4]
    found = {}
    for x, segments in fidelity.arctic_set(ref_dir, fidelity.id_range(ids)):
        for phone, start, end in segments:
            name = phone.lower()
            if name not in WANTED or end - start < 0.06:
                continue
            margin = 0.3 * (end - start)
            found.setdefault(name, []).extend(formants(x[int((start + margin) * fidelity.RATE):int((end - margin) * fidelity.RATE)]))
    table = {v: [float(numpy.median([f[i] for f in rows])) for i in range(3)] for v, rows in found.items() if len(rows) >= 5}
    for v in sorted(table):
        print("%-3s n %4d  %5.0f %5.0f %5.0f" % (v, len(found[v]), *table[v]))
    json.dump(table, open(out_path, "w"))


if __name__ == "__main__":
    main()
