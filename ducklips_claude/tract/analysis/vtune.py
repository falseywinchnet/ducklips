"""Vowel formants of the sheet in running speech under a phone table: plans
ARCTIC a0001-a0060 with say (--plan-only, pitch 120, the given table),
renders them through sheet_speak, and prints per vowel the median LPC F1/F2
(middle 40%) against bdl's own over the same sentences, with the distance
in semitones.
usage: vtune.py PHONES.txt NAME [VOWEL...]  (ENGINE=tract: the segment engine at 750 Pa)"""
import os, re, sys, subprocess, numpy
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, "/tmp")
from vw import load
from dips import phones_bdl, phones_ours
from scipy.signal import resample_poly, lfilter
from scipy.linalg import solve_toeplitz

HOME = os.path.expanduser("~")
REF = HOME + "/ducklips/voice/ref/arctic_bdl"
SAY = HOME + "/ducklips/voice/build/say"
SHEET = HOME + "/ducklips/tract/build/sheet_speak"


def lpc(y):
    out = []
    for i in range(0, len(y) - 250, 100):
        f = y[i:i + 250] * numpy.hamming(250); c = numpy.correlate(f, f, "full")[249:262]
        if c[0] <= 0:
            continue
        a = solve_toeplitz(c[:12], -c[1:13]); roots = numpy.roots(numpy.r_[1, a]); roots = roots[numpy.imag(roots) > 0]
        fr = numpy.angle(roots) * 10000 / (2 * numpy.pi); bw = -numpy.log(numpy.abs(roots)) * 10000 / numpy.pi
        q = sorted(a_ for a_, b_ in zip(fr, bw) if a_ > 150 and b_ < 800)[:3]
        if len(q) == 3:
            out.append(q)
    return numpy.median(out, axis=0) if out else None


def collect(source, want):
    got = {}
    for wav, phones in source:
        if not os.path.exists(wav):
            continue
        x = resample_poly(load(wav), 10000, 16000); x = lfilter([1, -0.63], [1], x)
        for name, s, e in phones:
            if name in want and e - s >= 0.06:
                r = lpc(x[int((s + 0.3 * (e - s)) * 10000):int((e - 0.3 * (e - s)) * 10000)])
                if r is not None:
                    got.setdefault(name, []).append(r)
    return {k: numpy.median(v, axis=0) for k, v in got.items()}, {k: len(v) for k, v in got.items()}


table, name = sys.argv[1], sys.argv[2]
want = set(sys.argv[3:]) or {"iy", "ih", "ey", "eh", "ae", "aa", "ah", "ao", "ow", "uh", "uw", "er"}
texts = {}
for line in open(REF + "/arctic.data"):
    m = re.match(r'\( (arctic_a00(?:[0-5]\d|60)) "(.*)" \)', line.strip())
    if m and m.group(1) != "arctic_a0000":
        texts[m.group(1)] = m.group(2)
out = "/tmp/es/vt_" + name
os.makedirs(out, exist_ok=True)


def job(item):
    key, text = item
    score, wav = os.path.join(out, key + ".score"), os.path.join(out, key + ".wav")
    if os.environ.get("ENGINE") == "tract":
        subprocess.run([SAY, text, wav, "--score", score, "--plan-only", "--pressure", "750", "--phones", table],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([HOME + "/ducklips/tract/build/speak", score, wav], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    else:
        subprocess.run([SAY, text, wav, "--score", score, "--plan-only", "--pitch", "120", "--phones", table],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([SHEET, score, wav], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


with ThreadPoolExecutor(8) as pool:
    list(pool.map(job, sorted(texts.items())))
cache = "/tmp/es/bdl_vowels60.npy"
if os.path.exists(cache):
    ref = numpy.load(cache, allow_pickle=True).item()
else:
    ref = collect(phones_bdl(REF + "/lab", REF + "/wav", 60), set(["iy", "ih", "ey", "eh", "ae", "aa", "ah", "ao", "ow", "uh", "uw", "er"]))[0]
    numpy.save(cache, ref)
ours, counts = collect(phones_ours(out, out), want)
total = []
for v in sorted(want):
    if v in ours and v in ref:
        d = numpy.hypot(12 * numpy.log2(ours[v][0] / ref[v][0]), 12 * numpy.log2(ours[v][1] / ref[v][1]))
        total.append(d)
        print("%-3s n %3d  ours %4.0f %4.0f  bdl %4.0f %4.0f  dist %.1f st" % (v, counts[v], ours[v][0], ours[v][1], ref[v][0], ref[v][1], d))
print("%s mean distance %.2f st" % (name, numpy.mean(total)))
