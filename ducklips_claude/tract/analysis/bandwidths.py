"""Formant bandwidths (Hz) by LPC (10 kHz, order 12, 25 ms frames, 10 ms
hop): for each frame the poles' frequencies and bandwidths; F1-F3 taken as
the three lowest poles above 150 Hz with bandwidth under 700 Hz; medians.
usage: bandwidths.py ours WAV:START:END ... | bandwidths.py arctic LABDIR WAVDIR VOWELS MAXFILES"""
import os, sys, numpy, soundfile, librosa
def frames(x, r):
    x = librosa.resample(x, orig_sr=r, target_sr=10000)
    x = numpy.append(x[0], x[1:] - 0.63 * x[:-1])
    out = []
    for i in range(0, len(x) - 250, 100):
        seg = x[i:i + 250] * numpy.hamming(250)
        if numpy.sqrt(numpy.mean(seg ** 2)) < 1e-4: continue
        a = librosa.lpc(seg, order=12)
        poles = [z for z in numpy.roots(a) if z.imag > 0]
        fb = sorted((numpy.angle(z) * 10000 / (2 * numpy.pi), -numpy.log(abs(z)) * 10000 / numpy.pi) for z in poles)
        fb = [(f, b) for f, b in fb if f > 150 and b < 700][:3]
        if len(fb) == 3: out.append([b for f, b in fb] + [f for f, b in fb])
    return out
def main():
    if sys.argv[1] == "ours":
        for spec in sys.argv[2:]:
            path, a, b = spec.split(":")
            x, r = soundfile.read(path)
            m = numpy.median(numpy.array(frames(x[int(float(a) * r):int(float(b) * r)], r)), axis=0)
            print("%-22s F %4.0f %4.0f %4.0f  bandwidth %3.0f %3.0f %3.0f Hz" % (path.split("/")[-1], m[3], m[4], m[5], m[0], m[1], m[2]))
        return
    labdir, wavdir, vowels, maxfiles = sys.argv[2], sys.argv[3], sys.argv[4].split(","), int(sys.argv[5])
    rows = {v: [] for v in vowels}
    for name in sorted(os.listdir(labdir))[:maxfiles]:
        lines = [l.split() for l in open(os.path.join(labdir, name)) if len(l.split()) == 3]
        wav = os.path.join(wavdir, name.replace(".lab", ".wav"))
        if not os.path.exists(wav): continue
        x, r = soundfile.read(wav); prev = 0.0
        for end, _, ph in lines:
            end = float(end)
            if ph in rows and end - prev >= 0.12:
                mid = (prev + end) / 2; half = (end - prev) * 0.3
                rows[ph].extend(frames(x[int((mid - half) * r):int((mid + half) * r)], r))
            prev = end
    for v in vowels:
        m = numpy.median(numpy.array(rows[v]), axis=0)
        print("bdl %-3s frames %-4d F %4.0f %4.0f %4.0f  bandwidth %3.0f %3.0f %3.0f Hz" % (v, len(rows[v]), m[3], m[4], m[5], m[0], m[1], m[2]))
if __name__ == "__main__":
    main()
