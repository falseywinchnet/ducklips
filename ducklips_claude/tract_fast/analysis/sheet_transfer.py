"""The sheet's own transfer for a held vowel: radiated output over glottal
flow, at the harmonics, over 1.0-2.0 s (trace at every 96 kHz sample, U in
column 5; output WAV at 48 kHz with fixed gain). Against the ideal chain
(ideal_level.transfer of the sheet's areas, times f for radiation).
Prints both, normalised at 600-1000 Hz, in bands.
usage: sheet_transfer.py TRACE WAV AREAS"""
import sys, numpy
from scipy.io import wavfile
from scipy.signal import resample_poly
sys.path.insert(0, "/tmp")
from ideal_level import transfer

trace, wav, areas_path = sys.argv[1:4]
u = numpy.loadtxt(trace, usecols=(4,))
u = resample_poly(u, 1, 2)[48000:96000]
r, y = wavfile.read(wav)
y = y.astype(float)[48000:96000]
n = len(y)
w = numpy.hanning(n)
U = numpy.fft.rfft((u - u.mean()) * w)
Y = numpy.fft.rfft(y * w)
f = numpy.fft.rfftfreq(n, 1 / 48000)
rows = [l.split() for l in open(areas_path)]
start = max(i for i, row in enumerate(rows) if row[2] == "05_Vocal_cord") + 1
areas = numpy.array([float(row[4]) for row in rows[start:]])
# harmonics: peaks of |U| near multiples of the F0 found from U's spectrum
f0 = f[numpy.argmax(numpy.abs(U[(f > 70) & (f < 300)])) + numpy.argmax(f > 70)]
bands = [(80, 300), (300, 600), (600, 1000), (1000, 2000), (2000, 3500), (3500, 5000)]
sheet, ideal, source = {}, {}, {}
for h in numpy.arange(f0, 5000, f0):
    sel = (f > h - 0.3 * f0) & (f < h + 0.3 * f0)
    k = numpy.argmax(numpy.abs(U[sel])) + numpy.argmax(sel)
    t = numpy.abs(Y[k]) / (numpy.abs(U[k]) + 1e-12)
    ti = transfer(areas, numpy.array([f[k]]))[0] * f[k]
    for a, b in bands:
        if a <= f[k] < b:
            sheet.setdefault((a, b), []).append(20 * numpy.log10(t + 1e-30))
            source.setdefault((a, b), []).append(20 * numpy.log10(numpy.abs(U[k]) + 1e-30))
            ideal.setdefault((a, b), []).append(20 * numpy.log10(ti))
ref_s = numpy.mean(sheet[(600, 1000)]); ref_i = numpy.mean(ideal[(600, 1000)])
print("F0 %.0f  band: sheet / ideal (dB re 600-1000 Hz) / sheet-ideal / absolute sheet, ideal" % f0)
for b in bands:
    if b in sheet:
        s, i = numpy.mean(sheet[b]) - ref_s, numpy.mean(ideal[b]) - ref_i
        print("  %4d-%4d  %+6.1f  %+6.1f  %+6.1f   abs %6.1f %6.1f  source %6.1f" % (b[0], b[1], s, i, s - i, numpy.mean(sheet[b]), numpy.mean(ideal[b]), numpy.mean(source[b])))
r, g = wavfile.read(wav)
print("  output rms %.1f dB, glottal AC rms %.1f dB" % (20 * numpy.log10(y.std()), 20 * numpy.log10(u.std())))
