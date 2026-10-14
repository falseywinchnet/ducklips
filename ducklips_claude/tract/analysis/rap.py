"""Cycle perturbation from speak's CYCLES logs (true closures): RAP (3-cycle
relative average perturbation, removes linear trends) median and p90 %, and
jumps per 100 cycles where a period differs >5% from the mean of its two
neighbours. Also deviation from the intended pitch (median |cents|)."""
import sys, numpy
for path in sys.argv[1:]:
    d = numpy.loadtxt(path, ndmin=2); rap, dev, n = [], [], 0
    run = []
    def flush(run):
        global n
        if len(run) < 5: return
        r = numpy.array(run); p = r[:, 1]; n += len(p) - 2
        m = (p[:-2] + p[1:-1] + p[2:]) / 3
        rap.extend(list(numpy.abs(p[1:-1] - m) / m * 100))
        dev.extend(list(numpy.abs(1200 * numpy.log2((1000 / p) / r[:, 4]))))
    for row in d:
        if 2.5 < row[1] < 15: run.append(row)
        else: flush(run); run = []
    flush(run)
    rap = numpy.array(rap)
    print("%-16s cycles %4d RAP median %.2f%% p90 %.2f%%  jumps(>5%%) %.1f/100  |dev| %.0f c" % (
        path.split("/")[-1], n, numpy.median(rap), numpy.percentile(rap, 90), 100 * numpy.mean(rap > 5 * 2 / 3), numpy.median(dev)))
