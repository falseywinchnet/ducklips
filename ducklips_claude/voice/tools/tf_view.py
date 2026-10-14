"""Shows the four time-frequency views (plain and enhanced) for one or two analyses.
usage: tf_view.py OUT.png PREFIX [PREFIX2] [--plain]"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy

VIEWS = [("stft", "STFT (fixed window)"), ("cq", "constant-Q Gabor (Garces)"),
         ("causal", "time-causal Gabor (Lindeberg)"), ("wigner", "smoothed pseudo-Wigner-Ville")]


def load(prefix, view, kind):
    field = numpy.load("%s_%s_%s.npy" % (prefix, view, kind))
    hz = numpy.load("%s_%s_%s_hz.npy" % (prefix, view, kind))[0]
    seconds = numpy.load("%s_seconds.npy" % prefix)[0]
    return field, hz, seconds


def main():
    output = sys.argv[1]
    prefixes = [a for a in sys.argv[2:] if not a.startswith("--")]
    kind = "plain" if "--plain" in sys.argv else "enhanced"
    figure, axes = plt.subplots(len(VIEWS), len(prefixes), figsize=(7 * len(prefixes), 11), squeeze=False)
    for column, prefix in enumerate(prefixes):
        for row, (view, title) in enumerate(VIEWS):
            field, hz, seconds = load(prefix, view, kind)
            level = 20 * numpy.log10(field + 1e-9)
            top = numpy.percentile(level, 99.7)
            axis = axes[row][column]
            axis.pcolormesh(seconds, hz, level, vmin=top - 60, vmax=top, cmap="magma", shading="auto")
            axis.set_yscale("log")
            axis.set_yticks([100, 200, 500, 1000, 2000, 4000, 7000])
            axis.set_yticklabels(["100", "200", "500", "1k", "2k", "4k", "7k"])
            axis.set_title("%s - %s (%s)" % (prefix.split("/")[-1], title, kind), fontsize=8)
    figure.tight_layout()
    figure.savefig(output, dpi=80)


main()
