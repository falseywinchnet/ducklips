"""One evaluation pass: render ARCTIC a0001-a0060 with every view analysed
(three seeds: the wander and turbulence are random), then report the four-view phone identification (tools/phone_test.py over
stft, cq, causal and wigner; higher is better, human ceiling about 55%) and
the four-view distance on a0001-a0020 and a0021-a0040 (lower is better).
Extra say arguments come from SAY_ARGS.
usage: evaluate.py"""
import contextlib
import io
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca


def distance(first, last):
    folder = os.path.join(ca.ROOT, "out", "evaluate_%d_%d" % (first, last))
    records = ca.compare(first, last, folder, os.environ.get("SAY_ARGS", "").split())
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink):
        total = ca.four_view_summary(records)
    return total, sink.getvalue()


SEEDS = (3, 2, 1)


def identification(tools, seed):
    """Renders every utterance at one seed and returns the phone test report."""
    extra = os.environ.get("SAY_ARGS", "").split() + ["--seed", str(seed)]
    environment = dict(os.environ, SAY_ARGS=" ".join(extra))
    subprocess.run([sys.executable, os.path.join(tools, "render_arctic.py"), "--analyse"], check=True, env=environment)
    return subprocess.run([sys.executable, os.path.join(tools, "phone_test.py"), "--views", "stft,cq,causal,wigner"],
                          check=True, capture_output=True, text=True).stdout


def main():
    tools = os.path.dirname(os.path.abspath(__file__))
    scores = []
    for seed in SEEDS:
        test = identification(tools, seed)
        scores.append(float(test.strip().splitlines()[-1].split()[-1]))
    # The last seed (1) stays in out/arctic for the diagnostic tools.
    open(os.path.join(ca.ROOT, "out", "phone_test_four_view.txt"), "w").write(test)
    train, train_text = distance(1, 20)
    held, held_text = distance(21, 40)
    print(train_text, end="")
    for line in test.splitlines():
        if line.startswith("ours:") or line.startswith("PHONE_TEST") or "affricate" in line:
            print(line)
    print("EVALUATE identification %.1f (seeds %s)  distance a1-20 %.3f  a21-40 %.3f" %
          (sum(scores) / len(scores), " ".join("%.1f" % v for v in scores), train, held))


if __name__ == "__main__":
    main()
