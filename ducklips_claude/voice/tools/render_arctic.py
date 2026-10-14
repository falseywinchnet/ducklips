"""Render ARCTIC transcripts a0001-a0060 in parallel at 16 kHz (no analysis),
for formant and duration surveys. Extra say arguments come from SAY_ARGS.
ENGINE=sheet renders through the sheet tube (../tract/build/sheet_speak); PHONES
replaces its phone table (../tract/data/phones_sheet.txt). ENGINE=tract renders
through the segment-module engine (../tract/build/speak) at 750 Pa with
../tract/data/phones_tract.txt (or PHONES).
--analyse also writes all four view fields that tools/phone_test.py reads.
usage: render_arctic.py [FIRST LAST] [--analyse]"""
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_arctic as ca


def render(name, extra, analyse):
    texts = ca.transcripts()
    wave48 = os.path.join(ca.OUT, name + "_ours48.wav")
    score = os.path.join(ca.OUT, name + ".score")
    if os.environ.get("ENGINE") == "sheet":
        # The sheet tube speaks the front end's score (base pitch in its range;
        # SHEET holds its own arguments, a --seed is passed on).
        seed = extra[extra.index("--seed") + 1] if "--seed" in extra else "1"
        plan = [a for i, a in enumerate(extra) if a != "--seed" and (i == 0 or extra[i - 1] != "--seed")]
        subprocess.run([os.path.join(ca.ROOT, "build", "say"), texts[name], wave48, "--score", score, "--plan-only",
                        "--pitch", os.environ.get("PITCH", "120"), "--phones",
                        os.environ.get("PHONES", os.path.join(ca.ROOT, "..", "tract", "data", "phones_sheet.txt"))] + plan,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        subprocess.run([os.path.join(ca.ROOT, "..", "tract", "build", "sheet_speak"), score, wave48, "--seed", seed] +
                       os.environ.get("SHEET", "").split(),
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    elif os.environ.get("ENGINE") == "tract":
        # The segment-module engine speaks the planner's score at 750 Pa (its
        # breath; 600 fluttered); PHONES optionally replaces the phone table.
        seed = extra[extra.index("--seed") + 1] if "--seed" in extra else "1"
        plan = [a for i, a in enumerate(extra) if a != "--seed" and (i == 0 or extra[i - 1] != "--seed")]
        table = ["--phones", os.environ.get("PHONES", os.path.join(ca.ROOT, "..", "tract", "data", "phones_tract.txt"))]
        subprocess.run([os.path.join(ca.ROOT, "build", "say"), texts[name], wave48, "--score", score, "--plan-only",
                        "--pressure", "750"] + table + plan,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        # SPEAK_ENV="VAR=value ..." passes engine switches to speak.
        speak_env = dict(os.environ)
        for item in os.environ.get("SPEAK_ENV", "").split():
            if "=" in item:
                key, value = item.split("=", 1)
                speak_env[key] = value
        subprocess.run([os.path.join(ca.ROOT, "..", "tract", "build", "speak"), score, wave48, seed], env=speak_env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    else:
        subprocess.run([os.path.join(ca.ROOT, "build", "say"), texts[name], wave48, "--score", score] + extra,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    data, rate = ca.read_wave(wave48)
    wave16 = os.path.join(ca.OUT, name + "_ours.wav")
    ca.write_wave(wave16, ca.resample(data, rate, 16000), 16000)
    if analyse:
        subprocess.run([os.path.join(ca.ROOT, "build", "tf_analyze"), wave16, os.path.join(ca.OUT, name + "_ours")],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return name


def render_job(arguments):
    name, extra, analyse = arguments
    return render(name, extra, analyse)


def render_all(first, last, extra, analyse):
    os.makedirs(ca.OUT, exist_ok=True)
    jobs = [("arctic_a%04d" % n, extra, analyse) for n in range(first, last + 1)]
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as pool:
        list(pool.map(render_job, jobs))


if __name__ == "__main__":
    numbers = [int(a) for a in sys.argv[1:] if not a.startswith("--")]
    first = numbers[0] if numbers else 1
    last = numbers[1] if len(numbers) > 1 else 60
    render_all(first, last, os.environ.get("SAY_ARGS", "").split(), "--analyse" in sys.argv)
