"""A video of the midsagittal mouth as the engine speaks: speak's GEOLOG
geometry (the hull, the organs, the physical tongue's sections) drawn frame
by frame beside the vocal tube's areas, with the render's audio.

usage: mouth_video.py SCORE GEOLOG WAV OUT.mp4 [SLOW]
SLOW: playback slower than real time by this factor (the GEOLOG must then
hold frames every 20 / SLOW ms); the audio is stretched to match (its pitch
kept)."""
import os
import re
import shutil
import subprocess
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pyplot  # noqa: E402
import numpy  # noqa: E402
import soundfile  # noqa: E402

FPS = 50
NODE_LENGTH = 0.3646  # cm, the vocal tube's node spacing (tract/tube.hpp)
FFMPEG = "/opt/homebrew/bin/ffmpeg"
TISSUE = "#e9b8a8"
TONGUE = "#c75b5b"
BONE = "#d8d2c4"
AIR = "#ffffff"


def read_geometry(path):
    frames = []
    for line in open(path):
        fields = [part.split() for part in line.split("|")]
        time = float(fields[0][0])
        hull = numpy.array(fields[1], dtype=float).reshape(-1, 2)
        organs = numpy.array(fields[2], dtype=float).reshape(-1, 2)
        tongue_points = int(fields[3][0])
        underside_points = int(fields[3][1]) if len(fields[3]) > 1 else 0
        areas = numpy.array(fields[-1], dtype=float)
        guard = None
        penetration = numpy.array(fields[7], dtype=float) if len(fields) >= 9 and fields[7] else None
        # The skeleton tongue's two supports (base, tip), each kBoneSegments + 1 points.
        bones = numpy.array(fields[8], dtype=float).reshape(2, -1, 2) if len(fields) == 10 and fields[8] else None
        if len(fields) >= 8 and fields[4]:
            # The airway guard: the planned surface, the guard surface, the intrusions.
            guard = (numpy.array(fields[4], dtype=float).reshape(-1, 2), numpy.array(fields[5], dtype=float).reshape(-1, 2),
                     numpy.array(fields[6], dtype=float))
        frames.append((time, hull, organs, tongue_points, areas, guard, penetration, underside_points, bones))
    return frames


def read_phones(score):
    phones = []
    for line in open(score):
        found = re.match(r"# (\S+) ([\d.]+) ([\d.]+)", line)
        if found:
            phones.append((found.group(1).split("-")[0].rstrip("012").lower(), float(found.group(2)),
                           float(found.group(3))))
    return phones


def draw_frame(axes_mouth, axes_tube, axes_wave, frame, phones, signal, rate, duration, largest_nodes):
    time, hull, organs, tongue_points, areas, guard, penetration, underside_points, bones = frame
    axes_mouth.clear()
    axes_mouth.set_facecolor(AIR)
    # The skull's side of the hull: the hull closed through the top of the frame.
    skull = numpy.vstack([hull, [[hull[-1, 0], 13.0], [-1.5, 13.0], [-1.5, hull[0, 1]]]])
    axes_mouth.fill(skull[:, 0], skull[:, 1], color=BONE, zorder=1)
    axes_mouth.plot(hull[:, 0], hull[:, 1], color="#6b6358", linewidth=1.5, zorder=2)
    axes_mouth.fill(organs[:, 0], organs[:, 1], color=TISSUE, zorder=3)
    axes_mouth.plot(numpy.append(organs[:, 0], organs[0, 0]), numpy.append(organs[:, 1], organs[0, 1]),
                    color="#8a5a4a", linewidth=1.0, zorder=4)
    if underside_points > 0:
        # A thick tongue: its whole body (surface, cap, underside, frenulum)
        # drawn over the floor, the air beneath its blade left white.
        body = organs[2:2 + tongue_points + underside_points]
        floor_line = organs[2 + tongue_points + underside_points - 1:2 + tongue_points + underside_points + 3]
        axes_mouth.fill(body[:, 0], body[:, 1], color="#d98a8a", zorder=4)
        axes_mouth.plot(numpy.append(body[:, 0], body[0, 0]), numpy.append(body[:, 1], body[0, 1]), color="#7a1f1f",
                        linewidth=1.0, zorder=4)
    if tongue_points > 0:
        # The physical tongue: the root's base, then its sections to the tip.
        surface = organs[2:3 + tongue_points - 1]
        axes_mouth.plot(surface[:, 0], surface[:, 1], color=TONGUE, linewidth=2.5, zorder=5)
        colours = ["#7a1f1f"] * (len(surface) - 1)
        if guard is not None:
            planned, wall, intrusion = guard
            # The guarded airway (the inverse tongue): between the guard
            # surface and the hull, shaded; the planned surface dashed.
            axes_mouth.fill(numpy.append(wall[1:, 0], planned[::-1][:-1, 0]),
                            numpy.append(wall[1:, 1], planned[::-1][:-1, 1]), color="#cfe3f7", alpha=0.0, zorder=2)
            axes_mouth.plot(wall[1:, 0], wall[1:, 1], color="#2e86de", linewidth=1.5, linestyle="-", zorder=7)
            axes_mouth.plot(planned[1:, 0], planned[1:, 1], color="#555555", linewidth=1.0, linestyle="--", zorder=7)
            colours = ["#f1c40f" if value > 0.02 else "#7a1f1f" for value in intrusion[1:]]
            worst = float(intrusion.max())
            axes_mouth.text(-1.2, -0.6, "guard (blue): deepest intrusion %.2f cm" % worst, fontsize=11,
                            color="#c0392b" if worst > 0.05 else "#2e86de")
        if penetration is not None:
            # Sections inside the hull (the palate, ridge or teeth): magenta.
            colours = ["#e100ff" if depth > 0.05 else colour for depth, colour in zip(penetration[1:], colours)]
            deepest = float(penetration.max())
            axes_mouth.text(-1.2, 0.0, "inside the palate: deepest %.2f cm" % deepest, fontsize=11,
                            color="#e100ff" if deepest > 0.05 else "#2e86de")
        axes_mouth.scatter(surface[1:, 0], surface[1:, 1], s=14, color=colours, zorder=8)
    if bones is not None:
        # The two supports inside the tongue: the base (red), the tip (green).
        for points, colour in zip(bones, ("#d62728", "#2ca02c")):
            axes_mouth.plot(points[:, 0], points[:, 1], color=colour, linewidth=4.0, solid_capstyle="round", zorder=9)
            axes_mouth.scatter(points[:, 0], points[:, 1], s=18, color="white", edgecolors=colour, zorder=10)
    axes_mouth.set_xlim(-1.5, 11.0)
    axes_mouth.set_ylim(-1.0, 12.0)
    axes_mouth.set_aspect("equal")
    axes_mouth.set_xticks([])
    axes_mouth.set_yticks([])
    sounding = [name for name, start, length in phones if start <= time < start + length]
    axes_mouth.set_title("/%s/   %.2f s" % (sounding[0] if sounding else "", time), fontsize=16, loc="left")
    axes_tube.clear()
    places = numpy.arange(len(areas)) * NODE_LENGTH
    axes_tube.bar(places, areas, width=NODE_LENGTH * 0.95, color=["#c0392b" if a < 0.05 else "#4a6fa5" for a in areas])
    axes_tube.set_xlim(0, largest_nodes * NODE_LENGTH)
    axes_tube.set_ylim(0, 6.0)
    axes_tube.set_xlabel("distance from the glottis (cm)")
    axes_tube.set_ylabel("area (cm²)")
    axes_tube.set_title("the tube the sound passes through (red: < 0.05 cm²)", fontsize=10, loc="left")
    axes_wave.clear()
    window = int(0.25 * rate)
    middle = int(time * rate)
    lo = max(0, middle - window)
    hi = min(len(signal), middle + window)
    axes_wave.plot((numpy.arange(lo, hi) / rate), signal[lo:hi], color="#333333", linewidth=0.6)
    axes_wave.axvline(time, color="#c0392b")
    axes_wave.set_xlim(time - 0.25, time + 0.25)
    axes_wave.set_ylim(-1.0, 1.0)
    axes_wave.set_yticks([])
    for name, start, length in phones:
        if start + length > time - 0.25 and start < time + 0.25:
            axes_wave.text(start + length / 2, 0.8, name, ha="center", fontsize=10)
            axes_wave.axvline(start, color="#bbbbbb", linewidth=0.5)


def main():
    score, geometry, wave, out = sys.argv[1:5]
    slow = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    frames = read_geometry(geometry)
    phones = read_phones(score)
    signal, rate = soundfile.read(wave)
    if signal.ndim > 1:
        signal = signal[:, 0]
    signal = signal / (numpy.abs(signal).max() + 1e-9)
    largest_nodes = max(len(frame[4]) for frame in frames)
    intrusions = [frame[5][2] for frame in frames if frame[5] is not None]
    if intrusions:
        deepest = numpy.array([values.max() for values in intrusions])
        print("intrusion into the guarded airway: median %.3f cm, 90th %.3f, largest %.3f; frames over 0.05: %d of %d" % (
            numpy.median(deepest), numpy.percentile(deepest, 90), deepest.max(), int((deepest > 0.05).sum()),
            len(deepest)))
    folder = tempfile.mkdtemp(prefix="mouth_")
    figure = pyplot.figure(figsize=(12.8, 7.2), dpi=100)
    axes_mouth = figure.add_axes([0.02, 0.05, 0.5, 0.88])
    axes_tube = figure.add_axes([0.58, 0.42, 0.40, 0.48])
    axes_wave = figure.add_axes([0.58, 0.08, 0.40, 0.22])
    for index, frame in enumerate(frames):
        draw_frame(axes_mouth, axes_tube, axes_wave, frame, phones, signal, rate, len(signal) / rate, largest_nodes)
        figure.savefig(os.path.join(folder, "frame_%05d.png" % index))
    audio = wave
    if slow != 1.0:
        audio = os.path.join(folder, "slow.wav")
        tempo = 1.0 / slow
        filters = []
        while tempo < 0.5:
            filters.append("atempo=0.5")
            tempo /= 0.5
        filters.append("atempo=%.4f" % tempo)
        subprocess.run([FFMPEG, "-y", "-loglevel", "error", "-i", wave, "-filter:a", ",".join(filters), audio],
                       check=True)
    subprocess.run([FFMPEG, "-y", "-loglevel", "error", "-framerate", str(FPS), "-i",
                    os.path.join(folder, "frame_%05d.png"), "-i", audio, "-c:v", "libx264", "-pix_fmt", "yuv420p",
                    "-c:a", "aac", "-shortest", out], check=True)
    shutil.rmtree(folder)
    print("wrote", out, len(frames), "frames")


if __name__ == "__main__":
    main()
