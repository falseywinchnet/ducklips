#!/usr/bin/env python3
"""Audit the spatial-pressure alphabet and its point-sampled airway defect.

This measures existing native renders. It neither edits audio nor rates speech.
The pad/floor integration is a controlled geometry reference, not a replacement
for the full tongue, blade, anatomy residuals, or acoustic solver.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def measurements(path):
    audio, rate = sf.read(path, always_2d=True)
    if rate != 48000 or len(audio) < round(12.65 * rate) or not np.isfinite(audio).all():
        raise ValueError("expected complete finite 48 kHz alphabet")
    results = {}
    for name, start, end in [("foxtrot", 4.55, 5.65), ("november", 11.70, 12.65),
                             ("foxtrot_ks", 4.952, 4.970), ("foxtrot_later", 5.023, 5.029)]:
        samples = audio[round(start * rate):round(end * rate)]
        jumps = np.abs(np.diff(samples[:, 0]))
        results[name] = dict(start=start, end=end, peak=float(np.max(np.abs(samples))),
                             rms=float(np.sqrt(np.mean(samples ** 2))),
                             largest_left_jump=float(jumps.max()),
                             jump_time=start + (int(jumps.argmax()) + 1) / rate)
    return dict(sha256=digest(path), frames=len(audio), rate=rate, intervals=results)


def pad_reference(geometry):
    row = geometry[round(5.025, 6)]
    section = 58
    width = 17.5 / 64
    centre = (section + .5) * width - 7 - np.pi + 2
    tip_y, jaw, radius = row["tip_y"], row["jaw"], .35

    def area(x, tip_x):
        roof = 9 + np.interp(x, [2, 4.5, 6, 7.2, 8.4, 9.4], [.9, 1.1, .9, .55, .5, .6])
        floor = 9 - (.35 + jaw if x >= 8 else 1 + 1.4 * jaw)
        gap = roof - floor
        dx = x - tip_x
        if abs(dx) <= radius:
            gap = min(gap, max(roof - tip_y - np.sqrt(max(0, radius ** 2 - dx ** 2)), 0))
        return 1.9 * max(gap, 0) ** 1.45

    results = []
    for epsilon in [1e-2, 1e-3, 1e-4, 1e-5]:
        samples = []
        for shift in [epsilon, -epsilon]:
            tip_x = centre - radius + shift
            left, right = centre - width / 2, centre + width / 2
            points = [p for p in [tip_x - radius, tip_x + radius, 8] if left < p < right]
            volume, error = quad(lambda x: area(x, tip_x), left, right,
                                 points=points, epsabs=1e-11, limit=200)
            samples.append((area(centre, tip_x) * width, volume, error))
        results.append(dict(tip_displacement_cm=2 * epsilon,
                            midpoint_volume_jump_mL=samples[1][0] - samples[0][0],
                            integrated_volume_change_mL=samples[1][1] - samples[0][1],
                            quadrature_error_mL=max(s[2] for s in samples)))
    return dict(model="controlled pad and floor only; no body, blade, anatomy residual or smoothing",
                section=section, tip_y_cm=tip_y, jaw=jaw, rows=results)


def records(path):
    with path.open() as stream:
        return {round(float(row["time"]), 6):
                {key: float(value) for key, value in row.items() if key != "phone"}
                for row in csv.DictReader(stream)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, default=ROOT / "out/nato-closed-contact.wav")
    parser.add_argument("--after", type=Path, default=ROOT / "out/nato-spatial-pressure.wav")
    parser.add_argument("--geometry", type=Path, default=ROOT / "out/nato-spatial-geometry.csv")
    parser.add_argument("--airway", type=Path, default=ROOT / "out/nato-spatial-airway.csv")
    parser.add_argument("--probe", type=Path, default=ROOT / "out/nato-spatial-pressure.csv")
    parser.add_argument("--launch", type=Path, default=ROOT / "out/nato-spatial-pressure-launch.json")
    parser.add_argument("--out", type=Path, default=ROOT / "reports/spatial_pressure_reference.json")
    args = parser.parse_args()
    launch = json.loads(args.launch.read_text())
    for name, expected in launch["sha256"].items():
        if digest(ROOT / name) != expected:
            raise RuntimeError(f"render source/binary changed since launch: {name}")
    geometry, airway, probe = records(args.geometry), records(args.airway), records(args.probe)
    errors = []
    for time, row in geometry.items():
        sample = probe[round(time + .001, 6)]
        for a, b in [("root_x", "blade_root_x_cm"), ("root_y", "blade_root_y_cm"),
                     ("attachment_motion", "attachment_motion_cm"), ("miss", "tip_target_error_cm")]:
            errors.append(abs(row[a] - sample[b]))
    if max(errors) != 0:
        raise RuntimeError("geometry replay does not exactly match printed speech probe")
    events = []
    for time, section in [(5.025, 58), (5.071, 52), (5.073, 51)]:
        previous = round(time - .001, 6)
        ray_x = (section + .5) * 17.5 / 64 - 7 - np.pi + 2
        before, after = airway[previous], airway[time]
        key, length = f"area_{section}", f"length_{section}"
        events.append(dict(time=time, section=section, ray_x_cm=ray_x,
                           before_tip_distance_cm=ray_x - geometry[previous]["tip_x"],
                           after_tip_distance_cm=ray_x - geometry[time]["tip_x"], tip_radius_cm=.35,
                           before_area_cm2=before[key] * 1e4, after_area_cm2=after[key] * 1e4,
                           section_volume_jump_mL=(after[key] * after[length] - before[key] * before[length]) * 1e6,
                           total_volume_jump_mL=(after["oral_volume"] - before["oral_volume"]) * 1e6))
    result = dict(use="experimental physical integration; not an intelligibility or realism score",
                  launch=launch, before=measurements(args.before), after=measurements(args.after),
                  geometry_rows=len(geometry), maximum_printed_geometry_error=max(errors),
                  silhouette_crossings=events, controlled_pad_reference=pad_reference(geometry),
                  inputs={str(p.relative_to(ROOT)): digest(p) for p in
                          [args.geometry, args.airway, args.probe, args.launch]},
                  audit_source_sha256=digest(Path(__file__)))
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ["before", "after", "silhouette_crossings"]}, indent=2))


if __name__ == "__main__":
    main()
