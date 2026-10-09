#!/usr/bin/env python3
"""Controlled pressure-topology experiment using a native speech geometry replay.

This does not render speech. All four runs receive the same prescribed larynx
and passage data; only front-volume motion and chamber topology differ.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(path):
    with path.open() as stream:
        return [{k: float(v) if k != "phone" else v for k, v in row.items()}
                for row in csv.DictReader(stream)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, default=ROOT / "build")
    parser.add_argument("--body", type=Path, default=ROOT / "bodies/reference.body")
    parser.add_argument("--dictionary", type=Path, default=ROOT / "samples/nato.dict")
    parser.add_argument("--text", default="Foxtrot.")
    parser.add_argument("--probe", required=True, type=Path, help="actual speech block probe")
    parser.add_argument("--trace", required=True, type=Path, help="actual per-sample source trace")
    parser.add_argument("--out", type=Path, default=ROOT / "out/cavity-partition")
    parser.add_argument("--partition", type=int, default=42, help="fixed geometry audit cut, zero based")
    parser.add_argument("--start", type=float, default=.870)
    parser.add_argument("--end", type=float, default=.880)
    args = parser.parse_args()
    if not (math.isfinite(args.start) and math.isfinite(args.end) and args.end > args.start):
        parser.error("end must follow start")
    args.out.mkdir(parents=True, exist_ok=True)
    geometry_file = args.out / "geometry.csv"
    command = [str(args.build / "ducklips_tongue_reach_audit"), str(args.body),
               str(args.dictionary), args.text, "--partition", str(args.partition)]
    with geometry_file.open("w") as stream:
        subprocess.run(command, stdout=stream, check=True)
    geometry = records(geometry_file)
    probes = {round(row["time"], 6): row for row in records(args.probe)}
    source = records(args.trace)
    dt = geometry[1]["time"] - geometry[0]["time"]
    match_fields = {"root_x": "blade_root_x_cm", "root_y": "blade_root_y_cm",
                    "attachment_motion": "attachment_motion_cm", "miss": "tip_target_error_cm"}
    maximum_geometry_difference = 0.0
    for row in geometry:
        probe = probes[round(row["time"] + dt, 6)]
        for geometry_key, probe_key in match_fields.items():
            maximum_geometry_difference = max(maximum_geometry_difference,
                                              abs(row[geometry_key] - probe[probe_key]))
    if maximum_geometry_difference != 0:
        raise RuntimeError("geometry replay does not match supplied speech probe")
    selected = [row for row in geometry if args.start <= row["time"] <= args.end]
    if len(selected) < 2 or any(row["neck_area"] >= 1e-12 for row in selected):
        raise RuntimeError("choose an interval with a sealed partition in every geometry row")
    # Native AirProperties defaults. These are a controlled reference, not fitted anatomy.
    density, sound_speed, viscosity = 1.14, 350., 1.86e-5
    wall_compliance = 1e-8
    initial_pressure = [600., 0.]
    results = {}
    for motion in ("moving", "frozen"):
        input_file = args.out / f"{motion}-input.csv"
        with input_file.open("w") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            for row in selected:
                # Each geometry is the control-block target. Use the matching
                # block's actual mean glottal area, not the block-lagged probe.
                time = row["time"]
                samples = [r for r in source if time <= r["time"] < round(time + dt, 6)]
                if len(samples) != 48:
                    raise RuntimeError(f"expected 48 source samples at {time}, got {len(samples)}")
                glottis = sum(r["area"] for r in samples) / len(samples)
                lung = probes[round(time + dt, 6)]["lung_pa"]
                rear = row["upstream_volume"]
                front = row["downstream_volume"] if motion == "moving" else selected[0]["downstream_volume"]
                areas = [glottis, row["neck_area"], row["front_area"]]
                areas = [a if a >= 1e-12 else 0. for a in areas]
                # Circular Poiseuille resistance is only an explicit port-law
                # assumption for this experiment. The front passage length is
                # the same one-cell proxy as the neck, not measured anatomy.
                lengths = [.003, row["neck_length"], row["neck_length"]]
                resistance = [8 * math.pi * viscosity * length / area**2 if area else 0.
                              for area, length in zip(areas, lengths)]
                caps = [v / (density * sound_speed**2) + wall_compliance / 2 for v in (rear, front)]
                writer.writerow([time, rear, front, *caps, *areas, *resistance, lung])
        for topology in ("separate", "lumped"):
            output_file = args.out / f"{motion}-{topology}.csv"
            command = [str(args.build / "ducklips_cavity_network_probe"), *map(str, initial_pressure)]
            if topology == "lumped":
                command.append("--lumped")
            with input_file.open() as stream, output_file.open("w") as output:
                subprocess.run(command, stdin=stream, stdout=output, check=True)
            results[f"{motion}_{topology}"] = records(output_file)
    separate_difference = max(abs(a["rear_pa"] - b["rear_pa"]) for a, b in
                              zip(results["moving_separate"], results["frozen_separate"]))
    bridge_flow = max(abs(row["neck_flow"]) for row in results["moving_separate"])
    if separate_difference >= 1e-6 or bridge_flow != 0:
        raise RuntimeError("sealed partition did not isolate front motion")
    interval = [r for r in source if args.start <= r["time"] <= args.end]
    receipt = {
        "scope": "prescribed two-chamber pressure reference; not production speech integration",
        "text": args.text, "interval_seconds": [args.start, args.end],
        "fixed_partition_section": args.partition, "control_dt_seconds": dt,
        "air": {"density_kg_m3": density, "sound_speed_m_s": sound_speed, "viscosity_pa_s": viscosity},
        "wall_compliance_m3_pa": wall_compliance, "wall_partition": "equal halves, illustrative",
        "initial_pressure_pa": initial_pressure,
        "geometry_replay_maximum_difference_cm": maximum_geometry_difference,
        "geometry_replay_compared_steps": len(geometry),
        "rear_volume_change_m3": selected[-1]["upstream_volume"] - selected[0]["upstream_volume"],
        "front_volume_change_m3": selected[-1]["downstream_volume"] - selected[0]["downstream_volume"],
        "separate_rear_motion_effect_maximum_pa": separate_difference,
        "sealed_neck_maximum_flow_m3_s": bridge_flow,
        "lumped_rear_motion_effect_final_pa": results["moving_lumped"][-1]["rear_pa"] - results["frozen_lumped"][-1]["rear_pa"],
        "baseline_production_interval_peak_pressure_pa": max(r["cavity_pressure"] for r in interval),
        "final": {key: rows[-1] for key, rows in results.items()},
        "maximum_mass_residual_m3_s": max(r["mass_residual"] for rows in results.values() for r in rows),
        "inputs": {str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path): sha256(path)
                   for path in (args.body.resolve(), args.dictionary.resolve(), args.probe.resolve(), args.trace.resolve())},
        "native_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in
                          (ROOT / "include/ducklips/physics.hpp", ROOT / "src/cavity_network.cpp",
                           ROOT / "include/ducklips/cavity_network.hpp", ROOT / "tools/cavity_network_probe.cpp",
                           ROOT / "tools/tongue_reach_audit.cpp")},
        "binary_sha256": {name: sha256(args.build / name) for name in
                          ("ducklips_cavity_network_probe", "ducklips_tongue_reach_audit")},
        "artifacts": {path.name: sha256(path) for path in sorted(args.out.glob("*.csv"))},
    }
    (args.out / "reference.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({key: receipt[key] for key in ("rear_volume_change_m3", "front_volume_change_m3",
                      "separate_rear_motion_effect_maximum_pa", "lumped_rear_motion_effect_final_pa",
                      "maximum_mass_residual_m3_s")}, indent=2))


if __name__ == "__main__":
    main()
