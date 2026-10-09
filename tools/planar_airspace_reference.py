#!/usr/bin/env python3
"""Independent polygon audit of native finite-volume air-space geometry.

Requires the locally available Shapely/GEOS for independent Boolean geometry.
No audio is synthesized, altered, recognized, or scored by this tool.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import shapely
from shapely.geometry import Polygon, box, LineString
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parts(geometry):
    if geometry.geom_type == "Polygon":
        return [geometry]
    return [piece for piece in getattr(geometry, "geoms", []) if piece.geom_type == "Polygon"]


def native_components(p):
    nx, ny = len(p["x_edges"]) - 1, len(p["y_edges"]) - 1
    seen, volumes = set(), []
    for seed, volume in enumerate(p["volume"]):
        if volume <= 0 or seed in seen:
            continue
        pending, total = [seed], 0.
        seen.add(seed)
        while pending:
            node = pending.pop()
            row, column = divmod(node, nx)
            total += p["volume"][node]
            links = []
            if column:
                links.append((node - 1, p["vertical_face"][row * (nx + 1) + column]))
            if column + 1 < nx:
                links.append((node + 1, p["vertical_face"][row * (nx + 1) + column + 1]))
            if row:
                links.append((node - nx, p["horizontal_face"][row * nx + column]))
            if row + 1 < ny:
                links.append((node + nx, p["horizontal_face"][(row + 1) * nx + column]))
            for other, aperture in links:
                if aperture > 0 and p["volume"][other] > 0 and other not in seen:
                    seen.add(other)
                    pending.append(other)
        volumes.append(total)
    return sorted(volumes)


def draw_snapshot(snapshot, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as PlotPath
    from matplotlib.patches import PathPatch, Patch
    from shapely.geometry.polygon import orient
    p = json.loads(snapshot.read_text())
    domain = Polygon(p["domain"])
    fluid = domain.difference(unary_union([Polygon(solid) for solid in p["solids"]]))
    fig, axis = plt.subplots(figsize=(8, 8.6))

    def patch(polygon, color):
        polygon = orient(polygon)
        rings = [polygon.exterior, *polygon.interiors]
        paths = [PlotPath(np.asarray(ring.coords) * 100, closed=True) for ring in rings]
        axis.add_patch(PathPatch(PlotPath.make_compound_path(*paths), facecolor=color, edgecolor="#263c46", linewidth=.8))

    patch(domain, "#df9d61")
    ordered = sorted(parts(fluid), key=lambda part: part.centroid.x)
    colors = ["#167fa0", "#43b5ac"]
    labels = ["Rear air", "Front and sublingual air"]
    legend = [Patch(color="#df9d61", label="Tissue inside the tract")]
    for index, region in enumerate(ordered):
        patch(region, colors[index])
        legend.append(Patch(color=colors[index], label=f"{labels[index]}: {region.area * p['depth'] * 1e6:.3f} mL"))
    axis.set(xlim=(-1, 10), ylim=(0, 10.5), xlabel="Forward position (cm)", ylabel="Height (cm)",
             title="Foxtrot /t/ closure — 5.040 s\nTwo separate air chambers in the solid reference")
    axis.set_aspect("equal")
    axis.legend(handles=legend, loc="lower right", frameon=False, fontsize=10)
    fig.text(.5, .025, "Geometry reference: 2.2 cm extrusion and an unvalidated blade model.\nThis is not a new speech rendering.", ha="center", fontsize=10)
    fig.tight_layout(rect=(0, .07, 1, 1))
    fig.savefig(output, dpi=150)
    plt.close(fig)


def audit(snapshot):
    p = json.loads(snapshot.read_text())
    domain = Polygon(p["domain"])
    solids = [Polygon(contour) for contour in p["solids"]]
    if not domain.is_valid or not all(solid.is_valid for solid in solids):
        raise RuntimeError("invalid native polygon; no automatic geometry repair allowed")
    fluid = domain.difference(unary_union(solids))
    x, y, depth = p["x_edges"], p["y_edges"], p["depth"]
    nx, ny = len(x) - 1, len(y) - 1
    volume_error, vertical_error, horizontal_error = 0., 0., 0.
    split = []
    for row in range(ny):
        for column in range(nx):
            clipped = fluid.intersection(box(x[column], y[row], x[column + 1], y[row + 1]))
            volume_error = max(volume_error, abs(clipped.area * depth - p["volume"][row * nx + column]))
            positive = [piece for piece in parts(clipped) if piece.area > 1e-18]
            if len(positive) > 1:
                split.append(dict(column=column, row=row, pieces=len(positive)))
    for row in range(ny):
        for face in range(1, nx):
            line = LineString([(x[face], y[row]), (x[face], y[row + 1])])
            vertical_error = max(vertical_error, abs(fluid.intersection(line).length * depth - p["vertical_face"][row * (nx + 1) + face]))
    for face in range(1, ny):
        for column in range(nx):
            line = LineString([(x[column], y[face]), (x[column + 1], y[face])])
            horizontal_error = max(horizontal_error, abs(fluid.intersection(line).length * depth - p["horizontal_face"][face * nx + column]))
    total_error = abs(fluid.area * depth - sum(p["volume"]))
    if volume_error > 1e-15 or total_error > 1e-14 or max(vertical_error, horizontal_error) > 1e-13:
        raise RuntimeError("native cells/faces disagree with independent Boolean geometry")
    native = native_components(p)
    exact = sorted(piece.area * depth for piece in parts(fluid))
    if split or len(native) != len(exact) or any(abs(a - b) > 1e-14 for a, b in zip(native, exact)):
        raise RuntimeError("one-pressure-per-cell graph does not preserve exact fluid connectivity")
    return dict(time=p["time"], maximum_cell_volume_error_m3=volume_error,
                total_volume_error_m3=total_error, maximum_vertical_face_error_m2=vertical_error,
                maximum_horizontal_face_error_m2=horizontal_error, split_pressure_cells=split,
                connected_air_volumes_mL=[piece.area * depth * 1e6 for piece in parts(fluid)],
                native_component_volumes_mL=[value * 1e6 for value in native],
                snapshot=str(snapshot.relative_to(ROOT)), sha256=digest(snapshot))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, default=ROOT / "out/nato-airspace-geometry.csv")
    parser.add_argument("--out", type=Path, default=ROOT / "out/planar-airspace")
    parser.add_argument("--receipt", type=Path, default=ROOT / "reports/planar_airspace_reference.json")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    # The added root-posture column must not change any existing replay value.
    with args.geometry.open() as f, (ROOT / "out/nato-spatial-geometry.csv").open() as old:
        new_rows, old_rows = list(csv.DictReader(f)), list(csv.DictReader(old))
    if len(new_rows) != len(old_rows) or any(a[key] != b[key] for a, b in zip(new_rows, old_rows) for key in b):
        raise RuntimeError("native geometry replay changed since the assessed rendering")
    binary = ROOT / "build/ducklips_airspace_geometry_probe"
    sweeps = {}
    commands = []
    for count in [128, 256, 512, 1024]:
        output = args.out / f"foxtrot-{count}.csv"
        command = [str(binary), str(args.geometry), "5.020", "5.085", str(count), str(output)]
        subprocess.run(command, check=True)
        commands.append(command)
        data = np.genfromtxt(output, delimiter=",", names=True)
        sweeps[count] = data
    details = {}
    for count, data in sweeps.items():
        events = {}
        for time in [5.025, 5.071, 5.073]:
            index = int(np.argmin(abs(data["time"] - time)))
            events[str(time)] = float((data["air_volume_m3"][index] - data["air_volume_m3"][index - 1]) * 1e6)
        details[count] = dict(rows=len(data), maximum_reconstructed_pose_error_cm=float(data["maximum_pose_error_cm"].max()),
                              largest_step_mL=float(abs(np.diff(data["air_volume_m3"])).max() * 1e6),
                              event_changes_mL=events,
                              difference_from_1024_mL=float(abs(data["air_volume_m3"] - sweeps[1024]["air_volume_m3"]).max() * 1e6))
    snapshots = []
    for time in [5.025, 5.040, 5.071, 5.073]:
        prefix = args.out / f"pose-{time:.3f}"
        command = [str(binary), str(args.geometry), str(time), str(time), "512", str(prefix) + ".csv", str(prefix) + ".json"]
        subprocess.run(command, check=True)
        commands.append(command)
        snapshots.append(audit(Path(str(prefix) + ".json")))
    closure = next(item for item in snapshots if item["time"] == 5.040)
    if len(closure["connected_air_volumes_mL"]) != 2:
        raise RuntimeError("diagnostic solid /t/ does not separate front and rear air")
    source = [ROOT / "include/ducklips/planar_airspace.hpp", ROOT / "src/planar_airspace.cpp",
              ROOT / "tools/airspace_geometry_probe.cpp", ROOT / "tools/tongue_reach_audit.cpp",
              ROOT / "tests/planar_airspace_tests.cpp", Path(__file__), binary, args.geometry]
    result = dict(use="extruded planar geometry reference, not newly rendered speech or validated anatomy",
                  assumptions=dict(depth_m=.022, blade="finite lower body-material attachment; tapered to pad", grid=[64, 64],
                                   body="native conservative material map", tongue_mass="combined blade/body union volume is not constrained"),
                  geometry_replay_rows=len(new_rows), existing_columns_identical=True,
                  shapely_version=shapely.__version__, geos_version=shapely.geos_version_string,
                  contour_refinement=details, independent_snapshots=snapshots, commands=commands,
                  hashes={str(path.relative_to(ROOT)): digest(path) for path in source})
    figure = args.out / "foxtrot-closure.png"
    draw_snapshot(args.out / "pose-5.040.json", figure)
    result["figure"] = dict(path=str(figure.relative_to(ROOT)), sha256=digest(figure))
    args.receipt.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(dict(refinement=details, snapshots=snapshots), indent=2))


if __name__ == "__main__":
    main()
