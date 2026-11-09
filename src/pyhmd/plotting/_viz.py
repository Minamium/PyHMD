#!/usr/bin/env python3
"""Scientific plots, ParaView export, animations, and decomposition comparisons.

All magnetic quantities use the solver's mu_0=1 convention.  Input data are
cell-centred CSV snapshots; coordinates and cell faces come from metadata.json.
This tool intentionally does not interpolate MPI pieces or conceal missing cells.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass

import numpy as np


FIELDS = ("rho", "vr", "vtheta", "vphi", "pressure", "Br", "Btheta", "Bphi", "psi", "divB")


@dataclass
class Snapshot:
    run_dir: Path
    path: Path
    metadata: dict
    time: float
    step: int
    data: dict[str, np.ndarray]

    @property
    def shape(self) -> tuple[int, int]:
        return int(self.metadata["nx2"]), int(self.metadata["nx1"])


def snapshot_files(run_dir: str | Path) -> list[Path]:
    paths = sorted(Path(run_dir).glob("snapshot_*.csv"))
    if not paths:
        raise ValueError(f"No snapshot_*.csv files found in {run_dir}")
    return paths


def _read_table(path: Path) -> tuple[np.ndarray, str]:
    comments: list[str] = []
    header_line = 0
    with path.open() as handle:
        for index, line in enumerate(handle):
            if not line.strip() or line.lstrip().startswith("#"):
                comments.append(line)
                continue
            header_line = index
            break
        else:
            raise ValueError(f"Empty data table: {path}")
    table = np.genfromtxt(path, delimiter=",", names=True, skip_header=header_line, encoding="utf-8")
    table = np.atleast_1d(table)
    if not table.dtype.names:
        raise ValueError(f"Missing CSV header: {path}")
    return table, " ".join(comments)


def load_snapshot(run_dir: str | Path, selection: str | Path = "latest") -> Snapshot:
    root = Path(run_dir).resolve()
    with (root / "metadata.json").open() as handle:
        metadata = json.load(handle)
    if metadata.get("geometry") not in ("spherical", "cartesian"):
        raise ValueError("metadata.geometry must be spherical or cartesian")
    ny, nx = int(metadata["nx2"]), int(metadata["nx1"])
    if nx < 1 or ny < 1:
        raise ValueError("Grid dimensions must be positive")
    selection = str(selection)
    if selection in ("latest", "first"):
        files = snapshot_files(root)
        path = files[-1] if selection == "latest" else files[0]
    else:
        path = Path(selection)
        if not path.is_absolute():
            path = root / path
    table, comment = _read_table(path)
    required = ("x1", "x2") + FIELDS
    missing = set(required) - set(table.dtype.names or ())
    if missing:
        raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
    if len(table) != nx * ny:
        raise ValueError(f"{path.name}: expected {nx * ny} cells, found {len(table)}")
    data = {field: np.asarray(table[field], dtype=np.float64).reshape(ny, nx) for field in required}
    for field in ("entropy", "total_energy", "internal_energy", "eta_p", "mu_poloidal"):
        if field in table.dtype.names:
            data[field] = np.asarray(table[field], dtype=np.float64).reshape(ny, nx)
    if not all(np.isfinite(values).all() for values in data.values()):
        raise ValueError(f"Non-finite values in {path}")
    if not np.allclose(data["x1"], data["x1"][0:1, :], rtol=1e-12, atol=1e-14):
        raise ValueError("x1 coordinates do not match a complete j-major structured grid")
    if not np.allclose(data["x2"], data["x2"][:, 0:1], rtol=1e-12, atol=1e-14):
        raise ValueError("x2 coordinates do not match a complete j-major structured grid")
    for key, count in (("x1_faces", nx + 1), ("x2_faces", ny + 1)):
        faces = np.asarray(metadata[key], dtype=float)
        if len(faces) != count or not np.isfinite(faces).all() or np.any(np.diff(faces) <= 0):
            raise ValueError(f"Invalid {key} in metadata")
    faces1, faces2 = np.asarray(metadata["x1_faces"]), np.asarray(metadata["x2_faces"])
    if not np.all((data["x1"] > faces1[:-1][None, :]) & (data["x1"] < faces1[1:][None, :])):
        raise ValueError("x1 cell centres lie outside their declared faces")
    if not np.all((data["x2"] > faces2[:-1][:, None]) & (data["x2"] < faces2[1:][:, None])):
        raise ValueError("x2 cell centres lie outside their declared faces")
    if metadata["geometry"] == "spherical" and (faces1[0] <= 0 or faces2[0] < 0 or faces2[-1] > np.pi):
        raise ValueError("Spherical faces require r>0 and 0<=theta<=pi")
    time_match = re.search(r"\btime\s*=\s*([+\-\d.eE]+)", comment)
    step_match = re.search(r"\bstep\s*=\s*(\d+)", comment)
    if time_match is None:
        raise ValueError(f"{path.name}: missing '# time=... step=...' comment")
    return Snapshot(root, path, metadata, float(time_match.group(1)),
                    int(step_match.group(1)) if step_match else -1, data)


def cell_volumes(snapshot: Snapshot) -> np.ndarray:
    """Return exact 3D spherical volumes per radian phi, or Cartesian areas."""
    faces1 = np.asarray(snapshot.metadata["x1_faces"], dtype=float)
    faces2 = np.asarray(snapshot.metadata["x2_faces"], dtype=float)
    if snapshot.metadata["geometry"] == "spherical":
        return (np.cos(faces2[:-1]) - np.cos(faces2[1:]))[:, None] * (np.diff(faces1**3) / 3)[None, :]
    return np.diff(faces2)[:, None] * np.diff(faces1)[None, :]


def _coordinates(snapshot: Snapshot, faces: bool = False) -> tuple[np.ndarray, np.ndarray]:
    if faces:
        first, second = np.meshgrid(snapshot.metadata["x1_faces"], snapshot.metadata["x2_faces"])
    else:
        first, second = snapshot.data["x1"], snapshot.data["x2"]
    if snapshot.metadata["geometry"] == "spherical":
        return first * np.sin(second), first * np.cos(second)
    return first, second


def _vector(snapshot: Snapshot, prefix: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if prefix == "v":
        radial, polar, azimuthal = (snapshot.data[name] for name in ("vr", "vtheta", "vphi"))
    else:
        radial, polar, azimuthal = (snapshot.data[name] for name in ("Br", "Btheta", "Bphi"))
    if snapshot.metadata["geometry"] == "spherical":
        theta = snapshot.data["x2"]
        return radial * np.sin(theta) + polar * np.cos(theta), azimuthal, radial * np.cos(theta) - polar * np.sin(theta)
    return radial, polar, azimuthal


def write_vtk(snapshot: Snapshot, output: str | Path) -> Path:
    """Write one 2D structured grid with cell data (ParaView-compatible ASCII)."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    ny, nx = snapshot.shape
    x, vertical = _coordinates(snapshot, faces=True)
    with output.open("w") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write(f"PyHMD {snapshot.metadata.get('case', '')} time={snapshot.time:.17g}\nASCII\n")
        handle.write(f"DATASET STRUCTURED_GRID\nDIMENSIONS {nx + 1} {ny + 1} 1\n")
        handle.write(f"POINTS {(nx + 1) * (ny + 1)} double\n")
        for xx, zz in zip(x.ravel(), vertical.ravel()):
            point = (xx, 0.0, zz) if snapshot.metadata["geometry"] == "spherical" else (xx, zz, 0.0)
            handle.write(" ".join(f"{value:.17g}" for value in point) + "\n")
        handle.write(f"CELL_DATA {nx * ny}\n")
        for field in FIELDS:
            handle.write(f"SCALARS {field} double 1\nLOOKUP_TABLE default\n")
            np.savetxt(handle, snapshot.data[field].ravel(), fmt="%.17g")
        beta = np.divide(2 * snapshot.data["pressure"], sum(snapshot.data[key]**2 for key in ("Br", "Btheta", "Bphi")),
                         out=np.full(snapshot.shape, np.nan),
                         where=sum(snapshot.data[key]**2 for key in ("Br", "Btheta", "Bphi")) > 0)
        handle.write("SCALARS plasma_beta double 1\nLOOKUP_TABLE default\n")
        np.savetxt(handle, beta.ravel(), fmt="%.17g")
        for name, prefix in (("velocity", "v"), ("magnetic_field", "B")):
            handle.write(f"VECTORS {name} double\n")
            np.savetxt(handle, np.column_stack([component.ravel() for component in _vector(snapshot, prefix)]), fmt="%.17g")
    return output


def _plot_imports():
    cache = Path(tempfile.gettempdir()) / "pymhd-matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titlesize": 11, "axes.labelsize": 10,
                         "figure.facecolor": "white", "savefig.facecolor": "white"})
    return plt, LogNorm, Normalize


def _log_limits(values: np.ndarray) -> tuple[float, float] | None:
    positive = values[np.isfinite(values) & (values > 0)]
    if not len(positive):
        return None
    low, high = float(positive.min()), float(positive.max())
    if low == high:
        return low / 1.1, high * 1.1
    return low, high


def render_snapshot(snapshot: Snapshot, output: str | Path, limits: dict | None = None,
                    window: dict | None = None) -> Path:
    plt, LogNorm, Normalize = _plot_imports()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    x, vertical = _coordinates(snapshot, faces=True)
    magnetic_squared = sum(snapshot.data[key] ** 2 for key in ("Br", "Btheta", "Bphi"))
    beta = np.divide(2 * snapshot.data["pressure"], magnetic_squared,
                     out=np.full(snapshot.shape, np.nan), where=magnetic_squared > 0)
    panels = [("rho", snapshot.data["rho"], "Density", "magma", True),
              ("pressure", snapshot.data["pressure"], "Gas pressure", "viridis", True),
              ("vphi", snapshot.data["vphi"], r"Azimuthal velocity $v_\phi$" if snapshot.metadata["geometry"] == "spherical" else r"Velocity $v_3$", "RdBu_r", False),
              ("beta", beta, r"Plasma $\beta=2p/B^2$", "cividis", True)]
    figure_height = 8.8
    if window and "xlim" in window and "ylim" in window:
        view_aspect = (window["ylim"][1] - window["ylim"][0]) / (window["xlim"][1] - window["xlim"][0])
        figure_height = min(11.8, max(5.8, 3.4 + 6.8 * view_aspect))
    figure, axes = plt.subplots(2, 2, figsize=(11.8, figure_height), constrained_layout=True)
    for axis, (key, values, title, cmap, logarithmic) in zip(axes.ravel(), panels):
        if logarithmic:
            domain = (limits or {}).get(key) or _log_limits(values)
            norm = LogNorm(*domain) if domain else None
            masked = np.ma.masked_where(~np.isfinite(values) | (values <= 0), values)
        else:
            maximum = max(float(np.max(np.abs(values))), np.finfo(float).eps)
            domain = (limits or {}).get(key, (-maximum, maximum))
            norm = Normalize(*domain)
            masked = np.ma.masked_invalid(values)
        artist = axis.pcolormesh(x, vertical, masked, cmap=cmap, norm=norm, shading="flat", rasterized=True)
        axis.set_aspect("equal", adjustable="box")
        if masked.count():
            figure.colorbar(artist, ax=axis, pad=0.02, shrink=0.9)
        else:
            axis.text(0.5, 0.5, "No magnetic field" if key == "beta" else "No positive values",
                      transform=axis.transAxes, ha="center", va="center")
        axis.set_title(title)
        axis.set_xlabel(r"Cylindrical radius $R$" if snapshot.metadata["geometry"] == "spherical" else r"$x_1$")
        axis.set_ylabel(r"Height $z$" if snapshot.metadata["geometry"] == "spherical" else r"$x_2$")
        if window and "xlim" in window:
            axis.set_xlim(*window["xlim"])
        if window and "ylim" in window:
            axis.set_ylim(*window["ylim"])
    # Unit-length arrows show direction, not field strength, on the density panel.
    bx, by, bz = _vector(snapshot, "B")
    horizontal_field = bx
    vertical_field = bz if snapshot.metadata["geometry"] == "spherical" else by
    magnitude = np.hypot(horizontal_field, vertical_field)
    if np.any(magnitude > 0):
        xx, zz = _coordinates(snapshot)
        # Select actual cell centres near a uniform physical sampling grid. A
        # simple index stride would pile up arrows at the inner log-r boundary.
        xmin, xmax = (window or {}).get("xlim", (float(x.min()), float(x.max())))
        zmin, zmax = (window or {}).get("ylim", (float(vertical.min()), float(vertical.max())))
        dx = max((xmax - xmin) / 26, 1e-12)
        dz = max((zmax - zmin) / 18, 1e-12)
        selected: dict[tuple[int, int], tuple[float, int]] = {}
        for index, (px, pz, strength) in enumerate(zip(xx.ravel(), zz.ravel(), magnitude.ravel())):
            if strength <= 0 or px < xmin or px > xmax or pz < zmin or pz > zmax:
                continue
            ix, iz = int((px - xmin) / dx), int((pz - zmin) / dz)
            distance = ((px - (xmin + (ix + 0.5) * dx)) / dx)**2 + ((pz - (zmin + (iz + 0.5) * dz)) / dz)**2
            previous = selected.get((ix, iz))
            if previous is None or distance < previous[0]:
                selected[(ix, iz)] = distance, index
        indices = np.asarray([entry[1] for entry in selected.values()], dtype=int)
        ux = np.divide(horizontal_field, magnitude, out=np.zeros_like(magnitude), where=magnitude > 0).ravel()[indices]
        uz = np.divide(vertical_field, magnitude, out=np.zeros_like(magnitude), where=magnitude > 0).ravel()[indices]
        span = max(xmax - xmin, zmax - zmin, 1e-12)
        if len(indices):
            axes[0, 0].quiver(xx.ravel()[indices], zz.ravel()[indices], ux, uz, color="white", alpha=0.7,
                              angles="xy", scale_units="xy", scale=35 / span, width=0.0022,
                              pivot="mid", headwidth=3.2)
            axes[0, 0].text(0.02, 0.02, "Arrows: poloidal B direction", transform=axes[0, 0].transAxes,
                            fontsize=8, color="white", bbox={"facecolor": "black", "alpha": 0.35, "edgecolor": "none"})
    figure.suptitle(f"PyHMD · {snapshot.metadata.get('case', 'MHD')} · t = {snapshot.time:.6g} · step {snapshot.step}", fontsize=14)
    figure.savefig(output, dpi=170)
    plt.close(figure)
    return output


def render_history(run_dir: str | Path, output: str | Path) -> Path | None:
    path = Path(run_dir) / "diagnostics.csv"
    if not path.exists():
        return None
    table, _ = _read_table(path)
    if not len(table):
        return None
    required = {"time", "mass", "energy", "min_rho", "min_pressure", "max_divB", "l1_divB", "dt", "floor_events"}
    if not required.issubset(table.dtype.names or ()):
        raise ValueError(f"diagnostics.csv missing columns {sorted(required - set(table.dtype.names or ())) }")
    plt, _, _ = _plot_imports()
    figure, axes = plt.subplots(2, 2, figsize=(11.8, 7.6), constrained_layout=True)
    time = table["time"]
    for field, label in (("mass", "Mass"), ("energy", "Total energy")):
        baseline = float(table[field][0])
        scale = abs(baseline) if baseline else 1.0
        axes[0, 0].plot(time, (table[field] - baseline) / scale, label=label)
    axes[0, 0].set_title("Change relative to initial integral")
    axes[0, 0].set_ylabel(r"$(Q-Q_0)/|Q_0|$")
    axes[0, 0].legend()
    for field, label in (("max_divB", r"max $|\nabla\cdot B|$"), ("l1_divB", r"mean $|\nabla\cdot B|$")):
        axes[0, 1].plot(time, np.maximum(table[field], 1e-30), label=label)
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_title("Magnetic divergence (code units)")
    axes[0, 1].legend()
    for field, label in (("min_rho", "Minimum density"), ("min_pressure", "Minimum pressure")):
        axes[1, 0].plot(time, np.maximum(table[field], 1e-30), label=label)
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_title("Positivity")
    axes[1, 0].legend()
    axes[1, 1].plot(time, table["dt"], color="#286090", label="dt")
    axes[1, 1].set_title("Time step and floor interventions")
    axes[1, 1].set_ylabel("dt")
    floors_axis = axes[1, 1].twinx()
    floors_axis.plot(time, table["floor_events"], color="#bf4b45", label="Floor events", alpha=0.8)
    floors_axis.set_ylabel("Floor events", color="#bf4b45")
    floors_axis.set_ylim(0, max(1.0, 1.1 * float(np.max(table["floor_events"]))))
    for axis in axes.ravel():
        axis.set_xlabel("Time (code units)")
        axis.grid(alpha=0.22)
    figure.suptitle("PyHMD diagnostics · open boundaries and gravity can change global integrals", fontsize=12)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=170)
    plt.close(figure)
    return output


def compare_runs(first: str | Path, second: str | Path, selection: str = "latest", rtol: float = 1e-9,
                 atol: float = 1e-12, time_atol: float = 1e-12) -> dict:
    left, right = load_snapshot(first, selection), load_snapshot(second, selection)
    if left.shape != right.shape or left.metadata["geometry"] != right.metadata["geometry"]:
        raise ValueError("Cannot compare different geometries or grid dimensions")
    if not math.isclose(left.time, right.time, rel_tol=0.0, abs_tol=time_atol):
        raise ValueError(f"Snapshot times differ: {left.time} versus {right.time}")
    for coordinate in ("x1", "x2"):
        if not np.allclose(left.data[coordinate], right.data[coordinate], rtol=1e-13, atol=1e-14):
            raise ValueError("Cannot compare different cell coordinates")
    for key in ("x1_faces", "x2_faces"):
        if not np.allclose(left.metadata[key], right.metadata[key], rtol=1e-13, atol=1e-14):
            raise ValueError("Cannot compare different cell faces")
    report = {"first": str(left.run_dir), "second": str(right.run_dir), "time": left.time,
              "rtol": rtol, "atol": atol, "passed": True, "fields": {}}
    settings = ("gamma", "eta", "eta_model", "cool_ohmic", "entropy_switch", "ct_energy_correction")
    report["settings_differences"] = {key: [left.metadata.get(key), right.metadata.get(key)]
                                      for key in settings if left.metadata.get(key) != right.metadata.get(key)}
    optional = tuple(field for field in ("entropy", "total_energy", "internal_energy")
                     if field in left.data and field in right.data)
    for field in FIELDS + optional:
        a, b = left.data[field], right.data[field]
        difference = np.abs(a - b)
        maximum = float(np.max(difference))
        denominator = max(float(np.max(np.abs(a))), np.finfo(float).tiny)
        passed = bool(np.all(difference <= atol + rtol * np.abs(a)))
        report["fields"][field] = {"max_abs": maximum,
                                    "normalized_linf": maximum / denominator,
                                    "rms_abs": float(np.sqrt(np.mean(difference**2))), "passed": passed}
        report["passed"] = report["passed"] and passed
    return report


def animate(run_dir: str | Path, output: str | Path, fps: int = 10, stride: int = 1,
            max_frames: int | None = None, window: dict | None = None) -> Path:
    run_dir = Path(run_dir).resolve()
    files = snapshot_files(run_dir)[::stride]
    if max_frames is not None:
        files = files[:max_frames]
    output = Path(output).resolve()
    if output.suffix.lower() not in (".gif", ".mp4"):
        raise ValueError("Animation output must end in .gif or .mp4")
    output.parent.mkdir(parents=True, exist_ok=True)
    # Fixed colour ranges prevent apparent fluctuations caused by autoscaling.
    limits: dict[str, tuple[float, float]] = {}
    velocity_max = 0.0
    for file in files:
        snap = load_snapshot(run_dir, file)
        b2 = sum(snap.data[key]**2 for key in ("Br", "Btheta", "Bphi"))
        arrays = {"rho": snap.data["rho"], "pressure": snap.data["pressure"],
                  "beta": np.divide(2 * snap.data["pressure"], b2, out=np.full(snap.shape, np.nan), where=b2 > 0)}
        for key, values in arrays.items():
            domain = _log_limits(values)
            if domain:
                previous = limits.get(key, domain)
                limits[key] = min(previous[0], domain[0]), max(previous[1], domain[1])
        velocity_max = max(velocity_max, float(np.max(np.abs(snap.data["vphi"]))))
    velocity_max = max(velocity_max, np.finfo(float).eps)
    limits["vphi"] = -velocity_max, velocity_max
    with tempfile.TemporaryDirectory(prefix="pymhd-frames-") as temporary:
        frame_paths = []
        for index, file in enumerate(files):
            frame = Path(temporary) / f"frame_{index:06d}.png"
            render_snapshot(load_snapshot(run_dir, file), frame, limits, window)
            frame_paths.append(frame)
        if output.suffix.lower() == ".gif":
            from PIL import Image
            frames = [Image.open(path).convert("RGB") for path in frame_paths]
            try:
                frames[0].save(output, save_all=True, append_images=frames[1:], duration=round(1000 / fps), loop=0)
            finally:
                for frame in frames:
                    frame.close()
        else:
            try:
                subprocess.run(["ffmpeg", "-hide_banner", "-loglevel", "error", "-y", "-framerate", str(fps),
                                "-i", str(Path(temporary) / "frame_%06d.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p",
                                "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", str(output)], check=True)
            except FileNotFoundError as error:
                raise ValueError("ffmpeg is required for MP4; use .gif or install ffmpeg") from error
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    render = subparsers.add_parser("render", help="Create a final-state PNG, diagnostics PNG, and VTK file")
    render.add_argument("run_dir", type=Path)
    render.add_argument("--snapshot", default="latest")
    render.add_argument("--output", type=Path, help="Output directory (default RUN_DIR/visualization)")
    render.add_argument("--xlim", nargs=2, type=float, metavar=("MIN", "MAX"), help="Physical horizontal limits (R for spherical)")
    render.add_argument("--ylim", nargs=2, type=float, metavar=("MIN", "MAX"), help="Physical vertical limits (z for spherical)")
    vtk = subparsers.add_parser("vtk", help="Export one snapshot for ParaView")
    vtk.add_argument("run_dir", type=Path)
    vtk.add_argument("--snapshot", default="latest")
    vtk.add_argument("--output", type=Path, required=True)
    animation = subparsers.add_parser("animate", help="Make a GIF or ffmpeg MP4 with fixed colour ranges")
    animation.add_argument("run_dir", type=Path)
    animation.add_argument("--output", type=Path, required=True)
    animation.add_argument("--fps", type=int, default=10)
    animation.add_argument("--stride", type=int, default=1)
    animation.add_argument("--max-frames", type=int)
    animation.add_argument("--xlim", nargs=2, type=float, metavar=("MIN", "MAX"))
    animation.add_argument("--ylim", nargs=2, type=float, metavar=("MIN", "MAX"))
    compare = subparsers.add_parser("compare", help="Compare CPU/GPU or MPI decompositions at the same time")
    compare.add_argument("first", type=Path)
    compare.add_argument("second", type=Path)
    compare.add_argument("--snapshot", default="latest")
    compare.add_argument("--rtol", type=float, default=1e-9)
    compare.add_argument("--atol", type=float, default=1e-12)
    compare.add_argument("--time-atol", type=float, default=1e-12)
    compare.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        window = {}
        for name in ("xlim", "ylim"):
            domain = getattr(args, name, None)
            if domain is not None:
                if not np.isfinite(domain).all() or domain[1] <= domain[0]:
                    raise ValueError(f"{name} requires finite increasing limits")
                window[name] = domain
        if args.command == "render":
            output_dir = args.output or args.run_dir / "visualization"
            snapshot = load_snapshot(args.run_dir, args.snapshot)
            paths = [render_snapshot(snapshot, output_dir / "overview.png", window=window),
                     render_history(args.run_dir, output_dir / "history.png"),
                     write_vtk(snapshot, output_dir / "final.vtk")]
            print(json.dumps({"outputs": [str(path.resolve()) for path in paths if path]}, indent=2))
        elif args.command == "vtk":
            print(write_vtk(load_snapshot(args.run_dir, args.snapshot), args.output))
        elif args.command == "animate":
            if args.fps < 1 or args.stride < 1 or (args.max_frames is not None and args.max_frames < 1):
                raise ValueError("fps, stride, and max-frames must be positive")
            print(animate(args.run_dir, args.output, args.fps, args.stride, args.max_frames, window))
        else:
            if min(args.rtol, args.atol, args.time_atol) < 0:
                raise ValueError("Comparison tolerances must be nonnegative")
            report = compare_runs(args.first, args.second, args.snapshot, args.rtol, args.atol, args.time_atol)
            content = json.dumps(report, indent=2, allow_nan=False)
            if args.output:
                args.output.parent.mkdir(parents=True, exist_ok=True)
                args.output.write_text(content + "\n")
            print(content)
            return 0 if report["passed"] else 1
        return 0
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
        print(f"pyhmd render: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
