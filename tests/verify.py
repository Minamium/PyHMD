#!/usr/bin/env python3
"""Run physical regression checks against analytic solutions and invariants.

The Brio-Wu and thin-disk checks are explicitly smoke/invariant checks, not
reference-solution verification.  No tolerance is inferred from the solver code.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from pymhd_viz import FIELDS, _read_table, cell_volumes, load_snapshot  # noqa: E402


def conserved(snapshot) -> dict[str, float]:
    volume = cell_volumes(snapshot)
    rho, pressure = snapshot.data["rho"], snapshot.data["pressure"]
    velocity2 = sum(snapshot.data[field]**2 for field in ("vr", "vtheta", "vphi"))
    magnetic2 = sum(snapshot.data[field]**2 for field in ("Br", "Btheta", "Bphi"))
    energy = pressure / (float(snapshot.metadata["gamma"]) - 1) + 0.5 * rho * velocity2 + 0.5 * magnetic2
    if "total_energy" in snapshot.data:
        energy = snapshot.data["total_energy"]
    return {"mass": float(np.sum(rho * volume)), "energy": float(np.sum(energy * volume))}


def relative_change(first: float, second: float) -> float:
    return abs(second - first) / max(abs(first), np.finfo(float).tiny)


def common_checks(run_dir: Path) -> tuple[dict, object, object]:
    initial, final = load_snapshot(run_dir, "first"), load_snapshot(run_dir, "latest")
    if final.time <= initial.time:
        raise AssertionError("No time advancement: initial and final snapshots have identical time")
    if np.min(final.data["rho"]) <= 0 or np.min(final.data["pressure"]) <= 0:
        raise AssertionError("Non-positive density or gas pressure")
    diagnostics, _ = _read_table(run_dir / "diagnostics.csv")
    maximum_floors = float(np.max(diagnostics["floor_events"]))
    if maximum_floors != 0:
        raise AssertionError(f"Positivity floors were triggered ({maximum_floors:g} events)")
    report = {"time": final.time, "step": final.step,
              "minimum_density": float(np.min(final.data["rho"])),
              "minimum_pressure": float(np.min(final.data["pressure"])),
              "maximum_divB": float(np.max(np.abs(final.data["divB"]))),
              "floor_events": maximum_floors}
    return report, initial, final


def uniform_check(run_dir: Path) -> dict:
    report, initial, final = common_checks(run_dir)
    expected = {"rho": 1.0, "pressure": 1.0, "vr": 0.0, "vtheta": 0.0, "vphi": 0.0,
                "Br": 0.0, "Btheta": 0.0, "Bphi": 0.0, "psi": 0.0}
    errors = {field: float(np.max(np.abs(final.data[field] - value))) for field, value in expected.items()}
    maximum = max(errors.values())
    if maximum > 1e-11:
        raise AssertionError(f"Static uniform spherical state changed by {maximum:.6g}; tolerance 1e-11")
    report.update({"kind": "analytic stationary spherical free-stream", "maximum_absolute_error": maximum,
                   "field_errors": errors})
    return report


def alfven_check(run_dir: Path) -> dict:
    report, initial, final = common_checks(run_dir)
    amplitude = float(final.metadata.get("wave_amplitude", 0.1))
    x = final.data["x1"]
    phase = 2 * np.pi * (x - final.time)
    expected_by, expected_bz = amplitude * np.sin(phase), amplitude * np.cos(phase)
    volume = cell_volumes(final)
    error2 = ((final.data["Btheta"] - expected_by)**2 + (final.data["Bphi"] - expected_bz)**2
              + (final.data["vtheta"] + expected_by)**2 + (final.data["vphi"] + expected_bz)**2)
    normalized_rms = float(np.sqrt(np.sum(error2 * volume) / (4 * np.sum(volume))) / amplitude)
    # 2.5 percent is deliberately looser than a high-resolution convergence run;
    # the second resolution in --suite full checks that error actually decreases.
    if normalized_rms > 0.025:
        raise AssertionError(f"Alfven wave normalized RMS error {normalized_rms:.6g} exceeds 0.025")
    start, end = conserved(initial), conserved(final)
    for integral in ("mass", "energy"):
        change = relative_change(start[integral], end[integral])
        if change > 1e-9:
            raise AssertionError(f"Periodic Alfven {integral} changed by {change:.6g}")
        report[f"relative_{integral}_change"] = change
    report.update({"kind": "analytic circularly polarized Alfven wave", "normalized_rms_error": normalized_rms,
                   "reference": "Bx=1,rho=1,vx=0; By=A sin(2pi(x-t)), Bz=A cos(2pi(x-t)), v_perp=-B_perp"})
    return report


def diffusion_check(run_dir: Path) -> dict:
    report, initial, final = common_checks(run_dir)
    amplitude = float(final.metadata.get("wave_amplitude", 1e-3))
    eta = float(final.metadata.get("eta", 0.01))
    x = final.data["x1"]
    sine = np.sin(2 * np.pi * x)
    volume = cell_volumes(final)
    projected_amplitude = float(np.sum(final.data["Bphi"] * sine * volume) / np.sum(sine**2 * volume))
    expected_amplitude = amplitude * np.exp(-eta * (2 * np.pi)**2 * final.time)
    relative_error = abs(projected_amplitude - expected_amplitude) / abs(expected_amplitude)
    if relative_error > 2e-3:
        raise AssertionError(f"Resistive decay amplitude error {relative_error:.6g} exceeds 0.002")
    if projected_amplitude >= amplitude:
        raise AssertionError("Resistive mode did not decay")
    start, end = conserved(initial), conserved(final)
    energy_error = relative_change(start["energy"], end["energy"])
    if energy_error > 1e-9:
        raise AssertionError(f"Periodic resistive total energy changed by {energy_error:.6g}")
    report.update({"kind": "weak-field resistive Fourier decay", "eta": eta,
                   "measured_amplitude": projected_amplitude, "analytic_amplitude": float(expected_amplitude),
                   "relative_amplitude_error": float(relative_error), "relative_energy_change": energy_error,
                   "limitation": "Finite-amplitude compressible response is present; pure diffusion is an asymptotic weak-field reference."})
    return report


def brio_wu_check(run_dir: Path) -> dict:
    report, initial, final = common_checks(run_dir)
    start, end = conserved(initial), conserved(final)
    for integral in ("mass", "energy"):
        change = relative_change(start[integral], end[integral])
        if change > 1e-8:
            raise AssertionError(f"Brio-Wu {integral} changed by {change:.6g} before waves should reach the boundary")
        report[f"relative_{integral}_change"] = change
    bx_error = float(np.max(np.abs(final.data["Br"] - 0.75)))
    if bx_error > 1e-11:
        raise AssertionError(f"1D Brio-Wu normal magnetic field changed by {bx_error:.6g}")
    # A nontrivial solution must develop transverse velocity from the initial jump.
    velocity_peak = float(np.max(np.abs(final.data["vtheta"])))
    if velocity_peak < 1e-3:
        raise AssertionError("Brio-Wu transverse response is missing")
    report.update({"kind": "shock-tube positivity and conservation", "normal_field_error": bx_error,
                   "transverse_velocity_peak": velocity_peak,
                   "limitation": "No exact or independently generated MHD Riemann profile is used; this is not a shock-profile accuracy certificate."})
    return report


def thin_disk_check(run_dir: Path) -> dict:
    report, initial, final = common_checks(run_dir)
    volume = cell_volumes(final)
    speed = np.sqrt(sum(final.data[field]**2 for field in ("vr", "vtheta", "vphi")))
    if not np.isfinite(np.sum(volume * speed)):
        raise AssertionError("Thin-disk diagnostic is non-finite")
    report.update({"kind": "spherical disk execution smoke test", "maximum_speed": float(np.max(speed)),
                   "limitation": "A short run checks execution and positivity only. Disk equilibrium, long-term accretion and winds require refinement and boundary studies."})
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build" / "pymhd_cpu")
    parser.add_argument("--output", type=Path, default=ROOT / "runs" / "verification")
    parser.add_argument("--suite", choices=("quick", "full"), default="quick")
    parser.add_argument("--cases", default="uniform,alfven,diffusion,brio_wu,thin_disk", help="Comma-separated case names")
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument("--no-run", action="store_true", help="Assess existing output directories")
    args = parser.parse_args(argv)
    args.binary = args.binary.resolve()
    args.output = args.output.resolve()
    if not args.no_run and not args.binary.is_file():
        parser.error(f"Solver binary does not exist: {args.binary}")
    quick = args.suite == "quick"
    specifications = {
        "uniform": [(32 if quick else 64, 16 if quick else 32, 0.02 if quick else 0.2)],
        "alfven": [(64, 4, 0.1 if quick else 0.2)] + ([] if quick else [(128, 4, 0.2)]),
        "diffusion": [(64 if quick else 128, 4, 0.2)],
        "brio_wu": [(128 if quick else 256, 4, 0.05 if quick else 0.1)],
        "thin_disk": [(32 if quick else 64, 16 if quick else 32, 0.005 if quick else 0.05)],
    }
    checkers = {"uniform": uniform_check, "alfven": alfven_check, "diffusion": diffusion_check,
                "brio_wu": brio_wu_check, "thin_disk": thin_disk_check}
    cases = [case.strip() for case in args.cases.split(",") if case.strip()]
    if not cases:
        parser.error("At least one verification case is required")
    unknown = set(cases) - set(specifications)
    if unknown:
        parser.error(f"Unknown cases: {sorted(unknown)}")
    args.output.mkdir(parents=True, exist_ok=True)
    report = {"suite": args.suite, "binary": str(args.binary), "passed": True, "checks": []}
    alfven_errors = []
    for case in cases:
        for nx, ny, end_time in specifications[case]:
            destination = args.output / f"{case}_{nx}x{ny}"
            record = {"case": case, "nx1": nx, "nx2": ny, "output": str(destination)}
            try:
                if not args.no_run:
                    destination.mkdir(parents=True, exist_ok=True)
                    # Refuse stale snapshots because 'latest' could otherwise select
                    # a previous longer run and silently invalidate the checks.
                    if list(destination.glob("snapshot_*.csv")):
                        raise ValueError(f"Output already contains snapshots: {destination}; choose a fresh --output or --no-run")
                    command = [str(args.binary), "--config", str(ROOT / "configs" / f"{case}.ini"),
                               "--output", str(destination), "--nx1", str(nx), "--nx2", str(ny),
                               "--t-end", str(end_time)]
                    record["command"] = command
                    started = time.monotonic()
                    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                            text=True, timeout=args.timeout, cwd=ROOT)
                    (destination / "execution.log").write_text(result.stdout)
                    record["wall_seconds"] = time.monotonic() - started
                    if result.returncode:
                        raise AssertionError(f"Solver exited {result.returncode}; see {destination / 'execution.log'}")
                result = checkers[case](destination)
                if not np.isclose(result["time"], end_time, rtol=0, atol=1e-11):
                    raise AssertionError(f"Solver ended at {result['time']}, requested {end_time}")
                record.update(result)
                record["passed"] = True
                if case == "alfven":
                    alfven_errors.append((nx, result["normalized_rms_error"]))
            except (AssertionError, ValueError, OSError, KeyError, subprocess.TimeoutExpired) as error:
                record.update({"passed": False, "error": str(error)})
                report["passed"] = False
            report["checks"].append(record)
            print(f"{'PASS' if record['passed'] else 'FAIL'} {case} {nx}x{ny}" + (f": {record['error']}" if not record["passed"] else ""), flush=True)
    if len(alfven_errors) > 1:
        alfven_errors.sort()
        low, high = alfven_errors[0], alfven_errors[-1]
        ratio = high[1] / max(low[1], np.finfo(float).tiny)
        convergence = {"case": "alfven_refinement", "coarse_nx": low[0], "fine_nx": high[0],
                       "error_ratio": ratio, "observed_order": float(np.log(max(low[1], 1e-300) / max(high[1], 1e-300)) / np.log(high[0] / low[0])),
                       "passed": ratio < 0.9}
        report["checks"].append(convergence)
        report["passed"] = report["passed"] and convergence["passed"]
        print(f"{'PASS' if convergence['passed'] else 'FAIL'} Alfven refinement: ratio={ratio:.6g}")
    report_path = args.output / "verification.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"Report: {report_path}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
