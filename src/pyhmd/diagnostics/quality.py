#!/usr/bin/env python3
"""Validate a completed PyHMD or legacy PyMHD run and compare a CPU reference.

Uses only the Python standard library, so no remote NumPy/venv is required.
"""
import argparse
from contextlib import contextmanager
import csv
from itertools import zip_longest
import json
import math
from pathlib import Path
import re
import sys

FIELDS = ("x1", "x2", "rho", "vr", "vtheta", "vphi", "pressure", "Br",
          "Btheta", "Bphi", "psi", "divB", "entropy", "total_energy", "internal_energy")
REQUIRED_META = ("case", "geometry", "nx1", "nx2", "gamma", "eta_model", "cool_ohmic",
                 "entropy_switch", "ct_energy_correction")
SETTINGS = REQUIRED_META + ("eta", "cfl", "entropy_ratio", "entropy_shock_threshold",
    "entropy_contact_threshold", "max_step_retries", "ohmic_cooling_discretization",
    "reconstruction", "riemann", "time_integrator", "divergence_control",
    "initial_equilibrium", "gm", "epsilon", "mu0", "mu_ref", "b0", "alpha_m",
    "chi", "corona_delta", "wave_amplitude", "x1_min", "x1_max", "x2_min", "x2_max")
SETTINGS += ("metric_ppm", "retry_flat_same_dt", "negative_pressure_policy", "cfl_definition", "diffusivity_normalization", "source_fingerprint")
SETTINGS += ("internal_energy_sync", "dual_energy_eta")
SETTINGS += ("ct_emf_scheme",)
SETTINGS += ("cfl_staggered",)
SETTINGS += ("inner_normal_field",)


class InvalidRun(Exception):
    pass


def require(condition, message):
    if not condition:
        raise InvalidRun(message)


def same_time(value, expected):
    return math.isfinite(value) and math.isclose(value, expected, rel_tol=1e-11, abs_tol=1e-12)


def number(value, where):
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError):
        raise InvalidRun("{}: invalid number {!r}".format(where, value))
    require(math.isfinite(result), "{}: nonfinite value {!r}".format(where, value))
    return result


def read_json(path):
    with path.open() as handle:
        result = json.load(handle)
    require(isinstance(result, dict), "{}: expected a JSON object".format(path))
    return result


@contextmanager
def snapshot_reader(path):
    with path.open(newline="") as handle:
        comments = []
        while True:
            start = handle.tell()
            line = handle.readline()
            require(bool(line), "{}: empty snapshot".format(path))
            if line.startswith("#") or not line.strip():
                comments.append(line)
            else:
                handle.seek(start)
                break
        reader = csv.DictReader(handle)
        headers = reader.fieldnames or []
        require(len(headers) == len(set(headers)), "{}: duplicate CSV columns".format(path))
        require(set(FIELDS).issubset(headers), "{}: missing NV10 columns {}".format(
            path, sorted(set(FIELDS) - set(headers))))
        yield reader, " ".join(comments)


def finite_row(row, where):
    require(None not in row and all(v is not None for v in row.values()),
            "{}: malformed CSV row".format(where))
    return {key: number(value, "{} {}".format(where, key)) for key, value in row.items()}


def latest_snapshot(root):
    choices = []
    for path in root.glob("snapshot_*.csv"):
        match = re.fullmatch(r"snapshot_(\d+)\.csv", path.name)
        if match:
            choices.append((int(match.group(1)), path))
    require(bool(choices), "{}: no completed snapshot_*.csv files".format(root))
    return max(choices, key=lambda item: item[0])[1]


def initial_energy_integrals(root, meta, nx, ny):
    """Physical energy scales for judging nonconservative correction ledgers."""
    xf, yf = meta.get("x1_faces"), meta.get("x2_faces")
    require(isinstance(xf, list) and len(xf) == nx+1 and
            isinstance(yf, list) and len(yf) == ny+1, "missing initial-energy grid faces")
    gamma = number(meta["gamma"], "metadata gamma")
    require(gamma > 1, "gamma must exceed one")
    if meta["geometry"] == "spherical":
        dx = [(b*b*b-a*a*a)/3 for a,b in zip(xf,xf[1:])]
        dy = [math.cos(a)-math.cos(b) for a,b in zip(yf,yf[1:])]
    else:
        dx = [b-a for a,b in zip(xf,xf[1:])]
        dy = [b-a for a,b in zip(yf,yf[1:])]
    require(all(v > 0 and math.isfinite(v) for v in dx+dy), "invalid cell volumes")
    result = dict(thermal=0., magnetic=0., kinetic=0.)
    with snapshot_reader(root / "snapshot_000000.csv") as (reader, _):
        count = 0
        for k,row in enumerate(reader):
            require(k < nx*ny, "initial snapshot has too many cells")
            q = {key:number(row[key], "initial energy "+key) for key in
                 ("rho","pressure","Br","Btheta","Bphi","vr","vtheta","vphi")}
            require(q["rho"] > 0 and q["pressure"] > 0, "invalid initial gas state")
            volume = dx[k % nx]*dy[k // nx]
            result["thermal"] += q["pressure"]/(gamma-1)*volume
            result["magnetic"] += .5*sum(q[c]**2 for c in ("Br","Btheta","Bphi"))*volume
            result["kinetic"] += .5*q["rho"]*sum(q[c]**2 for c in ("vr","vtheta","vphi"))*volume
            count += 1
    require(count == nx*ny, "initial snapshot cell count differs")
    return result


def validate(root, expected, progress=None):
    meta = read_json(root / "metadata.json")
    for key in REQUIRED_META:
        require(key in meta, "{}: metadata missing {}".format(root, key))
    nx = number(meta["nx1"], "metadata nx1")
    ny = number(meta["nx2"], "metadata nx2")
    require(nx > 0 and ny > 0 and nx.is_integer() and ny.is_integer(), "invalid grid dimensions")
    status = read_json(root / "status.json")
    require(status.get("complete") is True, "{}: status.complete is not true".format(root))
    require(status.get("reason") == "t_end", "{}: run did not finish by t_end".format(root))
    time = number(status.get("time"), "status time")
    steps = number(status.get("step"), "status step")
    require(same_time(time, expected), "{}: completed time {} != {}".format(root, time, expected))
    require(steps > 0 and steps.is_integer(), "{}: no valid completed steps".format(root))
    if progress is not None:
        progress.update(time=time, steps=int(steps))
    floors = max_div = unreliable_shocks = 0.0
    rows = 0
    previous_time = -math.inf
    with (root / "diagnostics.csv").open(newline="") as handle:
        reader = csv.DictReader(handle)
        headers = reader.fieldnames or []
        require(len(headers) == len(set(headers)), "{}: duplicate diagnostic columns".format(root))
        for rows, row in enumerate(reader, 1):
            q = finite_row(row, "{} diagnostics row {}".format(root, rows))
            if rows == 1:
                initial_diagnostics = q
            for key in ("time", "step", "min_rho", "min_pressure", "floor_events", "max_divB"):
                require(key in q, "diagnostics missing {}".format(key))
            require(q["min_rho"] > 0 and q["min_pressure"] > 0, "diagnostics has nonpositive rho/p")
            require(q["floor_events"] >= 0 and q["max_divB"] >= 0, "invalid negative diagnostic count")
            require(q["time"] >= previous_time, "diagnostic times are not monotone")
            previous_time = q["time"]
            floors = max(floors, q["floor_events"])
            max_div = max(max_div, q["max_divB"])
            unreliable_shocks = max(unreliable_shocks,q.get("dual_unreliable_shock_cells",0.0))
            final = q
    require(rows > 0, "{}: empty diagnostics".format(root))
    if progress is not None:
        progress.update(floors=floors, max_divB=max_div)
    require(floors == 0, "{}: positivity floors occurred ({:g})".format(root, floors))
    require(max_div <= 1e-10, "{}: maximum divB exceeds 1e-10 ({:.3g})".format(root, max_div))
    require(same_time(final["time"], expected) and final["step"] == steps,
            "{}: final diagnostics do not match completion status".format(root))
    snapshot = latest_snapshot(root)
    count = 0
    with snapshot_reader(snapshot) as (reader, comment):
        tm = re.search(r"\btime=([^\s]+)", comment)
        st = re.search(r"\bstep=(\d+)\b", comment)
        require(tm is not None and same_time(number(tm.group(1), "snapshot time"), expected),
                "{}: snapshot time mismatch".format(snapshot))
        require(st is not None and int(st.group(1)) == steps, "{}: snapshot step mismatch".format(snapshot))
        for count, row in enumerate(reader, 1):
            q = finite_row(row, "{} row {}".format(snapshot.name, count))
            require(q["rho"] > 0 and q["pressure"] > 0, "{}: nonpositive snapshot rho/p".format(snapshot))
            max_div = max(max_div, abs(q["divB"]))
    require(count == int(nx * ny), "{}: snapshot has {} cells, expected {}".format(root, count, int(nx * ny)))
    if progress is not None:
        progress["max_divB"] = max_div
    require(max_div <= 1e-10, "{}: snapshot divB exceeds 1e-10".format(root))
    initial_energy = initial_energy_integrals(root, meta, int(nx), int(ny))
    correction_keys = ("ct_energy", "ct_boundary_energy", "internal_energy", "entropy_energy", "floor_energy")
    corrections = {k:final.get(k,0.)-initial_diagnostics.get(k,0.) for k in correction_keys}
    # This is a screening scale, not a physical-accuracy acceptance criterion.
    # Absolute global terms still cannot reveal cancellation within each term.
    correction_scale = sum(abs(v) for v in corrections.values())/initial_energy["thermal"]
    warnings = []
    if unreliable_shocks:
        warnings.append("Modern dual energy was used in sampled shock-flagged low-thermal cells; shock heating needs separate physical validation.")
    if correction_scale > .1:
        warnings.append("Global nonconservative correction terms exceed 10% of initial thermal energy; assess convergence and local correction maps.")
    return meta, snapshot, {"directory": str(root), "time": time, "steps": int(steps),
        "floors": floors, "max_divB": max_div, "snapshot": snapshot.name, "cells": count,
        "validation_scope": "numerical_integrity_only", "scientific_reproduction_verified": False,
        "max_sampled_dual_unreliable_shock_cells": unreliable_shocks,
        "initial_energy_integrals": initial_energy, "interval_energy_corrections": corrections,
        "sum_abs_global_corrections_over_initial_thermal": correction_scale,
        "correction_screening_threshold": .1, "warnings": warnings}


def settings_equal(a, b):
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=1e-13, abs_tol=1e-14)
    return a == b


def compare(meta, snapshot, ref_meta, reference):
    differences = {key: [meta.get(key), ref_meta.get(key)] for key in SETTINGS
                   if not settings_equal(meta.get(key), ref_meta.get(key))}
    for key in ("x1_faces", "x2_faces"):
        a, b = meta.get(key), ref_meta.get(key)
        if not isinstance(a, list) or not isinstance(b, list) or len(a) != len(b) or not all(
                settings_equal(x, y) for x, y in zip(a, b)):
            differences[key] = "different face coordinates"
    require(not differences, "CPU/GPU settings differ: {}".format(json.dumps(differences, ensure_ascii=False)))
    tolerances = {}
    for field in FIELDS:
        rtol, atol = 1e-7, 1e-12
        if field in ("pressure", "internal_energy"):
            atol = 1e-20
        elif field == "rho":
            atol = 1e-18
        elif field in ("x1", "x2"):
            rtol, atol = 1e-13, 1e-14
        tolerances[field] = {"max_abs": 0.0, "max_scaled_error": 0.0,
                             "failures": 0, "rtol": rtol, "atol": atol}
    with snapshot_reader(snapshot) as (left, _), snapshot_reader(reference) as (right, _):
        for index, pair in enumerate(zip_longest(left, right), 1):
            require(None not in pair, "comparison snapshot lengths differ")
            a, b = (finite_row(row, "comparison row {}".format(index)) for row in pair)
            for field, metric in tolerances.items():
                error = abs(a[field] - b[field])
                limit = metric["atol"] + metric["rtol"] * abs(b[field])
                metric["max_abs"] = max(metric["max_abs"], error)
                metric["max_scaled_error"] = max(metric["max_scaled_error"], error / limit)
                metric["failures"] += int(error > limit)
    return tolerances


def json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--expect-time", type=float, required=True)
    parser.add_argument("--compare", type=Path, metavar="CPU_DIR")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    root = args.run_directory.resolve()
    destination = args.report or root / "validation.json"
    report = {"passed": False, "directory": str(root), "expected_time": args.expect_time,
              "time": None, "steps": None, "floors": None, "max_divB": None,
              "fielddiffs": {}, "errors": []}
    try:
        require(math.isfinite(args.expect_time) and args.expect_time > 0, "expected time must be finite and positive")
        meta, snapshot, summary = validate(root, args.expect_time, report)
        report.update(summary)
        if args.compare:
            ref_meta, reference, ref_summary = validate(args.compare.resolve(), args.expect_time)
            report["reference"] = ref_summary
            report["fielddiffs"] = compare(meta, snapshot, ref_meta, reference)
            bad = [name for name, result in report["fielddiffs"].items() if result["failures"]]
            require(not bad, "CPU/GPU field comparison failed: {}".format(", ".join(bad)))
        report["passed"] = True
    except (InvalidRun, OSError, ValueError, KeyError, TypeError, OverflowError, csv.Error) as error:
        report["errors"].append(str(error))
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(json_safe(report), indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    except OSError as error:
        print("FAIL: cannot write validation report: {}".format(error), file=sys.stderr)
        return 2
    if report["passed"]:
        print("PASS {} t={:g} steps={} floors=0; report={}".format(root.name, report["time"], report["steps"], destination))
        return 0
    print("FAIL {}: {}; report={}".format(root.name, "; ".join(report["errors"]), destination), file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
