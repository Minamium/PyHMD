---
name: pyhmd-use
description: Configure, build, run, resume, inspect, validate, and visualize PyHMD simulations using existing capabilities. Use for simulation preparation and saved-result analysis; use pyhmd-develop for changing implementation.
---

# PyHMD simulation use

Use the requested physical conditions and the current checkout or installed
package. Keep the user's response language. The [README](../../../README.md)
contains executable starting examples; load further references only as needed.

## Choose the operation

| Request | Read |
| --- | --- |
| Settings, grid, physics, IC/BC, recording | [Configuration](../../../docs/configuration.md), relevant preset in `src/pyhmd/presets/` |
| CPU/CUDA/MPI build or allocated HPC run | [HPC](../../../docs/hpc.md), selected example in `examples/` |
| Existing output, restart, stop/failure | [Results](../../../docs/results.md) |
| Paper comparison or physical interpretation | [Scientific status](../../../docs/scientific-status.md), [Validation](../../../docs/validation.md) |

## Prepare and run

- Normalize the request through `Config`/`RunSpec`, `RuntimeConfig`, and
  `ResultOptions`. Partial dictionaries merge into problem presets; a group
  dataclass replaces the full group. Resolve paths from the execution cwd.
- Check the actual supported capabilities before proposing settings. Current
  numerics are FP64 PPM/HLL/SSP-RK3/CT in 2D Cartesian or spherical axisymmetry.
  Standard SF2014 and strong resistivity are distinct experiments.
- Select snapshots, checkpoints, and optional transport diagnostics before
  launch. Missing ledgers cannot later be presented as measured transport.
- Build explicitly, then use `plan` for native validation and command inspection.
  `plan` does not launch or create a result directory. An absent binary needs a
  build; a backend mismatch is not an instruction to fall back to CPU.
- Use an absent or empty output directory. `run` blocks until the native process
  ends and records stdout/stderr in `run.log`; record the actual result path.

For an authorized small CPU run from a checkout:

```bash
pyhmd build --backend cpu
pyhmd plan examples/thin_disk_cpu.json --out runs/my_disk
pyhmd run examples/thin_disk_cpu.json --out runs/my_disk
pyhmd inspect runs/my_disk
pyhmd validate runs/my_disk
pyhmd render runs/my_disk
```

Choose a fresh name if that result already exists. Do not run a new simulation
merely to answer a saved-result or command-preparation request.

## HPC constraints

Use the user's existing execution scope and site allocation. Preserve
`CUDA_VISIBLE_DEVICES`; keep one native MPI rank per allocated GPU. Do not wrap
the Python Engine in MPI: it appends `-np N` to its launcher and starts native
ranks itself. Use a compatible MPI launcher, or direct native `srun` on a site
that supports it. Each rank needs at least three first-coordinate cells.

The CUDA+MPI build expects Open MPI-compatible `--showme` flags. Verify the actual
toolchain/SM target. Keep site credentials/accounts/hostnames out of public
artifacts. Preparing a run or loading this skill does not authorize acquiring
new resources, stopping unrelated jobs, or publishing changes.

## Inspect, resume, and explain

- Start with `result.json`, native `status.json`, `run.log`, and `diagnostics.csv`.
  Treat `incomplete`/exit 3 and `stopped`/4 as partial outcomes, not completed runs.
  `SimulationError.result` retains failed-run access.
- Load existing output through `Result.load()`, validate numerical integrity, and
  render offline. Legacy output needs `expect_time` for validation. Validation
  rejects floor events and does not certify a paper reproduction.
- Resume through `RuntimeConfig.restart` or INI `restart`, with a new output
  directory and compatible settings, rank decomposition, output interval, and
  source signatures. The Python CLI has no `--restart` switch. Do not modify
  checkpoint signatures or assume pre-migration checkpoints can resume.
- For a shrinking dt, inspect the active CFL restriction, evolving speeds,
  resistivity/grid scales, retries, floors, and correction ledgers before making
  a physical claim. Do not change floors, closure, boundaries, or numerical
  settings just to obtain a completed/validated run.
- Compare identical numerical/output settings, and report accepted steps,
  physical time advanced, wall time, and numerical quality together. Small divB
  or rapid execution alone is not evidence of the correct solution.

Return the requested configuration, commands, result analysis, or figures. State
what was prepared versus actually executed and which hardware was checked.
Quantitative SF2014 reproduction remains unverified; check the recorded fidelity
gaps before recommending a longer or higher-resolution run as the solution.
