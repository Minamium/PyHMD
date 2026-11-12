---
name: pymhd-use
description: Prepare, run, resume, or inspect existing PyMHD simulations. Use for job status, time-to-target estimates, diagnostics, validation, and plots. Use pymhd-develop for changes to the solver or API.
---

# PyMHD simulation use

## Status and runtime

- For an existing HPC job, use its job ID and run directory. Read the scheduler
  state, elapsed wall time, and exit code. Query finished-job history if the
  active queue is empty (PBS example: `qstat -x -f JOB_ID`).
- For a suite, read its progress file and the first failed stage log. For an
  integration, read `result.json` or native `status.json`, then the last accepted
  time, step, and solver wall time in `diagnostics.csv` or `run.log`. Keep queue
  wait, batch wall time, and solver wall time distinct.
- For time to a target, prefer the exact diagnostic row at that simulation time
  from a comparable run. Otherwise use the recent accepted Δt/Δwall for an
  active run or bracket the target with historical samples and label any
  interpolation or extrapolation. Report the grid, hardware, rank count, and
  relevant settings. Adaptive dt, retries, and output cadence can change the rate.
  For a stopped job, give its stop stage and separate any historical estimate.
  [Recorded experiments](../../../docs/experiment-results.md)

## Prepare, run, or resume

- Resolve `Config`/`RunSpec`, runtime, and output settings from the requested
  preset. Partial dictionaries merge into presets; group dataclasses replace
  the group. Standard SF2014 and strong diffusivity are separate experiments.
  [Configuration](../../../docs/configuration.md)
- Build the requested backend explicitly, then run `pymhd plan` to verify the
  resolved command and capabilities. Use a fresh output directory for each run.
  Select output and checkpoint cadence and enable transport ledgers before
  launch if needed. [Working examples](../../../README.md)
- In an existing allocation, preserve `CUDA_VISIBLE_DEVICES` and map one native
  MPI rank per GPU. The Python engine launches MPI itself. The CUDA/MPI build
  needs Open MPI-compatible `--showme` flags and the appropriate SM target.
  [HPC](../../../docs/hpc.md)
- Resume with `RuntimeConfig.restart` or INI `restart` and a fresh output
  directory. Preserve checkpoint signatures, rank decomposition, physics and
  scheduling. [Restart contract](../../../docs/results.md)

## Interpret saved results

Load with `Result.load()` or use `pymhd inspect`, `validate` and `render`. Check
terminal state and requested end time before reporting completion. For shrinking
dt, inspect the active CFL constraint, characteristic speeds, resistivity,
retries, floors and correction ledgers. Use accepted face fluxes for transport
when available. Report numerical integrity and paper agreement separately.
[Results](../../../docs/results.md) · [Scientific scope](../../../docs/scientific-status.md)
