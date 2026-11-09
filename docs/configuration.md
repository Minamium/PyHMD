# Configuration

[English](configuration.md) | [日本語](ja/configuration.md) | [Index](index.md)

`Config` describes the problem and numerics, `RuntimeConfig` describes execution,
and `ResultOptions` describes recording. All are frozen dataclasses. Unknown or
duplicate JSON keys, nonfinite numbers, and booleans passed as integers are rejected.

```json
{
  "schema_version": 1,
  "config": {
    "problem": "thin_disk",
    "grid": {"nx1": 64, "nx2": 32},
    "time": {"t_end": 0.1}
  },
  "runtime": {"backend": "cpu", "mpi": false, "ranks": 1},
  "result": {
    "out_dir": "runs/disk",
    "output_interval": 0.05,
    "checkpoint_interval": 0.05
  }
}
```

Paths resolve from the execution working directory, not the spec file's directory.
`--out` overrides the result path. A launcher is an argument array, not a shell
command string. Native INI limitations prohibit `#`, `;`, and line breaks in paths
passed to that parser.

## Configuration groups

| Group | Settings |
| --- | --- |
| `grid` | `geometry`, `nx1/nx2`, `radial_grid`, coordinate bounds |
| `physics` | `gamma`, `gm`, `eta_model/eta`, `h/corona`, `mu0/mu_ref/b0`, `bend`, `alpha_m/chi`, `cooling`, `initial_equilibrium`, `wave_amplitude` |
| `numerics` | Density/pressure floors, thermal recovery/sync, dual-energy and entropy thresholds, CT correction, retries, caches, metric PPM, CT EMF scheme, staggered CFL |
| `boundary` | `lower1/upper1/lower2/upper2`, `inner_normal_field` |
| `time` | `cfl`, `t_end`, `max_steps` |

The dataclasses in [config.py](../src/pyhmd/config.py) define supported fields and
choices. Only `resistive_mhd`, `ppm_hll_ssprk3`, and `float64` are registered for
model, scheme, and dtype. Existing options such as `thermal_recovery` remain;
HLLD, 3D, and FP32 are not selectable.

Partial dictionaries merge with the problem preset. Passing a group dataclass
instead supplies that entire group explicitly. Python checks types/basic ranges;
`Engine.plan()` invokes native cross-option checks. Invalid Cartesian axis
boundaries, one-sided periodic boundaries, and unsupported transport combinations
are rejected before launch.

## INI compatibility and recording

INI follows the native flat-key contract: sections do not create namespaces.
Unknown and duplicate settings are rejected; `entropy_switch` and `log_r` aliases
become canonical fields. `RunSpec.load()` also loads output/checkpoint intervals
and restart paths. Omitted `output_interval` and `log_interval` use problem presets.

Transport diagnostics have native restrictions, including the spherical/log
disk, `face_hll`, and `evolve` choices. Select recording before the run: missing
ledgers cannot be reconstructed as measured face-flux data. Keep output intervals
fixed in trajectory/performance comparisons because they clip time steps.

## Phase 0 measurement controls

See [Phase 0](phase0.md) for additional benchmark problems, `wave_mode`, disk-only `noise_amplitude`/`noise_seed`, and `ResultOptions.midplane_diagnostics`. The latter records raw-pressure midplane statistics at existing output/log samples and does not clip additional time steps.
