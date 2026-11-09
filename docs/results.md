# Results and execution states

[English](results.md) | [日本語](ja/results.md) | [Index](index.md)

| Artifact | Contents |
| --- | --- |
| `result.json` | Schema, config/runtime, native capabilities, binary SHA, status, artifact references |
| `inputs/run.json` | Normalized run specification |
| `inputs/native.ini` | Settings actually passed to native code |
| `inputs/pyhmd_*` | Frozen copy of the executed binary |
| `run.log` | Native/MPI stdout and stderr |
| `metadata.json`, `diagnostics.csv`, `status.json` | Native settings, diagnostics, completion/stop state |
| `snapshot_*.csv` | Cell data, time, and step |
| `checkpoint_step*/` | U/Ψ, clock, correction/face-flux ledgers; transport sidecar when enabled |

The Engine requires an absent or empty directory and claims it by exclusively
creating `result.json`. Later updates use a temporary file and atomic rename in
the same directory. It never appends a new run into nonempty results. On an
interruption/exception it cleans up its own process group and records the outcome.

## Terminal states

| Native outcome | Result status | CLI exit | Meaning |
| --- | --- | --- | --- |
| 0 and expected `t_end` reached | `completed` | 0 | Integration completed; scientific success is separate |
| 3 and `max_steps` | `incomplete` | 3 | Step limit before completion |
| 4 and `persistent_retries` | `stopped` | 4 | Intentional stop after repeated deep retries |
| Other failure or status mismatch | `failed` | 2 | Inspect the log; `SimulationError.result` exposes the result |

Preparation/interruption can also produce `created`, `running`, or `interrupted`.
`Result.load()` reads the manifest or legacy metadata/status output. Optional
array/plotting dependencies load only when requested.

`Result.validate()` is the strict numerical-integrity validator: floor events
cause rejection. It does not certify the physical solution or automatically mark
`scientific_validation` successful. Legacy runs need an explicit expected time,
for example `pyhmd validate runs/legacy --expect-time 1`.

## Restart

Grid, rank decomposition, output interval, physical/numerical settings, and
algorithm/source signatures must match. Set `RuntimeConfig.restart` to the
checkpoint directory and choose a new output directory. Clock-limit fields
excluded from the signature, such as `t_end`/`max_steps`, can extend beyond the
saved clock. The CLI has
no `--restart` switch: use Python or the existing INI `restart` setting.

Binary format version 1, legacy magic, and metadata version 2 are retained, but
the source fingerprint changes with this migration. Normal restart from an older
source generation is rejected. A compatibility conversion would need separate
validation; do not edit or bypass signatures to resume.

See the [README restart example](../README.md) and
[scientific status](scientific-status.md) before interpreting completion as reproduction.
