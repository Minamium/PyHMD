# PyHMD

[English](README.md) | [日本語](README.ja.md)

A finite-volume resistive MHD code with typed Python configuration and a
C++/CUDA time integrator. The current scope is **2D Cartesian and spherical
axisymmetric 2.5D**, FP64, PPM/HLL/SSP-RK3, and constrained transport (CT).

The **Config → Engine → Result** interface follows kineticEQ's separation of
physical configuration, execution policy, and saved results. Native code owns
the complete time integration; Python does not exchange state arrays at every
Runge–Kutta stage. Importing the package or creating a plan does not compile
code, submit a job, or select a GPU.

## Install and run

The CPU build needs Python 3.10+, a C++17 compiler, and GNU Make.

```bash
git clone https://github.com/Minamium/PyHMD.git
cd PyHMD
python3 -m pip install -e '.[viz]'
pyhmd build --backend cpu
pyhmd init --problem alfven --output run.json
pyhmd plan run.json --out runs/alfven_demo
pyhmd run run.json --out runs/alfven_demo
pyhmd inspect runs/alfven_demo
pyhmd validate runs/alfven_demo
pyhmd render runs/alfven_demo
```

Configuration, orchestration, and basic numerical validation use only the Python
standard library. The `analysis` extra adds NumPy; `viz` adds NumPy, Matplotlib,
and Pillow. Wheels include native sources and presets, so an installed package
can also build its native executable explicitly.

`plan` validates and prints the command without creating the result directory.
`run` executes synchronously and writes native output to the result's `run.log`.
Use a new output directory for each run; nonempty directories are rejected.

## Python API

```python
from pyhmd import Config, Engine, ResultOptions

cfg = Config(
    problem="thin_disk",
    grid={"nx1": 64, "nx2": 32},
    time={"t_end": 0.1},
)
result = Engine(
    cfg,
    result_out_dir="runs/disk_demo",
    result_options=ResultOptions(output_interval=0.05, checkpoint_interval=0.05),
).run()

print(result.status, result.complete)
result.validate()  # Numerical integrity, not scientific reproduction.
result.render(xlim=(0, 10), ylim=(0, 10))
```

Partial dictionaries merge into the problem preset and become frozen dataclasses.
Unknown keys and unsupported model/scheme/dtype names are rejected. Native
validation checks cross-option restrictions before execution. Existing INI files
remain usable through `RunSpec.load("configs/thin_disk.ini")` or
`pyhmd run configs/thin_disk.ini --out runs/my_disk`.

## Supported scope

| Area | Implemented |
| --- | --- |
| Physics | Classical compressible MHD, ideal gas, point-mass gravity, constant/standard SF2014/strong resistivity, Ohmic heat removal |
| Geometry | 2D Cartesian and spherical axisymmetric 2.5D, with restrictions on axis/equator boundary placement |
| Numerics | Metric PPM, HLL, explicit adaptive SSP-RK3, 2D flux-function CT |
| Robustness | CT energy correction, selectable thermal recovery, rejected-step rollback/retry, correction ledgers |
| CPU | Single process and MPI decomposition in the first coordinate |
| GPU | Single-GPU CUDA and MPI+CUDA with one rank per GPU, FP64 |
| Storage | CSV, binary checkpoints with configuration/source signatures and checksums, resolved settings and executed binary |
| Analysis | Result reload, Matplotlib, GIF/MP4, VTK, optional discrete transport diagnostics |

3D, relativity, AMR, HLLD, implicit integration, and GPU-aware MPI are not
implemented. Unsupported capabilities are rejected rather than silently
replaced. The structural migration preserves the existing numerical model;
**quantitative reproduction of SF2014 remains unverified**. See
[scientific status](docs/scientific-status.md) for known fidelity gaps.

## CUDA and MPI

CUDA builds use `nvcc`; MPI builds use `mpicxx`/`mpirun`. Follow the target site's
allocation and module rules. The current CUDA+MPI Make target obtains flags from
an Open MPI-compatible `mpicxx --showme` wrapper.

```bash
# Example architecture targets: V100=70, A30=80, GH200=90.
pyhmd build --backend cuda --cuda-sm 70
pyhmd run examples/thin_disk_cuda.json --out runs/gpu_demo

# The Engine launches native ranks; do not wrap Python itself in mpirun.
pyhmd build --backend cuda --mpi --cuda-sm 70
pyhmd run examples/thin_disk_cuda4.json --out runs/gpu4_demo
```

The native layer preserves `CUDA_VISIBLE_DEVICES` and maps local ranks to GPUs
inside that allocation. Decomposition is along the first coordinate; halos and
shared face tables pass through host memory. See the [HPC guide](docs/hpc.md)
for launchers, the PBS template, and direct execution.

`make cpu`, `make cuda-single CUDA_SM=70`, and `make cuda CUDA_SM=70` are also
supported. Canonical executables are `build/pyhmd_*`; the root Makefile also
creates legacy `pymhd_*` aliases. This layout has CPU validation; CUDA/MPI builds
and performance require separate hardware checks.

## Results and restart

`result.json` records resolved configuration, the source fingerprint, binary
SHA256, launch command, execution status, and artifact references. Native
metadata, diagnostics, status, and transport ledgers are preserved. The executed
binary is copied into `inputs/`, separating the run from subsequent rebuilds.
Manifest updates are atomic.

```python
from dataclasses import replace
from pyhmd import Result, RunSpec

result = Result.load("runs/disk_demo")
checkpoint = result.out_dir / result.checkpoints[-1]["path"]
spec = RunSpec.load("examples/thin_disk_cpu.json")
spec = replace(spec, config=replace(spec.config, time=replace(spec.config.time, t_end=0.2)))
continued = spec.with_runtime(restart=str(checkpoint)).engine(out_dir="runs/continued").run()
```

Reaching `max_steps` produces `incomplete`/CLI exit 3; persistent retries produce
`stopped`/4; failures produce `failed`/2. Completion and validator success are
distinct from scientific reproduction. Legacy PyMHD CSV results can be read;
normal restart from an older source generation is rejected by the source
signature. See [result contracts](docs/results.md) before resuming a run.

## Documentation and agent skills

English is the primary documentation language. Japanese counterparts cover the
same supported behavior and limitations. API names, configuration keys, and
machine-readable artifacts remain in English.

| Guide | English | 日本語 |
| --- | --- | --- |
| Documentation index | [Index](docs/index.md) | [一覧](docs/ja/index.md) |
| Architecture | [Modules and ownership](docs/architecture.md) | [設計](docs/ja/architecture.md) |
| kineticEQ design | [Design mapping](docs/kineticeq-design.md) | [設計の対応](docs/ja/kineticeq-design.md) |
| Configuration | [Python / JSON / INI](docs/configuration.md) | [設定](docs/ja/configuration.md) |
| Results | [States and restart](docs/results.md) | [結果と再開](docs/ja/results.md) |
| HPC | [Build and launch](docs/hpc.md) | [HPC実行](docs/ja/hpc.md) |
| Scientific status | [Scope and gaps](docs/scientific-status.md) | [検証範囲と既知の差](docs/ja/scientific-status.md) |
| Validation | [Checks and evidence](docs/validation.md) | [検証](docs/ja/validation.md) |
| Agent skills | [Development and use](docs/skills.md) | [開発用・利用用](docs/ja/skills.md) |

Two focused skills are included in `.agents/skills/`:

- [`$pyhmd-develop`](.agents/skills/pyhmd-develop/SKILL.md): implementation,
  numerical/API contracts, regression checks, packaging, and bilingual docs.
- [`$pyhmd-use`](.agents/skills/pyhmd-use/SKILL.md): supported configurations,
  explicit builds, runs, restart, diagnostics, and offline visualization.

## Development checks

```bash
python3 scripts/check_docs.py
make test
make test-python
python3 scripts/verify_migration.py --reference /path/to/baseline/pymhd_cpu \
  --reference-configs /path/to/frozen/configs \
  --candidate build/pyhmd_cpu --target-grid --report build/migration.json
```

Migration comparison requires a separately preserved baseline binary.
Documentation changes need link/language-pair checks; native/API changes need the
relevant solver/API tests. CPU results do not establish CUDA or MPI correctness.

License: MIT. This release preserves the numerical equations and supported scope
of the original PyMHD implementation.

Phase 0 benchmark definitions, strict-FP controls, passive midplane sampling, and pending GPU gates are described in [Phase 0](docs/phase0.md).
