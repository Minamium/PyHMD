# HPC and native execution

[English](hpc.md) | [日本語](ja/hpc.md) | [Index](index.md)

Choose CPU/CUDA and MPI explicitly at build time. GNU Make is the build system;
CMake is not required.

```bash
make cpu
make cpu-mpi MPICXX=mpicxx
make cuda-single CUDA_SM=70 NVCC=nvcc
make cuda CUDA_SM=70 NVCC=nvcc MPICXX=mpicxx
```

`build_native()` and `pymhd build` compile the same bundled targets. Set `CXX`,
`NVCC`, `MPICXX`, `CXXFLAGS`, or `CUDAFLAGS` for the toolchain. CUDA+MPI currently
uses Open MPI-compatible `mpicxx --showme:compile` / `--showme:link`; other wrappers
need a Makefile adaptation. Example SM targets: V100=70, A30=80, GH200=90.

Builds default to checkout `build/` or a fingerprint-specific XDG cache for a
wheel; override with `--build-dir` or `PYMHD_BUILD_DIR`. PyMHD does not automatically
compile, SSH, or call `qsub`/`sbatch`. Use the allocation permitted by the site.

## Launch native ranks from the Engine

```python
import os
from pymhd import RunSpec

spec = RunSpec.load("examples/thin_disk_cuda4.json")
spec = spec.with_runtime(launcher=("mpirun", "--hostfile", os.environ["PBS_NODEFILE"]))
result = spec.engine(out_dir="runs/my_allocated_job").run()
```

Do not wrap the Engine in `mpirun -np 4 python ...`. It launches native MPI ranks
once and appends `-np N` to its launcher array. Use a compatible MPI launcher;
`srun` belongs on the direct-native path below.

Native code owns MPI lifetime and device mapping, preserving scheduler-provided
`CUDA_VISIBLE_DEVICES`. Each rank uses its local rank inside that allocation.
The first-coordinate decomposition requires at least three owned cells per rank.
Halos/shared face tables are host-staged. GPU-aware MPI and communication/compute
overlap are not implemented.

## Direct native execution

This path writes native output without the Engine's `result.json` or frozen binary.

```bash
# Inside an existing PBS allocation.
mpirun --hostfile "$PBS_NODEFILE" -np 4 build/pymhd_cuda \
  --config configs/thin_disk.ini --t-end 1 \
  --checkpoint-interval 0.1 --output runs/native_example

# Inside an existing Slurm allocation with the site's MPI integration enabled.
srun --ntasks=4 build/pymhd_cuda \
  --config configs/thin_disk.ini --t-end 1 --output runs/slurm_native_example
```

[pbs_gpu4.pbs.example](../examples/pbs_gpu4.pbs.example) is a single-node template.
Supply the site's queue, account/group, and modules. These are launch examples,
not universal allocation requests. Public files contain no personal connections.

This layout has CPU validation. CUDA/MPI builds, correctness, profiling, and
scaling need actual hardware checks. Compare GPU model, ranks, grid, accepted
steps, physical time, retries/corrections, and output settings together.

## Strict-FP experimental controls

`pymhd build --backend cuda --mpi --cuda-sm 70 --fp-mode strict` produces `pymhd_cuda_strict`. Select it explicitly with `--binary`. The default compiler policy stays unchanged. See [Phase 0](phase0.md) for frozen-source controls and the allocated V100×4 suite.

## Rokko PBS

Submit from the PyMHD checkout on Rokko. The site environment loads
`python/3.12.1` and `nvhpc/23.1`, selects GCC 12.3 and Open MPI, and sets the
V100 target (SM 70) and single-threaded host math settings. It preserves the
scheduler's GPU visibility. [Environment script](../scripts/rokko/env.sh)

For a manual build in an allocation, run `source scripts/rokko/env.sh` in
Bash, then `PYTHONPATH=src python3 -m pymhd build --backend cuda --mpi --cuda-sm 70`.

```bash
# Four V100s, short default disk smoke case.
qsub scripts/rokko/disk.pbs

# One V100, with your own JSON/INI specification.
qsub -l select=1:ngpus=1:ncpus=8:mpiprocs=1:mem=64gb \
  -v PYMHD_SPEC=/path/to/run.json scripts/rokko/disk.pbs

# Phase 0, using an extracted frozen source beside its freeze.json and source.tar.
qsub -v PYMHD_FROZEN_SOURCE=/path/to/frozen/source scripts/rokko/phase0_gpu4.pbs
```

The [disk job](../scripts/rokko/disk.pbs) derives its MPI rank count from the
node file. For 2, 4 or 8 GPUs, request one node with matching `ngpus` and
`mpiprocs`. Set grid, end time, and checkpoint cadence in `PYMHD_SPEC`.
By default the job builds in a fresh directory. Set `PYMHD_BINARY` to use a
prebuilt CUDA+MPI executable, or `PYMHD_FP_MODE=strict` to build the strict
variant. For an existing allocation, run
`bash scripts/rokko/run_disk.sh` from the checkout.

The PBS log prints the job directory under `runs/` or `PYMHD_RUN_BASE`.
Disk output lives in `result/` there. The
[Phase 0 job](../scripts/rokko/phase0_gpu4.pbs) writes its current stage to
`suite/progress.json`. Its first V100 test has a known unresolved discrepancy.
[Experiment results](experiment-results.md)
