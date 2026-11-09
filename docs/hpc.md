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

`build_native()` and `pyhmd build` compile the same bundled targets. Set `CXX`,
`NVCC`, `MPICXX`, `CXXFLAGS`, or `CUDAFLAGS` for the toolchain. CUDA+MPI currently
uses Open MPI-compatible `mpicxx --showme:compile` / `--showme:link`; other wrappers
need a Makefile adaptation. Example SM targets: V100=70, A30=80, GH200=90.

Builds default to checkout `build/` or a fingerprint-specific XDG cache for a
wheel; override with `--build-dir` or `PYHMD_BUILD_DIR`. PyHMD does not automatically
compile, SSH, or call `qsub`/`sbatch`. Use the allocation permitted by the site.

## Launch native ranks from the Engine

```python
import os
from pyhmd import RunSpec

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
mpirun --hostfile "$PBS_NODEFILE" -np 4 build/pyhmd_cuda \
  --config configs/thin_disk.ini --t-end 1 \
  --checkpoint-interval 0.1 --output runs/native_example

# Inside an existing Slurm allocation with the site's MPI integration enabled.
srun --ntasks=4 build/pyhmd_cuda \
  --config configs/thin_disk.ini --t-end 1 --output runs/slurm_native_example
```

[pbs_gpu4.pbs.example](../examples/pbs_gpu4.pbs.example) is a single-node template.
Supply the site's queue, account/group, and modules. These are launch examples,
not universal allocation requests. Public files contain no personal connections.

This layout has CPU validation. CUDA/MPI builds, correctness, profiling, and
scaling need actual hardware checks. Compare GPU model, ranks, grid, accepted
steps, physical time, retries/corrections, and output settings together.

## Strict-FP experimental controls

`pyhmd build --backend cuda --mpi --cuda-sm 70 --fp-mode strict` produces `pyhmd_cuda_strict`. Select it explicitly with `--binary`. The default compiler policy stays unchanged. See [Phase 0](phase0.md) for frozen-source controls and the allocated V100×4 suite.
