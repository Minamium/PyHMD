#!/bin/bash
# Source from Bash on Rokko. Toolchain defaults match the V100 measurements.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    printf 'Use: source scripts/rokko/env.sh\n' >&2
    exit 2
fi
if ! type module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh || return 2
fi
module purge || return 2
module load python/3.12.1 nvhpc/23.1 || return 2

export LANG=C LC_ALL=C
export PYMHD_ROKKO_GCC_ROOT="${PYMHD_ROKKO_GCC_ROOT:-/apps/gcc/12.3.0}"
export PYMHD_ROKKO_PYTHON_ROOT="${PYMHD_ROKKO_PYTHON_ROOT:-/apps/python/3.12.1}"
export CUDA_HOME="${PYMHD_ROKKO_CUDA_ROOT:-/apps/nvidia/23.1/hpc_sdk/Linux_x86_64/23.1/cuda}"
export PYMHD_ROKKO_MPI_ROOT="${PYMHD_ROKKO_MPI_ROOT:-/apps/nvidia/23.1/hpc_sdk/Linux_x86_64/23.1/comm_libs/mpi}"

export CC="$PYMHD_ROKKO_GCC_ROOT/bin/gcc" CXX="$PYMHD_ROKKO_GCC_ROOT/bin/g++"
export NVCC="$CUDA_HOME/bin/nvcc" CUDAHOSTCXX="$CXX"
export MPICXX="$PYMHD_ROKKO_MPI_ROOT/bin/mpicxx" MPIEXEC="$PYMHD_ROKKO_MPI_ROOT/bin/mpirun"
export OMPI_CC="$CC" OMPI_CXX="$CXX"
export PATH="$PYMHD_ROKKO_GCC_ROOT/bin:$CUDA_HOME/bin:$PYMHD_ROKKO_MPI_ROOT/bin:$PYMHD_ROKKO_PYTHON_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$PYMHD_ROKKO_GCC_ROOT/lib64:$CUDA_HOME/lib64:$PYMHD_ROKKO_MPI_ROOT/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export CUDA_SM=70 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export CUDAFLAGS="${CUDAFLAGS:--O3 -lineinfo -Xptxas=-v -ccbin $CXX}"

for pymhd_rokko_tool in "$CC" "$CXX" "$NVCC" "$MPICXX" "$MPIEXEC" "$PYMHD_ROKKO_PYTHON_ROOT/bin/python3"; do
    if [[ ! -x "$pymhd_rokko_tool" ]]; then
        printf 'Missing Rokko tool: %s\n' "$pymhd_rokko_tool" >&2
        return 2
    fi
done
unset pymhd_rokko_tool
