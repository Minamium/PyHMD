#!/bin/bash
# Run from a single-node allocation with one PBS MPI slot per V100.
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$root"
source scripts/rokko/allocation.sh
ranks=$(pymhd_rokko_ranks)
source "${PYMHD_SITE_ENV:-scripts/rokko/env.sh}"
export PYTHONPATH="$root/src${PYTHONPATH:+:$PYTHONPATH}" PYTHONUNBUFFERED=1

spec=${PYMHD_SPEC:-examples/thin_disk_cuda4.json}
if [[ ! -r "$spec" ]]; then
    printf 'Cannot read run specification: %s\n' "$spec" >&2
    exit 2
fi
binary=${PYMHD_BINARY:-}
if [[ -n "$binary" && ! -x "$binary" ]]; then
    printf 'Cannot execute supplied CUDA/MPI binary: %s\n' "$binary" >&2
    exit 2
fi

run_base=${PYMHD_RUN_BASE:-$root/runs}
mkdir -p "$run_base"
job_dir=$(mktemp -d "$run_base/rokko_disk_${PBS_JOBID%%.*}_XXXXXX")
printf 'Job directory: %s\n' "$job_dir"
pymhd_rokko_verify_v100 "$job_dir" "$ranks"

if [[ -z "$binary" ]]; then
    fp_mode=${PYMHD_FP_MODE:-default}
    python3 -m pymhd build --backend cuda --mpi --cuda-sm 70 \
        --fp-mode "$fp_mode" --build-dir "$job_dir/build"
    binary="$job_dir/build/pymhd_cuda"
    if [[ "$fp_mode" == strict ]]; then
        binary="$job_dir/build/pymhd_cuda_strict"
    fi
fi

exec python3 - "$spec" "$binary" "$ranks" "$job_dir/result" <<'PY'
import os
from pathlib import Path
import sys

from pymhd import RunSpec, SimulationError

spec_path, binary, ranks, out_dir = sys.argv[1:]
spec = RunSpec.load(spec_path).with_runtime(
    backend="cuda",
    mpi=True,
    ranks=int(ranks),
    binary=Path(binary).resolve(),
    launcher=(os.environ["MPIEXEC"], "--hostfile", os.environ["PBS_NODEFILE"], "--bind-to", "core"),
)
try:
    result = spec.engine(out_dir=out_dir).run()
except SimulationError as error:
    print(f"pymhd: {error}", file=sys.stderr)
    raise SystemExit(2) from error
print(f"Status: {result.status}\nResults: {result.out_dir}")
raise SystemExit({"completed": 0, "incomplete": 3, "stopped": 4}.get(result.status, 2))
PY
