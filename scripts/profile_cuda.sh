#!/usr/bin/env bash
# Manual, single-GPU profiling inside an existing allocation.
# This script never submits a job or allocates nodes.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
MODE="${1:-timing}"
DESTINATION="${2:-runs/profiling}"
CONFIG="${3:-configs/thin_disk.ini}"
BINARY="${PYMHD_BINARY:-${PYHMD_BINARY:-build/pymhd_cuda_single}}"
STEPS="${PYMHD_PROFILE_STEPS:-${PYHMD_PROFILE_STEPS:-200}}"
if [[ "$MODE" != timing && "$MODE" != systems && "$MODE" != compute ]]; then
    printf 'Usage: bash scripts/profile_cuda.sh [timing|systems|compute] [OUTPUT_PARENT] [CONFIG]\n' >&2
    exit 2
fi
[[ -x "$BINARY" ]] || { printf 'CUDA binary not found: %s\n' "$BINARY" >&2; exit 2; }
[[ "$STEPS" =~ ^[0-9]+$ && "$STEPS" -ge 100 ]] || { printf 'PYMHD_PROFILE_STEPS must be at least 100.\n' >&2; exit 2; }
mkdir -p "$DESTINATION"
PROFILE_DIR="$(mktemp -d "${DESTINATION%/}/${MODE}_XXXXXX")"
printf 'Profile directory: %s\n' "$PROFILE_DIR"
command -v nvidia-smi >/dev/null
nvidia-smi -q > "$PROFILE_DIR/gpu.txt"
if command -v sha256sum >/dev/null; then sha256sum "$BINARY" "$CONFIG" > "$PROFILE_DIR/input_hashes.txt"; fi
LAUNCH=()
if [[ -n "${PBS_NODEFILE:-}" ]]; then LAUNCH=(mpirun --hostfile "$PBS_NODEFILE" -np 1); fi

run_profile() {
    local block="$1" destination="$2";shift 2
    mkdir -p "$destination"
    set +e
    PYMHD_GPU_BLOCK_SIZE="$block" "${LAUNCH[@]}" "$@" "$BINARY" \
        --config "$CONFIG" --t-end 1000000 --max-steps "$STEPS" \
        --output-interval 1000000 --log-interval 50 --output "$destination" \
        > "$destination/execution.log" 2>&1
    local result=$?
    set -e
    # A deliberately fixed-step run exits 3 (max_steps), not 0 (t_end).
    if [[ "$result" != 0 && "$result" != 3 ]]; then
        printf 'Profiling failed (%s): %s/execution.log\n' "$result" "$destination" >&2
        return "$result"
    fi
    python3 - "$destination" "$STEPS" <<'PY'
import json, sys
from pathlib import Path
p=Path(sys.argv[1]);s=json.loads((p/'status.json').read_text())
if s['reason'] != 'max_steps' or s['step'] != int(sys.argv[2]):
    raise SystemExit(f'Unexpected run termination: {s}')
PY
}

case "$MODE" in
timing)
    for block in ${PYMHD_PROFILE_BLOCKS:-64 128 256}; do
        for repetition in 1 2 3; do
            run_profile "$block" "$PROFILE_DIR/block${block}_repeat${repetition}"
        done
    done
    python3 - "$PROFILE_DIR" <<'PY'
import csv, json, re, statistics, sys
from pathlib import Path
root=Path(sys.argv[1]);rows=[]
for d in sorted(root.glob('block*_repeat*')):
    block=int(re.fullmatch(r'block(\d+)_repeat\d+', d.name).group(1))
    data=list(csv.DictReader((d/'diagnostics.csv').open()))
    # Exclude initialization and first 50 steps; compare identical step counts.
    a=next(r for r in data if int(r['step'])==50);b=data[-1]
    steps=int(b['step'])-50;seconds=float(b['wall_seconds'])-float(a['wall_seconds'])
    meta=json.loads((d/'metadata.json').read_text())
    row=dict(block=block,run=d.name,measured_steps=steps,seconds=seconds,
             milliseconds_per_step=1000*seconds/steps,
             million_cell_updates_per_second=steps*meta['nx1']*meta['nx2']/seconds/1e6)
    rows.append(row)
summary=[]
for block in sorted({r['block'] for r in rows}):
    group=[r['milliseconds_per_step'] for r in rows if r['block']==block]
    summary.append(dict(block=block,median_ms_per_step=statistics.median(group),min_ms_per_step=min(group),max_ms_per_step=max(group)))
(root/'timing.json').write_text(json.dumps(dict(runs=rows,summary=summary),indent=2)+'\n')
for r in summary:print(f"block={r['block']:3} median={r['median_ms_per_step']:.4f} ms/step range={r['min_ms_per_step']:.4f}..{r['max_ms_per_step']:.4f}")
PY
    ;;
systems)
    command -v nsys >/dev/null
    run_profile "${PYMHD_GPU_BLOCK_SIZE:-128}" "$PROFILE_DIR/run" \
        nsys profile --trace=cuda,mpi --mpi-impl=openmpi --sample=none --cpuctxsw=none \
        --output "$PROFILE_DIR/timeline"
    ;;
compute)
    command -v ncu >/dev/null
    run_profile "${PYMHD_GPU_BLOCK_SIZE:-128}" "$PROFILE_DIR/run" \
        ncu --section SpeedOfLight --section Occupancy --section LaunchStats \
        --section MemoryWorkloadAnalysis --launch-skip 100 --launch-count 40 \
        --export "$PROFILE_DIR/kernels"
    ;;
esac
