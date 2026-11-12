#!/bin/bash
# Helpers for single-node V100 PBS jobs. Source from the PyMHD checkout.
pymhd_rokko_ranks() {
    if [[ -z "${PBS_JOBID:-}" || ! -r "${PBS_NODEFILE:-}" ]]; then
        printf 'An existing PBS allocation is required (PBS_JOBID and PBS_NODEFILE).\n' >&2
        return 2
    fi

    local host
    host=$(hostname -s) || return
    awk -v current="$host" '
        NF {
            split($1, parts, ".");
            sub(/\[[0-9]+\]$/, "", parts[1]);
            if (parts[1] != current) other = 1;
            count++;
        }
        END {
            if (other || count < 1 || count > 8) {
                print "Expected one local physical node with 1 to 8 PBS MPI slots." > "/dev/stderr";
                exit 2;
            }
            print count;
        }
    ' "$PBS_NODEFILE"
}

pymhd_rokko_verify_v100() {
    local job_dir=$1 ranks=$2 verified_mask
    mkdir -p "$job_dir/bin" "$job_dir/p$ranks" || return
    "$NVCC" -O2 -ccbin "$CXX" -arch=sm_70 tests/phase0/device_probe.cpp \
        -o "$job_dir/bin/gpu_probe" || return
    verified_mask=$(python3 tests/phase0/v100_allocation.py "$job_dir" "$ranks") || return
    if [[ -n "$verified_mask" ]]; then
        export CUDA_VISIBLE_DEVICES="$verified_mask"
    fi
}
