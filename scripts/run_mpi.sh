#!/usr/bin/env bash
# Scheduler-independent launch; use inside the target cluster's allocation.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RANKS="${PYMHD_RANKS:-${PYHMD_RANKS:-1}}"
BINARY="${PYMHD_BINARY:-${PYHMD_BINARY:-${ROOT}/build/pymhd_cuda}}"
LAUNCHER="${PYMHD_LAUNCHER:-${PYHMD_LAUNCHER:-mpirun}}"
if [[ ! "${RANKS}" =~ ^[1-9][0-9]*$ ]]; then printf 'Invalid PYMHD_RANKS.\n' >&2; exit 2; fi
test -x "${BINARY}"
command -v "${LAUNCHER}" >/dev/null
cd "${ROOT}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
case "${LAUNCHER##*/}" in
    srun) exec "${LAUNCHER}" --ntasks="${RANKS}" "${BINARY}" "$@" ;;
    mpirun|mpiexec) exec "${LAUNCHER}" -n "${RANKS}" "${BINARY}" "$@" ;;
    *) printf 'Launcher must be mpirun, mpiexec, or srun; adapt this script for the site.\n' >&2; exit 2 ;;
esac
