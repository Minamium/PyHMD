#!/usr/bin/env python3
"""Identify only scheduler-assigned CUDA devices; never choose free GPUs."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys

experiment, rank_text = sys.argv[1:]
ranks = int(rank_text)
root = Path(experiment)
out = root/f'p{ranks}'
original = os.environ.get('CUDA_VISIBLE_DEVICES')
gpu_file = os.environ.get('PBS_GPUFILE')
evidence = dict(original_cuda_visible_devices=original, pbs_gpufile=gpu_file,
                requested_gpus=ranks, pbs_job_id=os.environ['PBS_JOBID'])
environment = dict(os.environ)
if not original and gpu_file and Path(gpu_file).is_file():
    entries = Path(gpu_file).read_text().split()
    evidence['pbs_gpufile_entries'] = entries
    numbers = [re.search(r'(?:-gpu|/gpu|gpu)(\d+)$', s) for s in entries]
    if len(entries) == ranks and all(numbers):
        environment['CUDA_VISIBLE_DEVICES'] = ','.join(m.group(1) for m in numbers)
        evidence['selection_source'] = 'PBS_GPUFILE'
result = subprocess.run([str(root/'bin/gpu_probe')], env=environment,
                        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
evidence['probe_returncode'] = result.returncode
evidence['probe_stderr'] = result.stderr
try:
    probe = json.loads(result.stdout)
    evidence['probe'] = probe
    if result.returncode or probe['count'] != ranks:
        raise RuntimeError('Visible GPU count does not equal requested allocation; refusing GPU computation')
    if not all(d['name'].startswith('Tesla V100') or d['name'].startswith('NVIDIA V100') or
               d['name'].startswith('V100') for d in probe['devices']):
        raise RuntimeError('Unexpected GPU model')
    if not all(d['major'] == 7 and d['minor'] == 0 for d in probe['devices']):
        raise RuntimeError('Unexpected GPU compute capability')
    evidence.setdefault('selection_source', 'scheduler CUDA mask' if original else 'CUDA runtime restricted visibility')
    evidence['passed'] = True
    evidence['cuda_visible_devices'] = environment.get('CUDA_VISIBLE_DEVICES')
except Exception as error:
    evidence.update(passed=False, error=str(error), probe_stdout=result.stdout)
(out/'gpu_allocation.json').write_text(json.dumps(evidence, indent=2)+'\n')
if not evidence['passed']:
    raise SystemExit(evidence['error'])
# stdout is only the verified mask; diagnostics stay in the evidence file.
print(evidence['cuda_visible_devices'] or '')
