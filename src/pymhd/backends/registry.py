"""Only real native backends are registered; no fallback or speculative plugins."""
from dataclasses import dataclass

@dataclass(frozen=True)
class NativeBackend:
    make_target: str
    executable: str

BACKENDS = {
    ('cpu',False): NativeBackend('cpu','pymhd_cpu'),
    ('cpu',True): NativeBackend('cpu-mpi','pymhd_cpu_mpi'),
    ('cuda',False): NativeBackend('cuda-single','pymhd_cuda_single'),
    ('cuda',True): NativeBackend('cuda','pymhd_cuda'),
}

def get_backend(backend, mpi=False):
    try: return BACKENDS[backend,mpi]
    except KeyError as e: raise ValueError(f'unsupported backend: {(backend,mpi)}') from e
