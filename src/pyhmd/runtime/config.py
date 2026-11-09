"""Scheduler-neutral execution policy; never chooses GPUs outside the allocation."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from pyhmd.config import _choice, _int

@dataclass(frozen=True)
class RuntimeConfig:
    backend: str = 'cpu'
    mpi: bool = False
    ranks: int = 1
    launcher: tuple[str, ...] = ('mpirun',)
    binary: str | Path | None = None
    build_dir: str | Path | None = None
    restart: str | Path | None = None

    def __post_init__(self):
        _choice(self.backend, ('cpu', 'cuda'), 'backend'); _int(self.ranks, 'ranks', 1)
        if not isinstance(self.mpi, bool): raise ValueError('mpi must be bool')
        if self.ranks != 1 and not self.mpi: raise ValueError('ranks > 1 requires mpi=true')
        if not isinstance(self.launcher, (tuple, list)) or not self.launcher or not all(isinstance(x,str) and x for x in self.launcher):
            raise ValueError('launcher must be a nonempty array of arguments, never a shell command string')
        object.__setattr__(self, 'launcher', tuple(self.launcher))
        for key in ('binary', 'build_dir', 'restart'):
            value = getattr(self,key)
            if value is not None and (not isinstance(value,(str,Path)) or not str(value)):
                raise ValueError(f'{key} must be a nonempty path')
