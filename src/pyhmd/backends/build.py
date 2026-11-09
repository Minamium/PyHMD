"""Explicit native builds from bundled source, with inspectable provenance."""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from pyhmd.native.fingerprint import source_fingerprint
from .registry import get_backend

NATIVE = Path(__file__).resolve().parents[1]/'native'

def default_build_dir():
    if os.environ.get('PYHMD_BUILD_DIR'): return Path(os.environ['PYHMD_BUILD_DIR']).expanduser().resolve()
    checkout = NATIVE.parents[2]
    if (checkout/'pyproject.toml').is_file(): return checkout/'build'
    cache = Path(os.environ.get('XDG_CACHE_HOME',Path.home()/'.cache'))
    return cache/'pyhmd'/'native'/source_fingerprint(NATIVE)[:16]

def sha256(path):
    with Path(path).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest() if sys.version_info >= (3,11) else hashlib.sha256(f.read()).hexdigest()

def describe(binary):
    p = subprocess.run([str(Path(binary).resolve()),'--describe'],text=True,capture_output=True,check=True)
    q = json.loads(p.stdout)
    if q.get('protocol_version') != 1 or q.get('application') != 'PyHMD': raise ValueError('unsupported native protocol')
    return q

@dataclass(frozen=True)
class BuildArtifact:
    binary: Path
    binary_sha256: str
    source_fingerprint: str
    description: dict
    manifest: Path

def build_native(backend='cpu', *, mpi=False, build_dir=None, cuda_sm='70', jobs=1, fp_mode='default'):
    target = get_backend(backend,mpi)
    if isinstance(jobs,bool) or not isinstance(jobs,int) or jobs < 1: raise ValueError('jobs must be positive integer')
    if not str(cuda_sm).isdigit(): raise ValueError('cuda_sm must contain digits such as 70 or 90')
    if fp_mode not in ('default','strict'): raise ValueError('fp_mode must be default or strict')
    directory = Path(build_dir).expanduser().resolve() if build_dir is not None else default_build_dir()
    directory.mkdir(parents=True,exist_ok=True)
    lock = directory/'.build.lock'
    try: fd = os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    except FileExistsError as e: raise RuntimeError(f'build is already locked: {lock}; inspect its owning process') from e
    try:
        with os.fdopen(fd,'w') as f: f.write(str(os.getpid())+'\n')
        executable=target.executable+('_strict' if fp_mode=='strict' else '')
        command = ['make','-B','-C',str(NATIVE),f'-j{jobs}',target.make_target,f'BUILD_DIR={directory}',f'CUDA_SM={cuda_sm}',f'FP_MODE={fp_mode}']
        log = directory/(executable+'.build.log')
        with log.open('w') as stream: process = subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
        if process.returncode: raise RuntimeError(f'native build failed ({process.returncode}); see {log}')
        binary = directory/executable; info = describe(binary)
        fingerprint = source_fingerprint(NATIVE)
        if info['source_fingerprint'] != fingerprint: raise RuntimeError('built source fingerprint differs')
        manifest = directory/(executable+'.build.json')
        q = dict(binary=str(binary),binary_sha256=sha256(binary),source_fingerprint=fingerprint,
                 description=info,command=command,fp_mode=fp_mode,cuda_sm=str(cuda_sm) if backend=='cuda' else None,
                 toolchain={name:os.environ.get(name) for name in ('CXX','NVCC','MPICXX','CXXFLAGS','CUDAFLAGS')},log=str(log))
        manifest.write_text(json.dumps(q,indent=2)+'\n')
        return BuildArtifact(binary,q['binary_sha256'],fingerprint,info,manifest)
    finally:
        lock.unlink()
