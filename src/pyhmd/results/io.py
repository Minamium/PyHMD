"""Atomic manifests and lazy references to native artifacts."""
from __future__ import annotations
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile

def jsonable(value):
    if is_dataclass(value): return jsonable(asdict(value))
    if isinstance(value,Path): return str(value)
    if isinstance(value,dict): return {str(k):jsonable(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [jsonable(v) for v in value]
    return value

def utc_now(): return datetime.now(timezone.utc).isoformat()

def write_json(path, value):
    path=Path(path)
    fd,name=tempfile.mkstemp(prefix='.'+path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as f:
            json.dump(jsonable(value),f,indent=2,allow_nan=False); f.write('\n')
        os.replace(name,path)
    finally:
        if os.path.exists(name): os.unlink(name)

def artifact_refs(root):
    root=Path(root);snapshots=[];checkpoints=[]
    for path in sorted(root.glob('snapshot_*.csv')):
        with path.open() as f: first=f.readline()
        clock=dict(item.split('=',1) for item in first.removeprefix('#').split() if '=' in item)
        snapshots.append(dict(path=path.name,time=float(clock['time']),step=int(clock['step'])))
    for path in root.glob('checkpoint_step*'):
        if not (path/'checkpoint.meta').is_file(): continue
        parts=(path/'checkpoint.meta').read_text().splitlines()[1].split()
        checkpoints.append(dict(path=path.name,time=float(parts[1]),step=int(parts[3])))
    return dict(snapshots=snapshots,checkpoints=sorted(checkpoints,key=lambda q:q['step']))
