"""Read saved runs without importing plotting libraries or advancing a solver."""
from __future__ import annotations
from dataclasses import dataclass
import json
from pathlib import Path
from .io import artifact_refs

@dataclass(frozen=True)
class Result:
    out_dir: Path

    def __post_init__(self):
        object.__setattr__(self,'out_dir',Path(self.out_dir).expanduser().resolve())
        if not (self.out_dir/'result.json').is_file() and not (self.out_dir/'metadata.json').is_file():
            raise FileNotFoundError(f'no PyMHD result or legacy metadata in {self.out_dir}')

    @classmethod
    def load(cls,path): return cls(Path(path))

    @property
    def manifest(self):
        path=self.out_dir/'result.json'
        if path.is_file(): return json.loads(path.read_text())
        status=self.native_status
        return dict(schema_version=0,legacy=True,status='completed' if status.get('complete') else 'incomplete',
                    scientific_validation='not_evaluated',native_status=status,artifacts=artifact_refs(self.out_dir))

    @property
    def native_status(self):
        path=self.out_dir/'status.json'
        return json.loads(path.read_text()) if path.is_file() else {}

    @property
    def complete(self): return self.status=='completed' and self.native_status.get('complete') is True
    @property
    def status(self): return self.manifest['status']
    @property
    def snapshots(self): return tuple(self.manifest.get('artifacts',{}).get('snapshots',[]))
    @property
    def checkpoints(self): return tuple(self.manifest.get('artifacts',{}).get('checkpoints',[]))

    def snapshot(self,selection='latest'):
        from pymhd.plotting._viz import load_snapshot
        return load_snapshot(self.out_dir,selection)

    def render(self,output_dir=None,**window):
        from pymhd.plotting import render_result
        return render_result(self.out_dir,output_dir=output_dir,**window)

    def validate(self,*,expect_time=None):
        from pymhd.diagnostics.quality import InvalidRun,validate
        expected=expect_time
        if expected is None:
            expected=self.manifest.get('config',{}).get('time',{}).get('t_end')
        if expected is None: raise ValueError('expect_time is required for legacy runs')
        try: _,_,summary=validate(self.out_dir,expected)
        except InvalidRun as error: raise ValueError(str(error)) from error
        return summary
