"""One strict run specification shared by Python and the scheduler-neutral CLI."""
from __future__ import annotations
from dataclasses import dataclass, field, replace
import json
from pathlib import Path
from pyhmd.config import Config, ResultOptions, _typed, config_from_ini
from pyhmd.runtime.config import RuntimeConfig
from pyhmd.engine import Engine
from pyhmd.results.io import jsonable

def _json_unique(pairs):
    result={}
    for key,value in pairs:
        if key in result: raise ValueError(f'duplicate JSON key: {key}')
        result[key]=value
    return result

@dataclass(frozen=True)
class RunSpec:
    config: Config = field(default_factory=Config)
    runtime: RuntimeConfig = field(default_factory=RuntimeConfig)
    result: ResultOptions = field(default_factory=ResultOptions)
    schema_version: int = 1

    def __post_init__(self):
        if type(self.schema_version) is not int or self.schema_version!=1: raise ValueError('unsupported run schema_version')
        object.__setattr__(self,'config',Config.from_dict(self.config) if not isinstance(self.config,Config) else self.config)
        object.__setattr__(self,'runtime',_typed(RuntimeConfig,self.runtime))
        object.__setattr__(self,'result',_typed(ResultOptions,self.result))

    @classmethod
    def from_dict(cls,value): return _typed(cls,value)
    @classmethod
    def load(cls,path):
        path=Path(path)
        if path.suffix.lower()=='.ini':
            config,options,restart=config_from_ini(path)
            return cls(config,RuntimeConfig(restart=restart),options)
        value=json.loads(path.read_text(),object_pairs_hook=_json_unique,
                         parse_constant=lambda x: (_ for _ in ()).throw(ValueError(f'nonfinite JSON number: {x}')))
        return cls.from_dict(value)

    def as_dict(self): return jsonable(self)
    def engine(self,*,out_dir=None):
        return Engine(self.config,runtime=self.runtime,result_options=self.result,result_out_dir=out_dir)

    def with_runtime(self,**overrides): return replace(self,runtime=replace(self.runtime,**overrides))
