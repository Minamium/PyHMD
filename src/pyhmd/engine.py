"""Config -> validated native plan -> execution -> reloadable Result."""
from __future__ import annotations
from dataclasses import asdict, dataclass, replace
import json
import os
from pathlib import Path
import signal
import shutil
import subprocess
import tempfile
from pyhmd.config import Config, ResultOptions, _typed, native_ini
from pyhmd.runtime.config import RuntimeConfig
from pyhmd.backends.registry import get_backend
from pyhmd.backends.build import default_build_dir, describe, sha256
from pyhmd.results.io import artifact_refs, jsonable, utc_now, write_json
from pyhmd.results.result import Result

@dataclass(frozen=True)
class RunPlan:
    command: tuple[str,...]
    output: Path
    native_settings: dict
    binary_description: dict
    binary_sha256: str

    def as_dict(self): return jsonable(asdict(self))

class SimulationError(RuntimeError):
    def __init__(self,message,result):
        super().__init__(message);self.result=result

def _stop_owned_process(process):
    if process.poll() is not None: return
    try:
        if os.name=='posix': os.killpg(process.pid,signal.SIGTERM)
        else: process.terminate()
    except ProcessLookupError:
        process.wait()
        return
    try: process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            if os.name=='posix': os.killpg(process.pid,signal.SIGKILL)
            else: process.kill()
        except ProcessLookupError:
            pass
        process.wait()

class Engine:
    def __init__(self,config:Config,*,runtime=None,result_options=None,result_out_dir=None):
        if not isinstance(config,Config): raise TypeError('config must be Config')
        self.config=config
        self.runtime=_typed(RuntimeConfig,runtime)
        options=_typed(ResultOptions,result_options)
        if result_out_dir is not None: options=replace(options,out_dir=result_out_dir)
        self.result_options=options.resolved(config.problem)

    def plan(self):
        runtime=self.runtime
        if any(int(os.environ.get(key,'1'))>1 for key in ('OMPI_COMM_WORLD_SIZE','PMI_SIZE','PMIX_SIZE')):
            raise RuntimeError('run the Python Engine outside MPI; it launches native ranks. Or invoke the native binary directly under MPI.')
        if runtime.ranks > self.config.grid.nx1//3:
            raise ValueError('radial decomposition requires at least three owned cells per rank')
        if self.result_options.out_dir is None: raise ValueError('result_out_dir or result_options.out_dir is required')
        output=Path(self.result_options.out_dir).expanduser().resolve()
        if output.exists() and (not output.is_dir() or any(output.iterdir())):
            raise FileExistsError(f'output must be absent or empty: {output}')
        backend=get_backend(runtime.backend,runtime.mpi)
        directory=Path(runtime.build_dir).expanduser().resolve() if runtime.build_dir is not None else default_build_dir()
        binary=Path(runtime.binary).expanduser().resolve() if runtime.binary is not None else directory/backend.executable
        if not binary.is_file(): raise FileNotFoundError(f'native binary missing: {binary}; build explicitly with pyhmd build')
        info=describe(binary)
        if info['backend']!=runtime.backend or info['mpi']!=runtime.mpi:
            raise ValueError('native binary backend/MPI capabilities do not match RuntimeConfig')
        if info.get('dtype')!=self.config.dtype or info.get('model')!=self.config.model or info.get('scheme')!=self.config.scheme:
            raise ValueError('native model/scheme/dtype does not match Config')
        if self.config.grid.geometry not in info.get('geometries',[]): raise ValueError('unsupported native geometry')
        settings=self.config.native_settings();options=asdict(self.result_options);options.pop('out_dir')
        settings.update(options,output=str(output))
        if runtime.restart is not None:
            checkpoint=Path(runtime.restart).expanduser().resolve()
            if not (checkpoint/'checkpoint.meta').is_file(): raise FileNotFoundError(f'checkpoint manifest missing: {checkpoint}')
            settings['restart']=str(checkpoint)
        # Native validation remains authoritative for cross-option restrictions.
        with tempfile.TemporaryDirectory(prefix='pyhmd_validate_') as folder:
            ini=Path(folder)/'native.ini';ini.write_text(native_ini(settings))
            process=subprocess.run([str(binary),'--validate-config','--config',str(ini)],text=True,capture_output=True)
            if process.returncode: raise ValueError('native configuration rejected: '+process.stderr.strip())
        command=([*runtime.launcher,'-np',str(runtime.ranks)] if runtime.mpi else [])
        command += [str(binary),'--config',str(output/'inputs/native.ini')]
        return RunPlan(tuple(command),output,settings,info,sha256(binary))

    def run(self):
        plan=self.plan();output=plan.output
        output.mkdir(parents=True,exist_ok=True)
        # O_EXCL is the ownership claim even if another process races mkdir.
        try: fd=os.open(output/'result.json',os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o644)
        except FileExistsError as e: raise FileExistsError(f'result already owned: {output}') from e
        manifest=dict(schema_version=1,application='PyHMD',status='created',created_at=utc_now(),
                      config=self.config.as_dict(),runtime=jsonable(self.runtime),result_options=jsonable(self.result_options),
                      execution=plan.as_dict(),scientific_validation='not_evaluated',artifacts={})
        with os.fdopen(fd,'w') as f: json.dump(manifest,f,indent=2,allow_nan=False)
        process=None
        try:
            (output/'inputs').mkdir()
            (output/'inputs/native.ini').write_text(native_ini(plan.native_settings))
            write_json(output/'inputs/run.json',dict(schema_version=1,config=self.config.as_dict(),runtime=jsonable(self.runtime),result=jsonable(self.result_options)))
            binary_index=len(plan.command)-3
            frozen=output/'inputs'/Path(plan.command[binary_index]).name
            shutil.copy2(plan.command[binary_index],frozen)
            if sha256(frozen)!=plan.binary_sha256: raise RuntimeError('native binary changed during run preparation')
            command=list(plan.command);command[binary_index]=str(frozen)
            manifest['execution']['executed_command']=command
            manifest.update(status='running',started_at=utc_now())
            write_json(output/'result.json',manifest)
            with (output/'run.log').open('w') as log:
                process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=os.name=='posix')
                code=process.wait()
            state=json.loads((output/'status.json').read_text()) if (output/'status.json').is_file() else {}
            valid_complete=code==0 and state.get('complete') is True and state.get('reason')=='t_end' and state.get('time')==self.config.time.t_end
            status='completed' if valid_complete else ('incomplete' if code==3 and state.get('reason')=='max_steps' else ('stopped' if code==4 and state.get('reason')=='persistent_retries' else 'failed'))
            manifest.update(status=status,returncode=code,native_status=state,finished_at=utc_now(),artifacts=artifact_refs(output))
            write_json(output/'result.json',manifest)
            result=Result(output)
            if status=='failed': raise SimulationError(f'native run failed ({code}); see {output / "run.log"}',result)
            return result
        except BaseException as error:
            if process is not None: _stop_owned_process(process)
            if manifest['status'] not in ('completed','incomplete','stopped','failed'):
                manifest.update(status='interrupted' if isinstance(error,KeyboardInterrupt) else 'failed',error=f'{type(error).__name__}: {error}',finished_at=utc_now())
                write_json(output/'result.json',manifest)
            raise
