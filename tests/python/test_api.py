"""End-to-end contracts at the Python/native boundary, not formula replicas."""
import dataclasses
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from pymhd import Config,Engine,Result,ResultOptions,RunSpec,RuntimeConfig,SimulationError

ROOT=Path(__file__).resolve().parents[2]
BINARY=ROOT/'build/pymhd_cpu'

class ConfigContracts(unittest.TestCase):
    def test_unknown_and_unsupported_are_rejected(self):
        for value in ({'model':'rmhd'},{'dtype':'float32'},{'scheme':'hlld'},
                      {'grid':{'nx3':16}},{'physics':{'alpha':.1}},{'typo':1},
                      {'grid':{'nx1':True}},{'time':{'cfl':float('nan')}}):
            with self.subTest(value=value),self.assertRaises((ValueError,TypeError)):
                Config.from_dict(value)

    def test_immutable_and_problem_defaults(self):
        cfg=Config(problem='alfven',grid={'nx1':32})
        self.assertEqual(cfg.grid.nx2,4)
        self.assertEqual(cfg.physics.eta_model,'constant')
        self.assertEqual(cfg.numerics.thermal_recovery,'off')
        with self.assertRaises(dataclasses.FrozenInstanceError):cfg.grid.nx1=50
        with self.assertRaises(ValueError):RuntimeConfig(ranks=2)
        with self.assertRaises(ValueError):RuntimeConfig(launcher='mpirun -np 4')

    def test_json_roundtrip_duplicate_and_nonfinite(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'run.json';spec=RunSpec(Config(problem='diffusion'))
            path.write_text(json.dumps(spec.as_dict()))
            self.assertEqual(RunSpec.load(path),spec)
            for bad in ('{"config":{},"config":{}}','{"schema_version":true}','{"config":{"time":{"cfl":NaN}}}'):
                path.write_text(bad)
                with self.assertRaises(ValueError):RunSpec.load(path)

@unittest.skipUnless(BINARY.exists(),'build the CPU native engine first')
class NativeIntegration(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='pymhd_api_');self.root=Path(self.temp.name)
        self.cfg=Config(problem='alfven',grid={'nx1':32},time={'t_end':.1})

    def tearDown(self):self.temp.cleanup()
    def engine(self,name,config=None,**runtime):
        return Engine(config or self.cfg,runtime=RuntimeConfig(binary=BINARY,**runtime),
                      result_options=ResultOptions(out_dir=self.root/name,checkpoint_interval=.02))

    def test_dry_run_has_no_output_and_native_combinations_reject(self):
        plan=self.engine('plan').plan()
        self.assertFalse(plan.output.exists())
        self.assertEqual(plan.binary_description['backend'],'cpu')
        bad=Config(problem='alfven',boundary={'lower1':'periodic','upper1':'outflow'})
        with self.assertRaises(ValueError):self.engine('bad',bad).plan()
        with self.assertRaises(ValueError):self.engine('badgpu',backend='cuda').plan()
        self.assertFalse((self.root/'badgpu').exists())

    def test_run_manifest_ownership_reload_and_exact_restart(self):
        full=self.engine('full').run()
        self.assertTrue(full.complete)
        self.assertEqual(full.validate()['validation_scope'],'numerical_integrity_only')
        self.assertEqual(Result.load(full.out_dir).manifest['scientific_validation'],'not_evaluated')
        self.assertTrue((full.out_dir/'inputs/pymhd_cpu').is_file())
        with self.assertRaises(FileExistsError):self.engine('full').run()
        split=self.engine('split',dataclasses.replace(self.cfg,time=dataclasses.replace(self.cfg.time,max_steps=3))).run()
        self.assertFalse(split.complete);self.assertEqual(split.status,'incomplete')
        resumed=self.engine('resume',restart=split.out_dir/split.checkpoints[-1]['path']).run()
        for key in ('snapshots','checkpoints'):
            a=getattr(full,key)[-1];b=getattr(resumed,key)[-1]
            pa=full.out_dir/a['path'];pb=resumed.out_dir/b['path']
            if key=='checkpoints':pa=pa/'rank000000.bin';pb=pb/'rank000000.bin'
            self.assertEqual(pa.read_bytes(),pb.read_bytes())
        self.assertEqual(full.validate()['time'],resumed.validate()['time'])

    def test_failed_restart_has_an_explicit_failed_result(self):
        split=self.engine('split',dataclasses.replace(self.cfg,time=dataclasses.replace(self.cfg.time,max_steps=3))).run()
        changed=dataclasses.replace(self.cfg,time=dataclasses.replace(self.cfg.time,cfl=.3))
        with self.assertRaises(SimulationError) as failure:
            self.engine('bad',changed,restart=split.out_dir/split.checkpoints[-1]['path']).run()
        self.assertEqual(failure.exception.result.status,'failed')
        self.assertFalse(failure.exception.result.complete)

    def test_all_existing_presets_validate_and_cli_reports_partial(self):
        for ini in (ROOT/'configs').glob('*.ini'):
            spec=RunSpec.load(ini)
            spec=spec.with_runtime(binary=BINARY)
            spec.engine(out_dir=self.root/ini.stem).plan()
        spec=RunSpec(dataclasses.replace(self.cfg,time=dataclasses.replace(self.cfg.time,max_steps=3)),
                     RuntimeConfig(binary=str(BINARY)),ResultOptions(checkpoint_interval=.02))
        path=self.root/'run.json';path.write_text(json.dumps(spec.as_dict()))
        p=subprocess.run([sys.executable,'-m','pymhd','run',str(path),'--out',str(self.root/'cli')],capture_output=True,text=True)
        self.assertEqual(p.returncode,3,p.stderr)
        self.assertFalse(json.loads(p.stdout)['complete'])

if __name__=='__main__':unittest.main()
