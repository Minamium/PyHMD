"""Runtime thermal-policy differences must invalidate a CPU/GPU comparison."""
from pathlib import Path
import importlib.util
import tempfile

root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('check_run',root/'scripts/check_run.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
meta={'internal_energy_sync':'modern','dual_energy_eta':.001,'x1_faces':[0,1],'x2_faces':[0,1]}
for key,value in [('internal_energy_sync','legacy'),('dual_energy_eta',.002),
                  ('ct_emf_scheme','corner_hll'),('cfl_staggered',0),('inner_normal_field','project')]:
    changed=dict(meta);changed[key]=value
    try:
        module.compare(meta,Path('/nonexistent-snapshot'),changed,Path('/nonexistent-reference'))
    except module.InvalidRun as e:
        assert key in str(e),str(e)
    else:
        raise AssertionError('comparison accepted changed '+key)
print('PASS CPU/GPU comparison rejects changed thermal synchronization, CT, staggered CFL or magnetic boundary')
