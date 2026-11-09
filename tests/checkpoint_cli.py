#!/usr/bin/env python3
"""Verify actual CLI interruption/restart follows the uninterrupted trajectory."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile

parser=argparse.ArgumentParser()
parser.add_argument('--binary',type=Path,default=Path('build/pymhd_cpu'))
args=parser.parse_args();binary=args.binary.resolve()
with tempfile.TemporaryDirectory(prefix='pymhd_checkpoint_cli_') as tmp:
    root=Path(tmp)
    common=[str(binary),'--case','alfven','--nx1','32','--nx2','4','--cfl','0.4',
            '--t-end','0.1','--output-interval','0.05','--checkpoint-interval','0.02']
    def run(name,extra=(),expected=0):
        options=dict(zip(common[1::2],common[2::2]));options.update(zip(extra[::2],extra[1::2]));options['--output']=str(root/name)
        command=[str(binary)]+[value for pair in options.items() for value in pair]
        result=subprocess.run(command,capture_output=True,text=True)
        if result.returncode!=expected:raise AssertionError(result.stdout+result.stderr)
        return root/name
    def final_checkpoint(run):
        return max(run.glob('checkpoint_step*'),key=lambda p:int(p.name.removeprefix('checkpoint_step')))
    def same_checkpoint(a,b):
        # Includes conservative/vertex bits, accepted thermal/flux ledgers,
        # counters, output schedule and the last accepted dt, not wall time.
        assert (final_checkpoint(a)/'rank000000.bin').read_bytes()==(final_checkpoint(b)/'rank000000.bin').read_bytes(), 'Restart changed the complete final binary checkpoint'
    full=run('full');split=run('split',['--max-steps','5'],3)
    status=json.loads((split/'status.json').read_text())
    assert status['reason']=='max_steps' and status['step']==5
    checkpoint=split/'checkpoint_step5'
    resumed=run('resumed',['--restart',str(checkpoint)])
    a=sorted(full.glob('snapshot_*.csv'))[-1].read_bytes()
    b=sorted(resumed.glob('snapshot_*.csv'))[-1].read_bytes()
    assert a==b,'Restart changed the final accepted snapshot'
    same_checkpoint(full,resumed)
    assert json.loads((resumed/'metadata.json').read_text())['restarted']
    assert json.loads((resumed/'status.json').read_text())['complete']
    run('bad',['--restart',str(checkpoint),'--cfl','0.3'],2)
    # Stopping between scheduled outputs must retain the next output time.
    half=run('half',['--t-end','0.03'])
    half_checkpoint=final_checkpoint(half)
    continuation=run('continued',['--restart',str(half_checkpoint)])
    times=[float(path.open().readline().split('time=')[1].split()[0]) for path in continuation.glob('snapshot_*.csv')]
    assert 0.05 in times,'A final snapshot before an output boundary skipped the next scheduled output after restart'
    # Cartesian radius() does not read mesh coordinates. Exercise the real
    # spherical disk too, including CUDA's device-resident geometry/cache.
    common=[str(binary),'--config',str(Path(__file__).resolve().parents[1]/'configs/thin_disk.ini'),
            '--nx1','64','--nx2','32','--t-end','0.04','--output-interval','0.01',
            '--checkpoint-interval','0.01']
    disk_full=run('disk_full')
    disk_split=run('disk_split',['--max-steps','1'],3)
    disk_resumed=run('disk_resumed',['--restart',str(disk_split/'checkpoint_step1')])
    assert sorted(disk_full.glob('snapshot_*.csv'))[-1].read_bytes()==sorted(disk_resumed.glob('snapshot_*.csv'))[-1].read_bytes(), 'Spherical restart changed the final accepted snapshot'
    same_checkpoint(disk_full,disk_resumed)
print('PASS CLI --cfl, exact Cartesian/spherical split continuation, changed-physics rejection, and between-output restart scheduling')
