#!/usr/bin/env python3
"""Compare independently frozen binaries/configs, all accepted state and ledgers."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from pymhd.results.checkpoint import read_checkpoint,compare_checkpoints

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def diagnostics(path):
    with Path(path).open() as stream:rows=list(csv.DictReader(stream))
    return [{k:v for k,v in row.items() if k!='wall_seconds'} for row in rows]

def compare_runs(left,right):
    differences=[]
    a,b=[json.loads((p/'status.json').read_text()) for p in (left,right)]
    for status in (a,b):status.pop('wall_seconds',None)
    if a!=b:differences.append('native_status')
    aa,bb=[sorted(p.name for p in root.glob('snapshot_*.csv')) for root in (left,right)]
    if aa!=bb:differences.append('snapshot_names')
    for name in set(aa)&set(bb):
        if (left/name).read_bytes()!=(right/name).read_bytes():differences.append(name)
    if diagnostics(left/'diagnostics.csv')!=diagnostics(right/'diagnostics.csv'):differences.append('diagnostics_except_wall_seconds')
    aa,bb=[sorted(p.name for p in root.glob('checkpoint_step*')) for root in (left,right)]
    if aa!=bb:differences.append('checkpoint_names')
    for name in set(aa)&set(bb):differences += [name+'/'+p for p in compare_checkpoints(left/name,right/name)]
    return differences

def invoke(binary,config,out,overrides,ranks,launcher,restart=None):
    command=([*launcher,'-np',str(ranks)] if launcher else [])+[str(binary),'--config',str(config),'--output',str(out)]
    for key,value in overrides.items():command+=['--'+key,str(value)]
    if restart is not None:command+=['--restart',str(restart)]
    out.mkdir();(out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    with (out/'run.log').open('w') as stream:process=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
    if process.returncode not in (0,3):raise RuntimeError(f'run failed ({process.returncode}): {out / "run.log"}')
    state=json.loads((out/'status.json').read_text())
    if (process.returncode==0)!=(state.get('reason')=='t_end' and state.get('complete') is True):raise RuntimeError('exit/status mismatch')
    return out/f'checkpoint_step{state["step"]}'

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True);p.add_argument('--candidate',type=Path,default=ROOT/'build/pymhd_cpu')
    p.add_argument('--reference-configs',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    p.add_argument('--work-dir',type=Path);p.add_argument('--target-grid',action='store_true');p.add_argument('--ranks',type=int,default=1)
    p.add_argument('--launcher',nargs='+');p.add_argument('--retry-reference',type=Path);p.add_argument('--retry-candidate',type=Path)
    args=p.parse_args(argv)
    if args.ranks<1 or (args.ranks>1 and not args.launcher):p.error('MPI ranks require a launcher')
    binaries={key:getattr(args,key).resolve() for key in ('reference','candidate')};configs=args.reference_configs.resolve()
    report=dict(schema_version=2,passed=False,rank_count=args.ranks,scope='same-platform/backend migration, not scientific validation',
        binary_sha256={k:sha(v) for k,v in binaries.items()},config_sha256={q.name:sha(q) for q in configs.glob('*.ini')},
        compared=['all U/Psi incl. ghosts','all correction and flux ledgers','transport sidecars','clock','trial counters','snapshots','diagnostics except wall_seconds'],cases=[])
    owner=tempfile.TemporaryDirectory(prefix='pymhd_migration_') if args.work_dir is None else None
    work=Path(owner.name) if owner else args.work_dir.resolve();work.mkdir(parents=True,exist_ok=True)
    cases=[('uniform',32,16,.05,{}),('alfven',64,4,.1,{}),('diffusion',64,4,.05,{}),('brio_wu',64,4,.03,{}),
           ('thin_disk',64,32,.1,{}),('alfven',64,4,.1,{'ct-emf-scheme':'face_hll'}),('thin_disk',64,32,.1,{'transport-diagnostics':'yes'})]
    if args.target_grid:cases.append(('thin_disk',600,128,100,{'max-steps':100,'transport-diagnostics':'yes'}))
    try:
        for index,(case,nx,ny,end,extra) in enumerate(cases):
            settings={'nx1':nx,'nx2':ny,'t-end':end,'output-interval':end/2,'checkpoint-interval':end/2,'log-interval':10,**extra}
            outputs={};finals={}
            for name,binary in binaries.items():
                out=work/f'{index}_{case}_{name}';outputs[name]=out
                finals[name]=invoke(binary,configs/f'{case}.ini',out,settings,args.ranks,args.launcher)
            delta=compare_runs(outputs['reference'],outputs['candidate'])
            item=dict(case=case,grid=[nx,ny],settings=settings,differences=delta,passed=not delta,step=read_checkpoint(finals['candidate'])[0].step)
            if case=='thin_disk' and nx==64:
                for name,binary in binaries.items():
                    checkpoint=invoke(binary,configs/f'{case}.ini',work/f'{index}_{name}_prefix',{**settings,'max-steps':3},args.ranks,args.launcher)
                    final=invoke(binary,configs/f'{case}.ini',work/f'{index}_{name}_resume',settings,args.ranks,args.launcher,restart=checkpoint)
                    delta=compare_checkpoints(finals[name],final);item[name+'_restart_differences']=delta;item['passed'] &= not delta
                delta=compare_runs(work/f'{index}_reference_resume',work/f'{index}_candidate_resume')
                item['resumed_cross_binary_differences']=delta;item['passed'] &= not delta
            report['cases'].append(item);print(json.dumps(item),flush=True)
        if bool(args.retry_reference)!=bool(args.retry_candidate):raise ValueError('both retry drivers are required')
        report['forced_retry_exercised']=False
        if args.retry_reference:
            outputs=[]
            for name,binary in [('reference',args.retry_reference),('candidate',args.retry_candidate)]:
                out=work/('forced_retry_'+name)
                command=([*args.launcher,'-np',str(args.ranks)] if args.launcher else [])+[str(binary.resolve()),str(out)]
                subprocess.run(command,check=True)
                if not all(q.rejected_steps>0 for q in read_checkpoint(out)):raise RuntimeError('retry driver did not reject a trial')
                outputs.append(out)
            report['forced_retry_differences']=compare_checkpoints(*outputs);report['forced_retry_exercised']=True
        report['passed']=all(q['passed'] for q in report['cases']) and not report.get('forced_retry_differences')
        report['migration_retry_coverage_complete']=report['passed'] and report['forced_retry_exercised']
    except Exception as error:report['error']=f'{type(error).__name__}: {error}'
    finally:
        args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(report,indent=2)+'\n')
        if owner:owner.cleanup()
    return 0 if report['passed'] else 1

if __name__=='__main__':raise SystemExit(main())
