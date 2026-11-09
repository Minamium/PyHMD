"""Shared Python/JSON entry points; numerical stepping stays in the native engine."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import subprocess
from pyhmd import Config,RunSpec,Result,SimulationError,__version__
from pyhmd.backends.build import build_native
from pyhmd.results.io import jsonable

def _print(value): print(json.dumps(jsonable(value),indent=2,allow_nan=False))

def main(argv=None):
    parser=argparse.ArgumentParser(prog='pyhmd')
    parser.add_argument('--version',action='version',version=__version__)
    commands=parser.add_subparsers(dest='command',required=True)
    for name in ('run','plan'):
        p=commands.add_parser(name);p.add_argument('spec',type=Path);p.add_argument('--out',type=Path)
        p.add_argument('--backend',choices=('cpu','cuda'));p.add_argument('--mpi',action='store_true',default=None)
        p.add_argument('--ranks',type=int);p.add_argument('--binary',type=Path);p.add_argument('--build-dir',type=Path)
    p=commands.add_parser('build');p.add_argument('--backend',choices=('cpu','cuda'),default='cpu')
    p.add_argument('--mpi',action='store_true');p.add_argument('--build-dir',type=Path);p.add_argument('--cuda-sm',default='70');p.add_argument('--jobs',type=int,default=1)
    p=commands.add_parser('init');p.add_argument('--problem',default='thin_disk');p.add_argument('--output',type=Path,required=True)
    p=commands.add_parser('inspect');p.add_argument('result',type=Path)
    p=commands.add_parser('validate');p.add_argument('result',type=Path);p.add_argument('--expect-time',type=float)
    p=commands.add_parser('render');p.add_argument('result',type=Path);p.add_argument('--output',type=Path)
    p.add_argument('--xlim',nargs=2,type=float);p.add_argument('--ylim',nargs=2,type=float)
    args=parser.parse_args(argv)
    try:
        if args.command in ('run','plan'):
            spec=RunSpec.load(args.spec)
            changes={k:getattr(args,k) for k in ('backend','mpi','ranks','binary','build_dir') if getattr(args,k) is not None}
            if changes: spec=spec.with_runtime(**changes)
            engine=spec.engine(out_dir=args.out)
            if args.command=='plan': _print(engine.plan().as_dict());return 0
            result=engine.run();_print(dict(out_dir=result.out_dir,status=result.status,complete=result.complete))
            return {'completed':0,'incomplete':3,'stopped':4}.get(result.status,2)
        if args.command=='build':
            artifact=build_native(args.backend,mpi=args.mpi,build_dir=args.build_dir,cuda_sm=args.cuda_sm,jobs=args.jobs)
            _print(artifact);return 0
        if args.command=='init':
            from dataclasses import replace
            spec=RunSpec(Config(problem=args.problem));spec=replace(spec,result=replace(spec.result,out_dir='runs/'+args.problem))
            with args.output.open('x') as f: json.dump(spec.as_dict(),f,indent=2);f.write('\n')
            _print(dict(spec=args.output));return 0
        result=Result.load(args.result)
        if args.command=='inspect': _print(result.manifest)
        elif args.command=='validate': _print(result.validate(expect_time=args.expect_time))
        else: _print(dict(outputs=result.render(args.output,xlim=args.xlim,ylim=args.ylim)))
        return 0
    except SimulationError as error:
        print(f'pyhmd: {error}',file=sys.stderr);return 2
    except (ValueError,TypeError,OSError,RuntimeError,ImportError,subprocess.SubprocessError) as error:
        print(f'pyhmd: {error}',file=sys.stderr);return 2

if __name__=='__main__':raise SystemExit(main())
