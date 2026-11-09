#!/usr/bin/env python3
"""Measure convergence, including initialization and staggered-face representation errors."""
import argparse
import json
import math
from pathlib import Path
import platform
import sys

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from pyhmd import Config,Engine,RuntimeConfig,ResultOptions
from pyhmd.backends.build import describe,sha256
from pyhmd.results.checkpoint import read_checkpoint
from analytic import references,assemble,errors,frame

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--resolutions',nargs='+',type=int,default=[16,32,64])
    p.add_argument('--amplitudes',nargs='+',type=float,default=[.1]);p.add_argument('--ct-schemes',nargs='+',default=['face_hll'])
    p.add_argument('--linear',action='store_true');p.add_argument('--ranks',type=int,default=1);p.add_argument('--launcher',nargs='+',default=['mpirun'])
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True);info=describe(args.binary)
    runtime=RuntimeConfig(backend=info['backend'],mpi=info['mpi'],ranks=args.ranks,launcher=tuple(args.launcher),binary=args.binary.resolve())
    records=[];passed=True
    jobs=[('oblique_alfven','fast',a,ct) for a in args.amplitudes for ct in args.ct_schemes]
    if args.linear:jobs += [('linear_wave',mode,1e-6,ct) for mode in ('fast','alfven','slow','entropy') for ct in args.ct_schemes]
    for problem,mode,amp,ct in jobs:
        group=[]
        for n in sorted(args.resolutions):
            cfg=Config(problem=problem,grid={'nx1':n,'nx2':n//2},physics={'wave_mode':mode,'wave_amplitude':amp},numerics={'ct_emf_scheme':ct})
            period=1/frame(cfg)[-1]
            from dataclasses import replace
            cfg=replace(cfg,time=replace(cfg.time,t_end=period))
            folder=args.output/f'{problem}_{mode}_a{amp:g}_{ct}_{n}'
            result=Engine(cfg,runtime=runtime,result_options=ResultOptions(out_dir=folder,output_interval=period/4,checkpoint_interval=period)).run()
            result.validate()
            state,psi=assemble(read_checkpoint(result.out_dir/result.checkpoints[-1]['path']))
            volume,discrete,bx,by=references(cfg,period);v0,d0,_,_=references(cfg,0)
            dx=(cfg.grid.x1_max-cfg.grid.x1_min)/n;dy=(cfg.grid.x2_max-cfg.grid.x2_min)/(n//2)
            face_errors={'Bx_L1':float(np.mean(np.abs(np.diff(psi,axis=0)/dy-bx))),
                         'Bx_Linf':float(np.max(np.abs(np.diff(psi,axis=0)/dy-bx))),
                         'By_L1':float(np.mean(np.abs(-np.diff(psi,axis=1)/dx-by))),
                         'By_Linf':float(np.max(np.abs(-np.diff(psi,axis=1)/dx-by)))}
            # A quarter-period check cannot be passed by a wave traveling backwards.
            quarter=result.snapshot(result.snapshots[1]['path'])
            x=quarter.data
            q=np.stack((x['rho'],x['rho']*x['vr'],x['rho']*x['vtheta'],x['rho']*x['vphi'],x['total_energy'],x['Br'],x['Btheta'],x['Bphi'],x['rho']*x['entropy'],x['internal_energy']),axis=-1)
            qref=references(cfg,quarter.time)[1]
            correct=errors(q,qref)['rms_L1_physical']
            backward=errors(q,references(cfg,-quarter.time)[1])['rms_L1_physical']
            stationary=errors(q,d0)['rms_L1_physical']
            direction=dict(correct_error=correct,backward_error=backward,stationary_error=stationary,
                           passed=correct<min(backward,stationary))
            if not direction['passed']:raise AssertionError(f'quarter-period propagation direction failed: {direction}')
            initial=result.snapshot('first').data
            initial_u=np.stack((initial['rho'],initial['rho']*initial['vr'],initial['rho']*initial['vtheta'],initial['rho']*initial['vphi'],initial['total_energy'],initial['Br'],initial['Btheta'],initial['Bphi'],initial['rho']*initial['entropy'],initial['internal_energy']),axis=-1)
            row=dict(problem=problem,mode=mode,amplitude=amp,ct_scheme=ct,resolution=[n,n//2],time=period,
                     cell_volume=errors(state,volume),cell_discrete=errors(state,discrete),initial_representation=errors(d0,v0),
                     initial_measured_volume=errors(initial_u,v0),initial_measured_discrete=errors(initial_u,d0),
                     face=face_errors,quarter_period=errors(q,qref),quarter_direction_check=direction,config=cfg.as_dict(),source=info,
                     linear_reference_approximate=problem=='linear_wave',linear_nonlinear_floor_scale=amp*amp if problem=='linear_wave' else 0.)
            if group:
                prior=group[-1];ratio=n/prior['resolution'][0]
                row['observed_order_L1']=math.log(prior['cell_volume']['rms_L1_physical']/row['cell_volume']['rms_L1_physical'])/math.log(ratio)
                row['observed_order_Linf']=math.log(prior['cell_volume']['max_Linf_physical']/row['cell_volume']['max_Linf_physical'])/math.log(ratio)
            group.append(row);records.append(row)
            print(json.dumps({k:row[k] for k in ('problem','mode','amplitude','resolution','observed_order_L1','observed_order_Linf') if k in row}),flush=True)
    report=dict(schema_version=1,measurement_completed=passed,scientific_reproduction_verified=False,binary_sha256=sha256(args.binary),
                environment={'platform':platform.platform(),'python':platform.python_version(),'numpy':np.__version__,'ranks':args.ranks},
                reference='Independent Gaussian volume averages and analytic sinc face averages; Felker & Stone 2018 definitions',
                limitation='Linear waves include O(amplitude^2) nonlinear error. Current CT cell-B representation is second-order at initialization. No high-order claim from scheme names.',results=records)
    (args.output/'orders.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return 0

if __name__=='__main__':raise SystemExit(main())
