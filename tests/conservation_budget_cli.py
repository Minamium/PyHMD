#!/usr/bin/env python3
"""Independently check CLI snapshot integrals and reported open/periodic budgets.

The C++ budget test exercises operators. This test additionally checks the
actual CSV columns, host output path and restarted segment accounting.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def table(path):
    with path.open() as handle:
        rows=list(csv.DictReader(line for line in handle if not line.startswith('#')))
    require(bool(rows), f'{path}: empty table')
    require(all(None not in r and all(v is not None for v in r.values()) for r in rows),
            f'{path}: inconsistent CSV column count')
    values=[{k:float(v) for k,v in r.items()} for r in rows]
    require(all(math.isfinite(v) for r in values for v in r.values()), f'{path}: nonfinite output')
    return values


def integrals(path, meta):
    rows=table(path)
    nx,ny=meta['nx1'],meta['nx2']
    require(len(rows)==nx*ny, f'{path}: incomplete grid')
    xf,yf=meta['x1_faces'],meta['x2_faces']
    spherical=meta['geometry']=='spherical'
    mass,energy,angular=[],[],[]
    for k,row in enumerate(rows):
        j,i=divmod(k,nx)
        if spherical:
            volume=(xf[i+1]**3-xf[i]**3)/3*(math.cos(yf[j])-math.cos(yf[j+1]))
            potential=-meta['gm']/row['x1']
            radius=row['x1']*math.sin(row['x2'])
        else:
            volume=(xf[i+1]-xf[i])*(yf[j+1]-yf[j]);potential=0;radius=1
        mass.append(row['rho']*volume)
        energy.append((row['total_energy']+row['rho']*potential)*volume)
        angular.append(radius*row['rho']*row['vphi']*volume)
        require(row['eta_p']>=0 and row['mu_poloidal']>=0, 'negative transport diagnostic')
    return tuple(math.fsum(v) for v in (mass,energy,angular))


def check(root, periodic=False):
    meta=json.loads((root/'metadata.json').read_text())
    snapshots=sorted(root.glob('snapshot_*.csv'))
    require(len(snapshots)>=2, 'no advanced snapshot')
    first,last=integrals(snapshots[0],meta),integrals(snapshots[-1],meta)
    diagnostics=table(root/'diagnostics.csv');a,b=diagnostics[0],diagnostics[-1]
    delta=lambda name:b[name]-a[name]
    corrections=sum(delta(name) for name in ('floor_energy','ct_boundary_energy','entropy_energy','ct_energy','internal_energy'))
    mass=last[0]-first[0]+delta('boundary_mass_out')-delta('floor_mass')
    energy=last[1]-first[1]+delta('boundary_energy_out')+delta('compatible_cooling')-corrections-delta('floor_potential_energy')
    angular=last[2]-first[2]+delta('boundary_angular_momentum_out')
    for name,error,scale in zip(('mass','energy','angular_momentum'),(mass,energy,angular),first):
        tolerance=2e-11*max(1,abs(scale))
        require(abs(error)<tolerance, f'{root.name}: independent {name} residual {error}')
        require(abs(b[name+'_balance_error']-error)<tolerance,
                f'{root.name}: reported {name} balance differs from snapshot integral')
    require(b['floor_events']==0, f'{root.name}: unexpected floor fixture')
    if periodic:
        for field in ('boundary_mass_out','boundary_energy_out','boundary_angular_momentum_out'):
            require(abs(delta(field))<1e-12, f'periodic net boundary transport: {field}')
        require(delta('compatible_cooling')>1e-8, 'periodic fixture has no Ohmic cooling')
    else:
        require(abs(delta('boundary_mass_out'))>1e-9, 'open fixture has no mass transport')
        require(abs(delta('boundary_angular_momentum_out'))>1e-12, 'open fixture has no angular-momentum transport')
        require(abs(delta('ct_boundary_energy'))+abs(delta('ct_energy'))+abs(delta('internal_energy'))>1e-12,
                'open fixture does not exercise energy corrections')
    return {'mass_residual':mass,'energy_residual':energy,'angular_residual':angular,
            'mass_out':delta('boundary_mass_out'),'energy_out':delta('boundary_energy_out'),
            'angular_out':delta('boundary_angular_momentum_out'),'cooling':delta('compatible_cooling'),
            'ct_and_thermal_corrections':corrections,'first':a,'last':b}


def run(binary, base, name, config, extra=()):
    cfg=base/(name+'.ini');cfg.write_text(config)
    root=base/name
    proc=subprocess.run([str(binary),'--config',str(cfg),'--output',str(root),*extra],
                        capture_output=True,text=True)
    (base/(name+'.log')).write_text(proc.stdout+proc.stderr)
    require(proc.returncode==0,f'{name} exited {proc.returncode}: {proc.stderr[-1000:]}')
    return root


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--output',type=Path);args=parser.parse_args()
    binary=args.binary.resolve()
    base=args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='pymhd_budget_cli_'))
    base.mkdir(parents=True,exist_ok=True)
    disk='''case=thin_disk
nx1=24
nx2=16
x1_min=1
x1_max=4
t_end=0.03
output_interval=0.02
log_interval=1
cfl=0.2
thermal_recovery=internal_energy
entropy_ratio=0.9
ct_energy_correction=yes
checkpoint_interval=0.03
'''
    periodic='''case=alfven
nx1=32
nx2=8
eta=0.02
cooling=yes
thermal_recovery=internal_energy
entropy_ratio=0.9
ct_energy_correction=yes
t_end=0.03
output_interval=0.01
log_interval=1
'''
    open_root=run(binary,base,'open_spherical',disk)
    result={'open_spherical':check(open_root)}
    periodic_root=run(binary,base,'periodic_resistive',periodic)
    result['periodic_resistive']=check(periodic_root,periodic=True)
    checkpoints=sorted(open_root.glob('checkpoint_step*'))
    require(len(checkpoints)==1,'expected final open-flow checkpoint')
    resumed_root=run(binary,base,'resumed_spherical',disk,
                     ('--restart',str(checkpoints[0]),'--t-end','0.05'))
    result['resumed_spherical']=check(resumed_root)
    for field in ('boundary_mass_out','boundary_energy_out','boundary_angular_momentum_out',
                  'floor_potential_energy','ct_energy','ct_boundary_energy','internal_energy','compatible_cooling'):
        require(result['resumed_spherical']['first'][field]==result['open_spherical']['last'][field],
                'restart changed cumulative ledger '+field)
    # t=0.03 was a final output between the t=.02 and .04 scheduled outputs.
    comments=[p.open().readline() for p in sorted(resumed_root.glob('snapshot_*.csv'))]
    require(any(line.startswith('# time=0.040000000000000001 ') or line.startswith('# time=0.04 ') for line in comments),
            'restart skipped the next scheduled snapshot at t=.04')
    (base/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS CLI open/periodic/restarted conservation budgets:',base)
    for name,r in result.items():
        print(name,{k:r[k] for k in ('mass_residual','energy_residual','angular_residual','mass_out','cooling')})


if __name__=='__main__':
    main()
