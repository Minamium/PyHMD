#!/usr/bin/env python3
"""Validate passive transport CSVs, control-volume closure and MPI duplicates."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

def rows(path):
    with Path(path).open() as handle:return list(csv.DictReader(handle))

def validate(root, full_snapshots=True, expect_time=None, expect_step=None, require_complete=False):
    root=Path(root);meta=json.loads((root/'metadata.json').read_text())
    assert meta['transport_diagnostics']==1
    table=rows(root/'transport_regions.csv');assert table
    status=json.loads((root/'status.json').read_text())
    last_step=max(int(q['step']) for q in table);last_time=max(float(q['time']) for q in table)
    assert status['step']==last_step and abs(status['time']-last_time)<1e-10
    if expect_time is not None:assert abs(last_time-expect_time)<1e-10
    if expect_step is not None:assert last_step==expect_step
    if require_complete:assert status['complete'] and status['reason']=='t_end'
    max_relative={k:0. for k in ('mass','energy','angular')}
    for row in table:
        q={k:float(v) for k,v in row.items() if k!='region'}
        assert all(np.isfinite(v) for v in q.values())
        for key in max_relative:
            scale=1+abs(q[key+'_initial'])+abs(q['radial_'+key+'_out'])+abs(q['polar_'+key+'_out'])
            relative=abs(q[key+'_residual'])/scale
            max_relative[key]=max(max_relative[key],relative)
            assert relative<1e-8,(row['region'],row['time'],key,relative)
    diagnostics={int(float(r['step'])):r for r in rows(root/'diagnostics.csv')}
    boundary_error=0.
    for row in table:
        if row['region']!='domain':continue
        q=diagnostics[int(row['step'])]
        for key,old in (('mass','boundary_mass_out'),('energy','boundary_energy_out'),('angular','boundary_angular_momentum_out')):
            actual=float(row['radial_'+key+'_out'])+float(row['polar_'+key+'_out'])
            error=abs(actual-float(q[old]))/(1+abs(float(q[old])))
            boundary_error=max(boundary_error,error);assert error<1e-9,(key,error)
    mid=rows(root/'transport_midplane.csv');assert mid
    max_psi=max(float(q['psi_change_closure_max']) for q in mid)
    assert max_psi<1e-7,('Psi change closure',max_psi)
    result=dict(passed=True,source_fingerprint=meta['source_fingerprint'],ranks=meta['ranks'],
                integration_complete=status['complete'],time=last_time,step=last_step,
                grid=[meta['nx1'],meta['nx2']],max_relative_control_residual=max_relative,
                max_relative_original_boundary_difference=boundary_error,max_abs_psi_change_residual=max_psi,
                scientific_reproduction_verified=False)
    if full_snapshots:
        directory=root/'transport_snapshots';assert directory.is_dir()
        psi_files=sorted(directory.glob('psi_step*_rank*.csv'))
        face_files=sorted(directory.glob('faces_step*_rank*.csv'));assert psi_files and face_files
        final=max(int(p.name.split('_step')[1].split('_')[0]) for p in psi_files)
        assert final==last_step,'Transport snapshot does not match final accepted state'
        final_psi=sorted(directory.glob(f'psi_step{final:010d}_rank*.csv'))
        final_faces=sorted(directory.glob(f'faces_step{final:010d}_rank*.csv'))
        assert len(final_psi)==len(final_faces)==meta['ranks']
        max_shared=0.;max_split=0.;previous_psi=None;previous_faces=None;max_scale=0.
        for pfile,ffile in zip(final_psi,final_faces):
            pv=np.genfromtxt(pfile,delimiter=',',names=True,skip_header=1)
            ff=np.genfromtxt(ffile,delimiter=',',names=True,skip_header=1)
            for data in (pv,ff):assert all(np.isfinite(data[k]).all() for k in data.dtype.names)
            total=sum(pv[k] for k in ('cumulative_advection','cumulative_uct_jump','cumulative_resistive','cumulative_operator_round'))
            residual=pv['psi']-pv['psi_initial']-total
            assert np.allclose(residual,pv['psi_change_residual'],rtol=0,atol=1e-14)
            max_scale=max(max_scale,float(np.max(np.abs(pv['psi_initial']))))
            assert np.max(np.abs(residual))<1e-9*(1+max_scale)
            for actual,parts in [('mass',('mass_adv','mass_hll_jump')),('angular',('angular_adv','angular_maxwell','angular_hll_jump'))]:
                field=lambda name:ff['cumulative_'+name+'_plus_coordinate']
                error=np.max(np.abs(field(actual)-sum(field(k) for k in parts))/(1+np.abs(field(actual))))
                max_split=max(max_split,float(error));assert error<1e-9,(actual,error)
            left=pv[pv['i']==pv['i'].min()];right=pv[pv['i']==pv['i'].max()]
            left_face=ff[(ff['dir']==0)&(ff['i']==ff['i'].min())]
            right_face=ff[(ff['dir']==0)&(ff['i']==ff['i'].max())]
            if previous_psi is not None:
                assert np.array_equal(left['j'],previous_psi['j'])
                for key in ('psi','psi_initial','cumulative_advection','cumulative_uct_jump','cumulative_resistive','cumulative_operator_round'):
                    error=float(np.max(np.abs(left[key]-previous_psi[key])/(1+np.abs(previous_psi[key]))))
                    max_shared=max(max_shared,error);assert error<1e-9,(key,error)
                assert np.array_equal(left_face['j'],previous_faces['j'])
                for key in left_face.dtype.names[4:]:
                    error=float(np.max(np.abs(left_face[key]-previous_faces[key])/(1+np.abs(previous_faces[key]))))
                    max_shared=max(max_shared,error);assert error<1e-9,(key,error)
            previous_psi,previous_faces=right,right_face
        result.update(final_step=final,max_relative_shared_interface_difference=max_shared,
                      max_relative_flux_decomposition_residual=max_split)
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('run',type=Path)
    parser.add_argument('--output',type=Path);parser.add_argument('--expect-time',type=float)
    parser.add_argument('--expect-step',type=int);parser.add_argument('--require-complete',action='store_true');args=parser.parse_args()
    result=validate(args.run,expect_time=args.expect_time,expect_step=args.expect_step,require_complete=args.require_complete)
    text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')

if __name__=='__main__':main()
