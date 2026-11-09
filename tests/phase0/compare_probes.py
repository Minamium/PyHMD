#!/usr/bin/env python3
"""Compare common-input probes using an input-stencil scale, not a small residual."""
import argparse,json
from pathlib import Path
import numpy as np

def load(root,name):
    meta=[json.loads(p.read_text()) for p in sorted(root.glob('rank*.json'))]
    if not meta:raise ValueError(f'No probe shards: {root}')
    ny,nx=meta[0]['ny'],meta[0]['global_nx'];vertex=name.endswith('Psi') or name=='emf'
    result=np.empty((ny+int(vertex),nx+int(vertex))) if vertex else np.empty((ny,nx,10))
    for item in meta:
        shape=(ny+1,item['nx']+1) if vertex else (ny,item['nx'],10)
        a=np.fromfile(root/f'{name}_rank{item["rank"]:06}.bin',dtype='<f8').reshape(shape)
        result[:,item['start']:item['start']+item['nx']+int(vertex)]=a
    if not np.all(np.isfinite(result)):raise ValueError('Nonfinite probe output')
    return result,meta

def compare(left,right,tolerance=1e-12):
    initial,ma=load(left,'input_U');other,mb=load(right,'input_U')
    psi,_=load(left,'input_Psi');other_psi,_=load(right,'input_Psi')
    same=initial.tobytes()==other.tobytes() and psi.tobytes()==other_psi.tobytes()
    scale=max(float(np.max(np.abs(initial))),1e-30)
    rate=max(max(q['rate'] for q in ma),max(q['rate'] for q in mb),1.)
    def settings(root):
        lines=(root/'parameters.signature.txt').read_text().splitlines()[2:]
        if lines[0]=='fp_contract_off':lines=lines[1:]
        return [[float(x) for x in line.split()] for line in lines]
    pa,pb=settings(left),settings(right)
    parameter_match=all(a[:min(len(a),len(b))]==b[:min(len(a),len(b))] and
                        all(x==0 for x in a[min(len(a),len(b)):]+b[min(len(a),len(b)):]) for a,b in zip(pa,pb)) and len(pa)==len(pb)
    result=dict(same_initial_bits=same,physics_settings_match=parameter_match,tolerance=tolerance,stencil_scale=scale*rate,fields={})
    for name in ('rhs','emf','final_U','final_Psi'):
        a,_=load(left,name);b,_=load(right,name);error=float(np.max(np.abs(a-b)))
        denominator=scale*rate if name in ('rhs','emf') else max(scale,float(np.max(np.abs(psi))),1e-30)
        result['fields'][name]=dict(max_abs=error,scaled=error/denominator,bits_equal=a.tobytes()==b.tobytes())
    result['same_final_clock']=ma[0]['time']==mb[0]['time']
    result['final_clock_abs_difference']=abs(ma[0]['time']-mb[0]['time'])
    result['integration_timing']={label:dict(ranks=len(metadata),
        max_rank_wall_seconds=max(q['integration_wall_seconds'] for q in metadata),
        accepted_steps=metadata[0]['steps'],physical_time=metadata[0]['time'])
        for label,metadata in [('left',ma),('right',mb)]}
    result['passed_rhs']=same and parameter_match and all(result['fields'][key]['scaled']<=tolerance for key in ('rhs','emf'))
    # The requested long-step/noise comparison is reported, not given a fake RHS bound.
    result['trajectory_within_initial_rhs_tolerance']=all(result['fields'][key]['scaled']<=tolerance for key in ('final_U','final_Psi'))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('left',type=Path);p.add_argument('right',type=Path);p.add_argument('--report',type=Path,required=True)
    args=p.parse_args();report=compare(args.left,args.right);args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));raise SystemExit(0 if report['passed_rhs'] else 1)
