#!/usr/bin/env python3
"""Compare a matched passive-transport resolution pair without changing physics."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from pyhmd.diagnostics.transport import validate

def table(path):
    with Path(path).open() as handle:return list(csv.DictReader(handle))
def numbers(row):return {k:float(v) for k,v in row.items() if k!='region'}

def analyze(root):
    meta=json.loads((root/'metadata.json').read_text());quality=validate(root,expect_time=1000,require_complete=True)
    mid=[numbers(q) for q in table(root/'transport_midplane.csv')]
    regions=table(root/'transport_regions.csv');disk=[numbers(q) for q in regions if q['region']=='disk_r10']
    mu_region=[numbers(q) for q in regions if q['region']=='mu_disk'][-1]
    theta_cut=disk[0]['theta_lower'];boundary=table(root/'transport_inner_boundary.csv')
    summed={}
    for row in boundary:
        q=numbers(row)
        if q['theta']<theta_cut:continue
        step=int(q['step']);item=summed.setdefault(step,dict(time=q['time'],mass=0.,adv=0.,jump=0.))
        for key,field in [('mass','mass'),('adv','mass_adv'),('jump','mass_hll_jump')]:
            item[key]+=q['cumulative_'+field+'_out']
    inner=list(summed.values());final_step=int(mid[-1]['step']);equator=[]
    for path in sorted((root/'transport_snapshots').glob(f'psi_step{final_step:010d}_rank*.csv')):
        a=np.genfromtxt(path,delimiter=',',names=True,skip_header=1)
        equator.extend(a[(a['j']==meta['nx2'])&(a['owner']==1)])
    equator=np.array(equator,dtype=a.dtype);equator.sort(order='i')
    rf=np.asarray(meta['x1_faces']);ilo=int(np.argmin(abs(rf-mu_region['r_inner'])));ihi=int(np.argmin(abs(rf-mu_region['r_outer'])))
    fields=['psi','psi_initial','cumulative_advection','cumulative_uct_jump','cumulative_resistive','cumulative_operator_round','psi_change_residual']
    at=lambda i:equator[equator['i']==i][0]
    flux={key:float(at(ihi)[key]-at(ilo)[key]) for key in fields}
    flux.update(r_inner=float(rf[ilo]),r_outer=float(rf[ihi]),definition='Psi(r_outer,equator)-Psi(r_inner,equator), per radian')
    return dict(metadata=meta,quality=quality,midplane=mid,disk=disk,inner_boundary=inner,annular_poloidal_flux=flux)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('coarse',type=Path);parser.add_argument('fine',type=Path)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--legacy-coarse',type=Path);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    data=[analyze(args.coarse),analyze(args.fine)]
    settings=('gamma','gm','epsilon','mu0','alpha_m','chi','eta_model','cool_ohmic','cfl','initial_equilibrium','ct_emf_scheme','inner_normal_field','internal_energy_sync','dual_energy_eta','metric_ppm')
    assert all(data[0]['metadata'][k]==data[1]['metadata'][k] for k in settings)
    for i in (0,-1):
        for key in ('r_inner','r_outer','theta_lower','theta_upper'):
            assert np.isclose(data[0]['disk'][i][key],data[1]['disk'][i][key],rtol=1e-13,atol=1e-14),(key,i)
    summary=dict(cases=[],same_physics=True,matched_control_faces=True,target_time=1000,
                 scientific_reproduction_verified=False,causal_identification_complete=False,
                 limitation='Two resolutions at t=1000 test early resolution dependence, not asymptotic convergence or full paper reproduction')
    for root,case in zip((args.coarse,args.fine),data):
        summary['cases'].append(dict(run=str(root),grid=[case['metadata']['nx1'],case['metadata']['nx2']],quality=case['quality'],
            initial_mu=case['midplane'][0]['mu_r1p1_r1p5'],final_mu=case['midplane'][-1]['mu_r1p1_r1p5'],
            initial_disk_mass=case['disk'][0]['mass'],final_disk_mass=case['disk'][-1]['mass'],
            final_inner_boundary=case['inner_boundary'][-1],annular_poloidal_flux=case['annular_poloidal_flux']))
    summary['fine_over_coarse_final_mu']=summary['cases'][1]['final_mu']/summary['cases'][0]['final_mu']
    if args.legacy_coarse:
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
        from pyhmd.plotting._viz import load_snapshot
        legacy={}
        for path in args.legacy_coarse.glob('snapshot_*.csv'):
            line=path.open().readline();time=float(line.split('time=')[1].split()[0]);legacy[time]=path
        comparisons=[]
        for path in sorted(args.coarse.glob('snapshot_*.csv')):
            current=load_snapshot(args.coarse,path)
            nearest=min(legacy,key=lambda time:abs(time-current.time)) if legacy else None
            old_path=legacy[nearest] if nearest is not None and abs(nearest-current.time)<1e-10*max(1,current.time) else None
            if old_path is None:continue
            old=load_snapshot(args.legacy_coarse,old_path)
            fields={k:dict(max_abs=float(np.max(np.abs(current.data[k]-old.data[k]))),
                           passed=bool(np.allclose(current.data[k],old.data[k],rtol=1e-9,atol=1e-12))) for k in old.data}
            comparisons.append(dict(time=current.time,passed=all(q['passed'] for q in fields.values()),fields=fields))
        summary['legacy_coarse_comparison']=comparisons
        summary['legacy_coarse_match_verified']=len(comparisons)>=11 and all(q['passed'] for q in comparisons)
    (args.output/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11.5,8),layout='constrained')
    for case,label,color,ax in zip(data,('600 x 128','1200 x 256'),('#2474b7','#df7127'),axes[1]):
        mid=case['midplane'];disk=case['disk'];inner=case['inner_boundary']
        axes[0,0].semilogy([q['time'] for q in mid],[q['mu_r1p1_r1p5'] for q in mid],label=label,color=color)
        axes[0,1].plot([q['time'] for q in disk],[q['mass']/disk[0]['mass'] for q in disk],label=label,color=color)
        for key,ls in [('mass','-'),('adv','--'),('jump',':')]:
            ax.plot([q['time'] for q in inner],[q[key]*4*np.pi for q in inner],ls=ls,label=key,lw=1.4)
        ax.axhline(0,color='gray',lw=.8);ax.set(title=label+' inner disk boundary',ylabel='Cumulative mass OUT [bipolar]')
        ax.legend(fontsize=8)
    axes[0,0].set(title='Authoritative equatorial-face magnetization',ylabel='mu, r=1.1..1.5')
    axes[0,1].set(title='Mass in matched disk control volume',ylabel='M / M(initial)')
    axes[0,0].legend();axes[0,1].legend()
    for ax in axes.flat:ax.set_xlabel('Normalized time');ax.grid(alpha=.22)
    fig.suptitle('Same physics, different resolution | actual accepted FV/CT transport')
    fig.supxlabel('Inner boundary: negative OUT means injection. HLL advective/jump parts are numerical-flux components.',fontsize=9)
    fig.savefig(args.output/'transport_comparison.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4.7),layout='constrained')
    labels=['advection','UCT jump','physical resistance','operator round']
    for offset,case,label in zip((-.18,.18),data,('600 x 128','1200 x 256')):
        q=case['annular_poloidal_flux'];values=[q[k]/q['psi_initial'] for k in ('cumulative_advection','cumulative_uct_jump','cumulative_resistive','cumulative_operator_round')]
        ax.bar(np.arange(4)+offset,values,width=.36,label=label)
    ax.axhline(0,color='gray',lw=.8);ax.set_xticks(range(4),labels);ax.set(ylabel='Change / initial annular Psi difference',title='Equatorial magnetic flux transport over matched radial interval')
    ax.legend();ax.grid(axis='y',alpha=.22)
    fig.supxlabel('The UCT jump term does not represent all numerical error. Actual Psi-change closure is checked separately.',fontsize=8.5)
    fig.savefig(args.output/'magnetic_transport.png',dpi=170);plt.close(fig)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
