"""Conditional perturbation sensitivity, not an unconditional numerical error bound."""
import csv,json,os,signal,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from pyhmd import Config,RunSpec,RuntimeConfig,ResultOptions,Result
from pyhmd.backends.build import describe,sha256
from pyhmd.engine import _stop_owned_process

def bounded_run(command,env,log,seconds=1800):
    process=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=os.name=='posix')
    try:return process.wait(timeout=seconds)
    except subprocess.TimeoutExpired:
        try:
            if os.name=='posix':os.killpg(process.pid,signal.SIGINT)
            else:process.send_signal(signal.SIGINT)
            process.wait(timeout=45)
        except (ProcessLookupError,subprocess.TimeoutExpired):_stop_owned_process(process)
        return 124

def table(path):
    with Path(path).open() as file:return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(file)]

def run_ensemble(binary,output,profile,ct,members,end,grid,launcher):
    output=Path(output);output.mkdir(parents=True,exist_ok=False);info=describe(binary)
    report=dict(schema_version=1,status='running',profile=profile,ct_scheme=ct,grid=list(grid),t_end=end,
                density_perturbation=1e-13,source=info,binary_sha256=sha256(binary),members=[],
                scientific_reproduction_verified=False,noise_floor_qualified=False,
                method='Uniform [-1,1] deterministic global-index density perturbation at fixed initial pressure; conditional on grid/scheme/FP profile',
                pressure_convention='raw total-energy pressure',band=[1.1,1.5],normalization='one simulated hemisphere, per radian azimuth')
    def save(): (output/'ensemble.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    curves=[];baseline=None
    for index in range(members+1):
        amplitude=0 if index==0 else 1e-13
        config=Config(problem='thin_disk',grid={'nx1':grid[0],'nx2':grid[1]},physics={'noise_amplitude':amplitude,'noise_seed':index},
                      numerics={'ct_emf_scheme':ct},time={'t_end':end})
        folder=output/f'member_{index:02}'
        runtime=RuntimeConfig(backend=info['backend'],mpi=info['mpi'],ranks=4 if info['mpi'] else 1,launcher=(launcher,),binary=binary)
        spec=RunSpec(config,runtime,ResultOptions(out_dir=folder,output_interval=end,log_interval=1000,checkpoint_interval=end,midplane_diagnostics=True))
        path=output/f'member_{index:02}.json';path.write_text(json.dumps(spec.as_dict(),indent=2)+'\n')
        save();env=dict(os.environ);env['PYTHONPATH']=str(ROOT/'src')+os.pathsep+env.get('PYTHONPATH','')
        with (output/f'member_{index:02}.log').open('w') as log:
            code=bounded_run([sys.executable,'-m','pyhmd','run',str(path)],env,log)
        item=dict(seed=index,amplitude=amplitude,returncode=code,valid=False)
        try:
            result=Result.load(folder);item['status']=result.status;item['integrity']=result.validate()
            if not result.complete:raise ValueError('integration incomplete')
            midplane=table(folder/'midplane_metrics.csv');diagnostics=table(folder/'diagnostics.csv')
            by_step={int(q['step']):q for q in diagnostics}
            curve=[]
            for sample in midplane:
                diag=by_step[int(sample['step'])]
                curve.append(dict(time=sample['time'],mu=sample['mu_r1p1_r1p5'],delta_psi=sample['delta_psi_r1p1_r1p5'],
                                  mass=diag['mass'],boundary_mass_out=diag['boundary_mass_out']))
            item['final']=curve[-1];item['valid']=True
            item['native_sample_count']=len(curve)
            (folder/'sampled_metrics.json').write_text(json.dumps(curve,indent=2)+'\n')
            if index==0:baseline=curve
            else:curves.append(curve)
        except Exception as error:item['error']=f'{type(error).__name__}: {error}'
        report['members'].append(item);save()
        if code in (124,137):raise RuntimeError('ensemble member timed out; stop before another GPU case')
        if index==0 and not item['valid']:
            report.update(status='baseline_invalid',noise_floor_qualified=False);save();return report
    if len(curves)==members:
        times=np.array(sorted({0.,end,*[t for t in (1.,2.,5.,10.,20.,50.,100.) if t<end]}))
        quantiles={};final_quantiles={}
        for key in ('mu','delta_psi','mass','boundary_mass_out'):
            data=np.array([np.interp(times,[x['time'] for x in curve],[x[key] for x in curve]) for curve in curves])
            quantiles[key]={str(q):np.quantile(data,q,axis=0).tolist() for q in (0.,.05,.5,.95,1.)}
            final=np.array([curve[-1][key] for curve in curves])
            final_quantiles[key]=dict(mean=float(np.mean(final)),std=float(np.std(final,ddof=1)),
                quantiles={str(q):float(np.quantile(final,q)) for q in (0.,.05,.5,.95,1.)},unperturbed=baseline[-1][key])
        report.update(status='completed',noise_floor_qualified=True,final_time_quantiles_measured=final_quantiles,
                      time_grid=times.tolist(),time_quantiles_interpolated=quantiles,
                      interpolation_note='Only final t_end values are sampled at identical physical time. Interior time-grid quantiles use linear interpolation of accepted-step samples; they are not direct measurements.')
    else:report.update(status='members_invalid',noise_floor_qualified=False)
    save();return report
