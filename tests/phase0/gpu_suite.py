#!/usr/bin/env python3
"""Allocated single-node Phase 0: frozen controls, RHS/MPI, convergence, noise.

No scheduler submission, SSH, global environment modification, or signature
rewriting. Site setup and allocation validation are the batch wrapper's job.
"""
import argparse,csv,hashlib,json,os,platform,shlex,subprocess,sys,time
from dataclasses import replace
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tests/convergence')]
import numpy as np
from pymhd import Config,Engine,RuntimeConfig,ResultOptions
from pymhd.config import native_ini
from pymhd.results.checkpoint import read_checkpoint,compare_checkpoints
from pymhd.backends.build import describe
from compare_probes import compare

def native_layout(source):
    for package in ('pymhd', 'pyhmd'):
        directory = source / 'src' / package / 'native'
        if (directory / 'fingerprint.py').is_file():
            return directory, package
    raise FileNotFoundError(f'Bundled native sources not found under {source}')

class Suite:
    def __init__(self,args):
        self.a=args;self.out=args.output.resolve();self.out.mkdir(parents=True,exist_ok=False)
        frozen=args.reference_source.resolve().parent
        proof=json.loads((frozen/'freeze.json').read_text())
        if hashlib.sha256((frozen/'source.tar').read_bytes()).hexdigest()!=proof['archive_sha256']:
            raise ValueError('frozen source archive checksum differs')
        native, _ = native_layout(args.reference_source.resolve())
        fingerprint=subprocess.check_output([sys.executable,str(native/'fingerprint.py')],text=True).strip()
        if fingerprint!=proof['description']['source_fingerprint']:raise ValueError('reference source differs from frozen binary generation')
        for name,expected in proof['configs'].items():
            if hashlib.sha256((args.reference_source.resolve()/'configs'/name).read_bytes()).hexdigest()!=expected:
                raise ValueError(f'frozen configuration differs: {name}')
        for name in ('logs','bin','configs','reports','cases'): (self.out/name).mkdir(exist_ok=True)
        self.report=dict(schema_version=1,status='running',phase0_gate_complete=False,checks=[],
            environment=dict(platform=platform.platform(),python=platform.python_version(),numpy=np.__version__,
                             pbs_job_id=os.environ.get('PBS_JOBID'),cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES')))
        self.binaries={}
    def save(self): (self.out/'progress.json').write_text(json.dumps(self.report,indent=2)+'\n')
    def run(self,command,label,accepted=(0,)):
        self.report['stage']=label;self.save();start=time.monotonic()
        with (self.out/'logs'/f'{label}.log').open('w') as log:
            log.write(json.dumps(list(map(str,command)))+'\n');log.flush()
            process=subprocess.run(list(map(str,command)),stdout=log,stderr=subprocess.STDOUT,timeout=1800)
        record=dict(check=label,returncode=process.returncode,seconds=time.monotonic()-start,accepted=process.returncode in accepted)
        self.report['checks'].append(record);self.save()
        if process.returncode not in accepted:raise RuntimeError(f'{label} failed; inspect logs/{label}.log')
        return process.returncode
    def mpi(self,binary,ranks,args=()):
        return [os.environ.get('MPIEXEC','mpirun'),'--bind-to','core','-np',str(ranks),str(binary),*map(str,args)]
    def config(self,name,problem,nx,ny,**overrides):
        cfg=Config(problem=problem,grid={'nx1':nx,'nx2':ny})
        values=cfg.native_settings();values.update(output_interval=1000000,log_interval=1000,checkpoint_interval=1000000)
        values.update(overrides)
        # Old controls do not implement new optional fields. Values are their
        # disabled defaults; common INI bytes must be identical for both sources.
        for key in ('wave_mode','noise_seed','noise_amplitude'):values.pop(key,None)
        path=self.out/'configs'/f'{name}.ini';path.write_text(native_ini(values));return path
    def build(self):
        cxx=os.environ.get('CXX','c++');mpi=os.environ.get('MPICXX','mpicxx');nvcc=os.environ.get('NVCC','nvcc')
        mpi_compile=shlex.split(subprocess.check_output([mpi,'--showme:compile'],text=True))
        mpi_link=[x for x in shlex.split(subprocess.check_output([mpi,'--showme:link'],text=True)) if x.startswith(('-L','-l'))]
        for tool in (cxx,mpi,nvcc):
            self.run([tool,'--version'],'version_'+Path(tool).name)
        for generation,source in [('frozen_control',self.a.reference_source.resolve()),('candidate',ROOT)]:
            native, binary_prefix = native_layout(source)
            fingerprint=subprocess.check_output([sys.executable,str(native/'fingerprint.py')],text=True).strip()
            for mode in ('default','strict'):
                directory=self.out/'bin'/generation/mode;directory.mkdir(parents=True)
                strict=mode=='strict';cpu_flags='-O3 -std=c++17 -Wall -Wextra -Wpedantic'+(' -ffp-contract=off' if strict else '')
                cuda_flags=f'-O3 -lineinfo -Xptxas=-v -ccbin {cxx}'+(' -fmad=false' if strict else '')
                command=['make','-B','-C',native,'cpu-mpi','cuda',f'BUILD_DIR={directory}',f'CXX={cxx}',f'MPICXX={mpi}',f'NVCC={nvcc}',
                         'CUDA_SM=70',f'CXXFLAGS={cpu_flags}',f'CUDAFLAGS={cuda_flags}']
                if generation=='candidate':command += [f'FP_MODE={mode}']
                self.run(command,'build_'+generation+'_'+mode)
                suffix='_strict' if strict and generation=='candidate' else ''
                entry=dict(cpu=directory/(binary_prefix+'_cpu_mpi'+suffix),cuda=directory/(binary_prefix+'_cuda'+suffix))
                for backend,compiler in [('cpu',mpi),('cuda',nvcc)]:
                    for probe in ('operator_probe','periodic_probe','retry_checkpoint'):
                        output=directory/(probe+'_'+backend)
                        common=['-O3','-std=c++17','-DPYMHD_MPI',f'-DPYMHD_SOURCE_FINGERPRINT="{fingerprint}"','-I'+str(native/'include')]
                        if strict:common += ['-DPYMHD_STRICT_FP']
                        if backend=='cuda':
                            command=[compiler,'-x','cu','--extended-lambda','-arch=sm_70','-DPYMHD_CUDA','-ccbin',cxx,*common,*mpi_compile]
                            if strict:command += ['-fmad=false']
                        else:command=[compiler,*common,*(['-ffp-contract=off'] if strict else [])]
                        command += [ROOT/'tests/phase0'/f'{probe}.cpp','-o',output]
                        if backend=='cuda':command += mpi_link
                        self.run(command,f'compile_{generation}_{mode}_{probe}_{backend}')
                        entry[probe+'_'+backend]=output
                self.binaries[generation,mode]=entry
                if generation=='candidate' and mode=='default':
                    output=directory/'cached_rhs_cuda'
                    self.run([nvcc,'-x','cu','-O3','-std=c++17','--extended-lambda','-arch=sm_70','-DPYMHD_CUDA','-DPYMHD_MPI','-ccbin',cxx,
                              '-I'+str(native/'include'),*mpi_compile,ROOT/'tests/gpu_rhs_equivalence.cpp','-o',output,*mpi_link],
                             'compile_cached_rhs_cuda')
                    entry['cached_rhs_cuda']=output
        files=list((self.out/'bin').rglob('*'))
        hashes={str(p.relative_to(self.out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.is_file()}
        (self.out/'reports/frozen_binary_manifest.json').write_text(json.dumps(dict(binaries=hashes,reference_origin='new platform controls from checksum-verified archived pre-change source, not historical SM70 binaries'),indent=2)+'\n')
    def verify(self):
        self.run(self.mpi(self.binaries['candidate','default']['cached_rhs_cuda'],4,[1200,256]),'cached_uncached_cuda_rhs')
        for ranks in (1,2,4):
            self.run(self.mpi(self.binaries['frozen_control','default']['periodic_probe_cpu'],ranks),'periodic_baseline_'+str(ranks),accepted=(0,) if ranks==1 else (1,))
            for backend in ('cpu','cuda'):
                self.run(self.mpi(self.binaries['candidate','default']['periodic_probe_'+backend],ranks),'periodic_candidate_'+backend+'_'+str(ranks))
        for mode in ('default','strict'):
            ref=self.binaries['frozen_control',mode];candidate=self.binaries['candidate',mode]
            for ranks in (1,2,4):
                self.run([sys.executable,ROOT/'scripts/verify_migration.py','--reference',ref['cpu'],'--candidate',candidate['cpu'],
                          '--reference-configs',self.a.reference_source/'configs','--retry-reference',ref['retry_checkpoint_cpu'],
                          '--retry-candidate',candidate['retry_checkpoint_cpu'],'--ranks',ranks,'--report',self.out/'reports'/f'migration_cpu_{mode}_{ranks}.json',
                          '--launcher',os.environ.get('MPIEXEC','mpirun')],f'migration_cpu_{mode}_{ranks}')
            for name,problem,nx,ny,ct in [('disk','thin_disk',600,128,'face_hll'),('corner','thin_disk',64,32,'corner_hll'),('oblique','oblique_alfven',64,32,'face_hll')]:
                cfg=self.config(name+'_'+mode,problem,nx,ny,ct_emf_scheme=ct)
                if problem=='oblique_alfven':
                    # Legacy controls cannot initialize a newly added benchmark.
                    control=candidate
                else:control=ref
                initial=self.out/'cases'/f'{name}_{mode}_control_cpu_1'
                self.run(self.mpi(control['operator_probe_cpu'],1,[cfg,initial,'-',100]),f'probe_{name}_{mode}_control_cpu_1')
                for generation,entry in [('candidate',candidate)]:
                    for ranks in (1,2,4):
                        cpu=self.out/'cases'/f'{name}_{mode}_cpu_{ranks}'
                        gpu=self.out/'cases'/f'{name}_{mode}_cuda_{ranks}'
                        for backend,path in [('cpu',cpu),('cuda',gpu)]:
                            self.run(self.mpi(entry['operator_probe_'+backend],ranks,[cfg,path,initial,100]),f'probe_{name}_{mode}_{backend}_{ranks}')
                        for label,left,right in [('cpu_vs_cuda',cpu,gpu),('partition',initial,cpu),('generation',initial,gpu)]:
                            data=compare(left,right);data.update(label=label,mode=mode,problem=problem,ranks=ranks)
                            (self.out/'reports'/f'{name}_{mode}_{ranks}_{label}.json').write_text(json.dumps(data,indent=2)+'\n')
                            if not data['passed_rhs']:raise RuntimeError(f'common-input RHS gate failed: {name}/{mode}/{ranks}/{label}')
        # Actual GPU checkpoint resume and rejected-step coverage, same generation.
        candidate=self.binaries['candidate','default']
        self.run([sys.executable,ROOT/'scripts/verify_migration.py','--reference',candidate['cuda'],'--candidate',candidate['cuda'],
                  '--reference-configs',self.a.reference_source/'configs','--target-grid','--retry-reference',candidate['retry_checkpoint_cuda'],
                  '--retry-candidate',candidate['retry_checkpoint_cuda'],'--ranks',4,'--report',self.out/'reports/gpu_repeat_restart_retry.json',
                  '--launcher',os.environ.get('MPIEXEC','mpirun')],'gpu_restart_retry')
    def convergence(self):
        binary=self.binaries['candidate','default']['cuda']
        self.run([sys.executable,ROOT/'tests/convergence/run.py','--binary',binary,'--output',self.out/'convergence',
                  '--resolutions','16','32','64','128','256','512','--amplitudes','.2','.1','.01','--ct-schemes','face_hll','corner_hll',
                  '--linear','--ranks','4','--launcher',os.environ.get('MPIEXEC','mpirun')],'convergence')
    def noise(self):
        from noise_ensemble import run_ensemble
        for mode in ('default','strict'):
            for ct in ('face_hll','corner_hll'):
                self.report['stage']='noise_'+mode+'_'+ct;self.save()
                outcome=run_ensemble(self.binaries['candidate',mode]['cuda'],self.out/f'noise_{mode}_{ct}',
                             mode,ct,self.a.ensemble_members,self.a.ensemble_time,self.a.ensemble_grid,
                             os.environ.get('MPIEXEC','mpirun'))
                self.report['checks'].append(dict(check=self.report['stage'],accepted=outcome['noise_floor_qualified']));self.save()
    def execute(self):
        try:
            self.build();self.verify();self.convergence();self.noise()
            self.report.update(status='completed' if all(c['accepted'] for c in self.report['checks']) else 'completed_with_unqualified_groups',phase0_gate_complete=False,
                pending=['Independent Brio-Wu/RJ2a accuracy references','WP0.5 component elision/communication deletion deliberately deferred until measured',
                         'Human assessment of noise floor and cross-generation fixed-step trajectory differences'])
        except BaseException as error:
            self.report.update(status='failed',error=f'{type(error).__name__}: {error}');raise
        finally:self.save()

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--reference-source',type=Path,required=True);parser.add_argument('--ensemble-members',type=int,default=16)
    parser.add_argument('--ensemble-time',type=float,default=200);parser.add_argument('--ensemble-grid',nargs=2,type=int,default=[600,128]);parser.add_argument('--plan',action='store_true')
    args=parser.parse_args()
    if args.ensemble_members<16:parser.error('Phase0 noise floor requires at least 16 perturbed members')
    if args.plan:
        print(json.dumps(dict(steps=['build frozen-source controls/default+strict','reproduce periodic baseline defect','CPU migration 1/2/4 ranks',
              'common-input CPU/CUDA+partition+generation RHS and 100 steps','GPU repeats/restart/forced retry','2D oblique+linear convergence','4 conditional noise ensembles'],
              gpu_count=4,physical_nodes=1,ensemble_members=args.ensemble_members,ensemble_time=args.ensemble_time,ensemble_grid=args.ensemble_grid,phase0_gate_complete=False),indent=2));return 0
    Suite(args).execute();return 0

if __name__=='__main__':raise SystemExit(main())
