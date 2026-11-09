#include "pymhd/setup.hpp"
#include "pymhd/checkpoint.hpp"
#include <chrono>
#include <filesystem>
#include <iostream>

static void require(bool ok,const char* what){if(!ok)throw std::runtime_error(what);}
static void initialize_wave(Solver& s) {
    const auto m=s.m;const auto p=s.p;std::vector<State> u(s.u.n);std::vector<double>a(s.a.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=m.yf[j];
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){auto v=initial(m.xc[i],m.yc[j],p);u[m.cell(i,j)]=conserved(v,1,p);}
    s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p);
}
static void equal(const Solver& a,const Solver& b) {
    auto x=a.u.download(),y=b.u.download();require(x.size()==y.size(),"state layout");
    for(std::size_t n=0;n<x.size();++n)for(int k=0;k<NV;++k)require(x[n][k]==y[n][k],"checkpoint changed conserved bits");
    require(a.a.download()==b.a.download(),"checkpoint changed vertex bits");
    require(a.floors.download()==b.floors.download()&&a.floor_mass.download()==b.floor_mass.download()
            &&a.floor_energy.download()==b.floor_energy.download()&&a.boundary_energy.download()==b.boundary_energy.download(),"checkpoint changed accepted repair accounting");
    auto ta=a.thermal.download(),tb=b.thermal.download();
    for(std::size_t n=0;n<ta.size();++n)require(std::memcmp(&ta[n],&tb[n],sizeof(ThermalAccounting))==0,"checkpoint changed thermal accounting");
    auto fa=a.flux_accounting.download(),fb=b.flux_accounting.download();
    for(std::size_t n=0;n<fa.size();++n)require(std::memcmp(&fa[n],&fb[n],sizeof(FluxAccounting))==0,"checkpoint changed numerical boundary flux accounting");
    require(a.rejected_steps==b.rejected_steps&&a.flat_steps==b.flat_steps,"checkpoint changed counters");
}
int main(int argc,char** argv) {
    try {
        World w(argc,argv);require(w.size==1,"checkpoint test uses one rank");
        auto p=setup_case_defaults(ALFVEN);p.nx=32;p.ny=4;p.eta=.002;p.ct_energy_correction=1;p.entropy_switch=3;
        Mesh m(p,w);Solver uninterrupted(p,m.view,w),resumed(p,m.view,w);initialize_wave(uninterrupted);initialize_wave(resumed);
        double t=0,dt=0;for(int n=0;n<5;++n){dt=uninterrupted.step(uninterrupted.timestep());t+=dt;}
        const auto dir=std::filesystem::temp_directory_path()/("pymhd_checkpoint_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        pymhd_checkpoint::Clock clock{t,dt,t+1,1,5};pymhd_checkpoint::write(uninterrupted,dir,clock);
        const auto restored=pymhd_checkpoint::read(resumed,dir,1);
        require(restored.time==t&&restored.dt==dt&&restored.step==5&&restored.next_output==t+1,"checkpoint changed integrator clock");equal(uninterrupted,resumed);
        for(int n=0;n<15;++n){auto a=uninterrupted.timestep(),b=resumed.timestep();require(a==b,"restart CFL differs");require(uninterrupted.step(a)==resumed.step(b),"restart accepted dt differs");}
        equal(uninterrupted,resumed);std::cout<<"PASS exact restart and continuation of resistive Alfven state/accounting\n";
        // Signature mismatch and corrupted bytes must reject before mutating s.
        auto saved=resumed.u.download();resumed.p.cfl=.3;bool rejected=false;
        try{pymhd_checkpoint::read(resumed,dir,1);}catch(const std::exception&){rejected=true;}
        require(rejected,"changed physics signature was accepted");resumed.p.cfl=p.cfl;
        auto shard=dir/pymhd_checkpoint::shard_name(0);{std::fstream f(shard,std::ios::in|std::ios::out|std::ios::binary);f.seekg(-16,std::ios::end);char ch=0;f.read(&ch,1);ch^=0x20;f.seekp(-16,std::ios::end);f.write(&ch,1);}
        rejected=false;try{pymhd_checkpoint::read(resumed,dir,1);}catch(const std::exception&){rejected=true;}
        require(rejected,"corrupted checkpoint was accepted");auto after=resumed.u.download();
        for(std::size_t n=0;n<saved.size();++n)for(int k=0;k<NV;++k)require(saved[n][k]==after[n][k],"failed load partially changed state");
        std::cout<<"PASS incompatible/corrupt checkpoint rejection is transactional\n";
        // Exercise the actual rejected-trial hook: save BASE while trial U is
        // invalid, reload the exact previous accepted state, and retain cause.
        auto rp=setup_case_defaults(UNIFORM);rp.nx=64;rp.ny=4;rp.spherical=rp.log_r=0;
        rp.xlo=rp.ylo=0;rp.xhi=1;rp.yhi=.125;rp.gm=rp.eta=0;rp.eta_model=rp.cooling=0;
        rp.lower1=rp.upper1=OUTFLOW;rp.lower2=rp.upper2=PERIODIC;rp.entropy_switch=rp.ct_energy_correction=0;
        Mesh rm(rp,w);Solver trial(rp,rm.view,w),from_reject(rp,rm.view,w);
        std::vector<State> expanding(trial.u.n);std::vector<double> zero(trial.a.n,0);
        for(int j=0;j<trial.m.ty;++j)for(int i=0;i<trial.m.tx;++i){Primitive v;v[RHO]=1;v[ENE]=.1;v[M1]=trial.m.xc[i]<.5?-2:2;expanding[trial.m.cell(i,j)]=conserved(v,1,rp);}
        trial.u.upload(expanding);trial.a.upload(zero);trial.exchange(trial.u.p,trial.a.p);const auto accepted_before=trial.u.download();
        bool captured=false;const auto retry_dir=dir/"retry";
        trial.rejected_trial_hook=[&](const Solver& rejected_solver,const StageFailure& failure,int attempt){if(!captured){captured=true;pymhd_checkpoint::write(rejected_solver,retry_dir,{0,0,1,1,0},true,&failure,attempt);}};
        const double oversized=8*trial.timestep();trial.step(oversized);require(captured,"retry fixture never failed");
        pymhd_checkpoint::read(from_reject,retry_dir,1);auto loaded=from_reject.u.download();
        for(std::size_t n=0;n<loaded.size();++n)for(int k=0;k<NV;++k)require(loaded[n][k]==accepted_before[n][k],"retry checkpoint stored failed trial instead of accepted base");
        require(std::filesystem::file_size(retry_dir/"rejected_stage_rank0.txt")>0,"successful retry lost failed-stage record");
        from_reject.rejected_trial_hook=[](const Solver&,const StageFailure&,int){throw std::runtime_error("simulated checkpoint I/O error");};
        rejected=false;try{from_reject.step(oversized);}catch(const std::exception&){rejected=true;}
        require(rejected&&from_reject.p.flat_reconstruction==0,"checkpoint I/O error did not abort/reset reconstruction policy");loaded=from_reject.u.download();
        for(std::size_t n=0;n<loaded.size();++n)for(int k=0;k<NV;++k)require(loaded[n][k]==accepted_before[n][k],"checkpoint I/O error did not rollback physical state");
        std::cout<<"PASS successful-retry checkpoint records accepted base and cause; I/O exception rolls back\n";
        std::filesystem::remove_all(dir);return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL checkpoint: "<<e.what()<<'\n';return 1;}
}
