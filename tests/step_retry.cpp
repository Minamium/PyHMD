// A rejected large step must equal a clean flat/HLL step at the accepted dt.
// Build: c++ -O3 -std=c++17 -Iinclude tests/step_retry.cpp -o build/step_retry
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <chrono>
#include <filesystem>
#include <iomanip>

static void require(bool condition,const char* message) {
    if(!condition)throw std::runtime_error(message);
}

static void initialize_expansion(Solver& s) {
    auto m=s.m;auto p=s.p;
    std::vector<double> a(s.a.n,0);std::vector<State> u(s.u.n);
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i) {
        Primitive v;v[RHO]=1;v[ENE]=.1;v[M1]=m.xc[i]<.5?-2:2;
        u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
    }
    s.a.upload(a);s.u.upload(u);s.exchange(s.u.p,s.a.p);
}

static void check_no_budget(Solver& s) {
    for(auto* buffer:{&s.floors,&s.floor_mass,&s.floor_energy,&s.boundary_energy})
        for(double x:buffer->download())require(x==0,"Rejected trial contaminated accepted repair accounting");
    for(auto t:s.thermal.download())
        require(t.entropy_events==0&&t.entropy_energy==0&&t.internal_events==0&&t.internal_energy==0&&t.ct_energy==0
                &&t.compatible_cooling==0&&t.physical_ohmic==0&&t.negative_exchange==0,
                "Rejected trial contaminated accepted thermal accounting");
}

static double compare_physical(const Solver& a,const Solver& b) {
    auto x=a.u.download(),y=b.u.download();auto m=a.m;double maximum=0;
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i)
        for(int c=0;c<NV;++c)maximum=std::max(maximum,std::abs(x[m.cell(i,j)][c]-y[m.cell(i,j)][c]));
    return maximum;
}

int main(int argc,char** argv) {
    try {
        World w(argc,argv);require(w.size==1,"step_retry is single process");
        auto p=setup_case_defaults(UNIFORM);p.nx=64;p.ny=4;p.spherical=p.log_r=0;
        p.xlo=p.ylo=0;p.xhi=1;p.yhi=.125;p.gm=p.eta=0;p.eta_model=p.cooling=0;
        p.lower1=p.upper1=OUTFLOW;p.lower2=p.upper2=PERIODIC;
        p.entropy_switch=p.ct_energy_correction=0;
        Mesh mesh(p,w);Solver retried(p,mesh.view,w);initialize_expansion(retried);
        double requested=8*retried.timestep(),accepted=retried.step(requested);
        require(accepted>0&&accepted<requested,"Large expansion step did not reduce its timestep");
        require(retried.rejected_steps>0&&retried.flat_steps==1,"Retry/flat accepted-step counters are incorrect");
        require(retried.p.flat_reconstruction==0,"Retry reconstruction mode leaked into the next step");
        auto state=retried.u.download();auto m=retried.m;
        for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i) {
            auto q=state[m.cell(i,j)];
            require(q[RHO]>0&&raw_pressure(q,m.radius(i,j,p),p)>0,"Accepted retry has invalid raw density/pressure");
        }
        check_no_budget(retried);
        p.flat_reconstruction=1;Solver clean(p,mesh.view,w);initialize_expansion(clean);
        double clean_dt=clean.step(accepted);
        require(clean_dt==accepted&&clean.rejected_steps==0,"Clean flat reference could not accept the retry dt");
        double error=compare_physical(retried,clean);
        require(error<1e-13,"Retry differs from a clean flat/HLL step; rollback is incomplete");
        check_no_budget(clean);
        std::cout<<std::setprecision(17)<<"PASS retry requested="<<requested<<" accepted="<<accepted
                 <<" rejected_trials="<<retried.rejected_steps<<" rollback_reference_max_error="<<error<<'\n';

        // Exhausted retries must restore the original state and emit stage data.
        p.flat_reconstruction=0;p.max_step_retries=0;Solver exhausted(p,mesh.view,w);initialize_expansion(exhausted);
        auto initial=exhausted.u.download();auto initial_a=exhausted.a.download();
        auto temp=std::filesystem::temp_directory_path()/("pymhd_step_retry_failure_test_"+
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(temp);exhausted.failure_directory=temp.string();
        bool threw=false;try{exhausted.step(requested);}catch(const std::exception&){threw=true;}
        require(threw,"Exhausted invalid trial was silently accepted");
        require(std::filesystem::file_size(temp/"failed_stage_rank0.txt")>0,"Failed-stage diagnostic is missing");
        auto restored=exhausted.u.download();auto restored_a=exhausted.a.download();
        for(size_t n=0;n<initial.size();++n)for(int c=0;c<NV;++c)
            require(restored[n][c]==initial[n][c],"Exhausted step failed to restore the full conserved state");
        require(initial_a==restored_a,"Exhausted step failed to restore the flux function");
        check_no_budget(exhausted);
        std::filesystem::remove(temp/"failed_stage_rank0.txt");std::filesystem::remove(temp);
        std::cout<<"PASS exhausted-step rollback and failed-stage diagnostic\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL step retry: "<<error.what()<<'\n';return 1;}
}
