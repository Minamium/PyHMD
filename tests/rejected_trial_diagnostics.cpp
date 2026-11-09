// A successful retry must retain its rejected-stage diagnostics and call the
// collective hook against an unchanged accepted checkpoint, without accepting
// invalid intermediate states or polluting the accepted correction budget.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>

static void require(bool value,const char* message) {
    if(!value)throw std::runtime_error(message);
}
static void initialize(Solver& s) {
    auto m=s.m;auto p=s.p;std::vector<State> u(s.u.n);std::vector<double>a(s.a.n,0);
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){Primitive v;v[RHO]=1;v[ENE]=.1;v[M1]=m.xc[i]<.5?-2:2;u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);}
    s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p);
}
int main(int argc,char** argv) {
    try {
        World w(argc,argv);require(w.size==1,"This diagnostic regression requires one CPU process");
        auto p=setup_case_defaults(UNIFORM);p.nx=64;p.ny=4;p.spherical=p.log_r=0;
        p.xlo=p.ylo=0;p.xhi=1;p.yhi=.125;p.gm=p.eta=0;p.eta_model=p.cooling=0;
        p.lower1=p.upper1=OUTFLOW;p.lower2=p.upper2=PERIODIC;p.entropy_switch=p.ct_energy_correction=0;
        Mesh mesh(p,w);Solver s(p,mesh.view,w);initialize(s);
        auto u=s.u.download();auto a=s.a.download();int hooks=0;int first_stage=0;
        s.rejected_trial_hook=[&](const Solver& x,const StageFailure& f,int attempt) {
            require(f.active==1&&f.stage>=1&&f.stage<=3,"Hook lost its failed RK stage");
            require(f.reason!=THERMAL_OK,"Hook lost the failure reason");
            require(x.last_trial_stages==f.stage,"Solver evaluated a later stage after rejection");
            require(x.rejected_steps==hooks,"Hook did not run before rejected_steps increment");
            require(attempt==hooks,"Hook attempt index is inconsistent");
            require(x.last_step_retries==attempt+1,"Per-step retry diagnostic is inconsistent");
            const auto b=x.base.download();const auto ba=x.abase.download();
            require(ba==a,"Hook checkpoint Psi differs from accepted state");
            for(size_t k=0;k<b.size();++k)for(int c=0;c<NV;++c)require(b[k][c]==u[k][c],"Hook checkpoint U differs from accepted state");
            for(const auto& t:x.thermal.download())require(t.ct_energy==0&&t.internal_energy==0&&t.entropy_energy==0,"Failed trial was committed before hook");
            if(hooks==0)first_stage=f.stage;
            ++hooks;
        };
        double requested=8*s.timestep(),accepted=s.step(requested);
        require(accepted<requested&&hooks>0,"Expansion did not exercise successful retry diagnostics");
        require(hooks==s.rejected_steps&&hooks==s.last_step_retries,"Accepted-step rejection counters disagree");
        require(s.first_rejected_stage.active&&s.first_rejected_stage.stage==first_stage,"First rejected trial was overwritten or lost");
        require(s.last_rejected_stage.active&&s.last_trial_stages==3,"Last rejected trial or accepted stage count is missing");
        require(first_stage<3,"Test did not exercise early-stage short circuit");
        double next=s.step(.1*s.timestep());(void)next;
        require(!s.first_rejected_stage.active&&!s.last_rejected_stage.active&&s.last_step_retries==0,"Rejected-stage diagnostics leaked into a clean next step");
        std::cout<<std::setprecision(17)<<"PASS successful retry hooks="<<hooks<<" first_stage="<<first_stage
                 <<" requested="<<requested<<" accepted="<<accepted<<"; clean next-step diagnostics reset\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL rejected-trial diagnostics: "<<error.what()<<'\n';return 1;}
}
