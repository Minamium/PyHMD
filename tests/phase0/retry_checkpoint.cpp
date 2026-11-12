// Force a rejected trial using the public stepper; write all accepted ledgers.
#if __has_include("pymhd/io/checkpoint.hpp")
#include "pymhd/problems/initialize.hpp"
#include "pymhd/io/checkpoint.hpp"
#else
#include "pyhmd/problems/initialize.hpp"
#include "pyhmd/io/checkpoint.hpp"
#endif
int main(int argc,char**argv) {
    World world(argc,argv);
    try {
        if(argc!=2)throw std::runtime_error("retry_checkpoint OUTPUT");
        auto p=setup_case_defaults(UNIFORM);p.nx=64;p.ny=4;p.spherical=p.log_r=0;
        p.xlo=p.ylo=0;p.xhi=1;p.yhi=.125;p.gm=p.eta=0;p.eta_model=p.cooling=0;
        p.lower1=p.upper1=OUTFLOW;p.lower2=p.upper2=PERIODIC;p.entropy_switch=p.ct_energy_correction=0;
        Mesh mesh(p,world);HostMesh h(mesh);Solver s(p,mesh.view,world);
        std::vector<State>u(s.u.n);std::vector<double>a(s.a.n,0);
        for(int j=0;j<h.view.ty;++j)for(int i=0;i<h.view.tx;++i){Primitive v;v[RHO]=1;v[ENE]=.1;v[M1]=h.xc[i]<.5?-2:2;u[h.view.cell(i,j)]=conserved(v,1,p);}
        s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p);
        double requested=8*s.timestep(),accepted=s.step(requested);
        if(!(s.rejected_steps>0&&accepted<requested&&s.flat_steps==1))throw std::runtime_error("retry not exercised");
        pymhd_checkpoint::write(s,argv[1],{accepted,accepted,1,1,1});
        if(world.rank==0)std::cout<<"PASS forced retry checkpoint; rejected="<<s.rejected_steps<<'\n';
        return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';world.abort(1);return 1;}
}
