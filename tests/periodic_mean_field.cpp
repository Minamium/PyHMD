// Periodic flux functions carry arbitrary net B through quasi-periodic jumps.
// Neither jump may depend on the selected named test problem.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
int main(int argc,char**argv) {
    try{World w(argc,argv);auto p=setup_case_defaults(UNIFORM);p.nx=24;p.ny=16;p.spherical=p.log_r=0;
        p.xlo=p.ylo=0;p.xhi=p.yhi=1;p.gm=p.eta=0;p.lower1=p.upper1=p.lower2=p.upper2=PERIODIC;
        p.entropy_switch=3;p.internal_energy_sync=1;p.ct_energy_correction=1;
        Mesh mesh(p,w);auto m=mesh.view;Solver s(p,m,w);std::vector<State>u(s.u.n);std::vector<double>a(s.a.n);
        constexpr double bx=2.3,by=-.7;
        for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=bx*m.yf[j]-by*m.xf[i];
        Primitive v;v[RHO]=1;v[ENE]=.3;v[M1]=.2;v[M2]=-.13;v[B1]=bx;v[B2]=by;
        std::fill(u.begin(),u.end(),conserved(v,1,p));s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p);
        for(int n=0;n<20;++n)s.step(.5*s.timestep());u=s.u.download();double error=0;
        for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i){auto q=primitive(u[m.cell(i,j)],1,p);
            for(int k=0;k<NV;++k)if(k!=ENTR&&k!=EINT)error=std::max(error,std::abs(q[k]-v[k]));}
        if(error>1e-11)throw std::runtime_error("Arbitrary periodic mean field drifted");
        std::cout<<"PASS arbitrary periodic mean magnetic flux, max primitive drift="<<error<<'\n';return 0;
    }catch(const std::exception&e){std::cerr<<"FAIL periodic mean field: "<<e.what()<<'\n';return 1;}
}
