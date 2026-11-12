// A single exchange must preserve arbitrary mean flux, before time integration.
#if __has_include("pymhd/problems/initialize.hpp")
#include "pymhd/problems/initialize.hpp"
#else
#include "pyhmd/problems/initialize.hpp"
#endif
#include <iostream>
int main(int argc,char**argv) {
    World world(argc,argv);
    try {
        auto p=setup_case_defaults(UNIFORM);
        p.nx=24;p.ny=16;p.spherical=p.log_r=0;p.xlo=p.ylo=0;p.xhi=p.yhi=1;p.gm=p.eta=0;
        p.lower1=p.upper1=p.lower2=p.upper2=PERIODIC;
        Mesh mesh(p,world);HostMesh host(mesh);auto m=host.view;Solver s(p,mesh.view,world);
        std::vector<double>a(s.a.n);std::vector<State>u(s.u.n);
        for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=2.3*m.yf[j]+.7*m.xf[i];
        Primitive v;v[RHO]=1;v[ENE]=.3;v[B1]=2.3;v[B2]=-.7;
        std::fill(u.begin(),u.end(),conserved(v,1,p));s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p);
        u=s.u.download();double error=0;
        for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i)
            error=std::max(error,std::max(fabs(u[m.cell(i,j)][B1]-2.3),fabs(u[m.cell(i,j)][B2]+.7)));
        error=world.maximum(error);
        if(world.rank==0)std::cout<<std::setprecision(17)<<"{\"ranks\":"<<world.size<<",\"max_mean_field_error\":"<<error<<",\"passed\":"<<(error<1e-12?"true":"false")<<"}\n";
        return error<1e-12?0:1;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';world.abort(2);return 2;}
}
