// A ghost magnetic boundary may not reset any physical CT flux or inject
// magnetic energy into active cells outside Faraday's update.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
int main(int argc,char**argv){try{World w(argc,argv);auto p=setup_case_defaults(DISK);p.nx=64;p.ny=32;p.inner_flux_project=0;
    Mesh mesh(p,w);auto m=mesh.view;Solver s(p,m,w);std::vector<double>a(s.a.n);std::vector<State>u(s.u.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=initial_potential(m.xf[i],m.yf[j],p);
    // The axis gauge is a legitimate physical constraint, separate from the
    // radial ghost boundary. Start in that constraint subspace exactly.
    for(int i=0;i<=m.tx;++i)a[m.vertex(i,NG)]=0;
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){auto v=initial(m.xc[i],m.yc[j],p);
        v[B1]=.5*(field1(a.data(),i,j,m,p)+field1(a.data(),i+1,j,m,p));v[B2]=.5*(field2(a.data(),i,j,m,p)+field2(a.data(),i,j+1,m,p));u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);}
    auto direct=a;paper_inner_potential(0,direct.data(),m,p);auto once=direct;paper_inner_potential(0,direct.data(),m,p);
    if(direct!=once)throw std::runtime_error("Ghost magnetic fill is not idempotent");
    s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p,1);auto out=s.a.download();auto state=s.u.download();
    for(int j=NG;j<=NG+m.ny;++j)for(int i=NG;i<=NG+m.nx;++i)if(out[m.vertex(i,j)]!=a[m.vertex(i,j)])throw std::runtime_error("Ghost boundary overwrote physical CT flux");
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i)for(int c=0;c<NV;++c)if(state[m.cell(i,j)][c]!=u[m.cell(i,j)][c])throw std::runtime_error("Ghost fill changed active conserved state");
    // The existing difference-of-squares ledger can contain roundoff even
    // when the actual conserved energy above is bitwise unchanged.
    double residual=0;for(double e:s.boundary_energy.download())residual+=fabs(e);
    if(residual>1e-15)throw std::runtime_error("Ghost fill injected resolved active magnetic projection energy");
    std::cout<<"PASS ghost boundary preserves every active Psi/state bit and is idempotent; absolute ledger roundoff="<<residual<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<"FAIL CT boundary contract: "<<e.what()<<'\n';return 1;}}
