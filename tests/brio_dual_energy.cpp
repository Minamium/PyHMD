// Brio-Wu waveform comparison: modern dual energy must preserve the
// conservative gas heating in this reliable-thermal regime. The independent
// analytic shock test is parallel_shock.cpp; this is an operator comparison.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
static std::vector<State> evolve(const World&w,int n,bool modern) {
    auto p=setup_case_defaults(BRIO_WU);p.nx=n;p.ct_energy_correction=1;
    p.ct_emf_scheme=1;
    p.entropy_switch=modern?3:0;p.internal_energy_sync=modern?1:0;
    Mesh mesh(p,w);Solver s(p,mesh.view,w);auto m=s.m;
    std::vector<double>a(s.a.n);std::vector<State>u(s.u.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=initial_potential(m.xf[i],m.yf[j],p);
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i)u[m.cell(i,j)]=conserved(initial(m.xc[i],m.yc[j],p),1,p);
    s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p);double t=0;
    while(t<.1)t+=s.step(std::min(s.timestep(),.1-t));
    for(double f:s.floors.download())if(f!=0)throw std::runtime_error("Brio-Wu required a floor");
    u=s.u.download();std::vector<State>row(n);for(int i=0;i<n;++i)row[i]=u[m.cell(i+NG,NG+1)];return row;
}
int main(int argc,char**argv) {
    try{World w(argc,argv);for(int n:{128,256}){auto a=evolve(w,n,false),b=evolve(w,n,true);double error=0;
        for(int k=0;k<8;++k){double scale=1e-30,e=0;for(int i=0;i<n;++i){scale=std::max(scale,std::abs(a[i][k]));e=std::max(e,std::abs(a[i][k]-b[i][k]));}error=std::max(error,e/scale);}
        if(error>2e-10)throw std::runtime_error("Modern synchronization changed reliable Brio-Wu shock heating");
        std::cout<<"PASS Brio-Wu full waveform modern/conservative n="<<n<<" max component-normalized difference="<<error<<'\n';}
        return 0;
    }catch(const std::exception&e){std::cerr<<"FAIL Brio dual energy: "<<e.what()<<'\n';return 1;}
}
