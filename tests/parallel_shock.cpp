// Exact stationary Rankine-Hugoniot gas shock with a constant normal B field.
// Normal B must not change its physical rho2=3, p2=11 at upstream Mach 3.
// B=200 is deliberately outside the modern method's magnetic reliability gate;
// report its error as a limitation instead of calling a positive run a pass.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
static void require(bool x,const char* message){if(!x)throw std::runtime_error(message);}
struct Result {double l1p=0,post_p=0,floor=0,aux=0;int steps=0;};
static Result run(const World&w,int n,double bn,int policy) {
    auto p=setup_case_defaults(UNIFORM);p.nx=n;p.ny=3;p.spherical=p.log_r=0;
    p.xlo=p.ylo=0;p.xhi=1;p.yhi=.1;p.gm=p.eta=0;p.entropy_switch=3;p.internal_energy_sync=policy;
    p.lower1=p.upper1=OUTFLOW;p.lower2=p.upper2=PERIODIC;p.ct_energy_correction=1;
    p.ct_emf_scheme=1;
    Mesh mesh(p,w);Solver s(p,mesh.view,w);auto m=mesh.view;
    std::vector<State> u(s.u.n);std::vector<double>a(s.a.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=bn*m.yf[j];
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){Primitive v;bool upstream=m.xc[i]<.7;
        v[RHO]=upstream?1:3;v[ENE]=upstream?1:11;v[M1]=sqrt(15.)/v[RHO];
        v[B1]=.5*(field1(a.data(),i,j,m,p)+field1(a.data(),i+1,j,m,p));
        u[m.cell(i,j)]=conserved(v,1,p);}
    s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p);double t=0;Result out;
    while(t<.08){t+=s.step(std::min(s.timestep(),.08-t));++out.steps;}
    u=s.u.download();double count=0;
    for(int i=NG;i<NG+m.nx;++i){double pressure=raw_pressure(u[m.cell(i,NG+1)],1,p),x=m.xc[i];
        require(std::isfinite(pressure)&&pressure>0,"Parallel shock has invalid pressure");
        // Lower OUTFLOW is a diode and cannot prescribe the analytic inflow.
        // Its disturbance remains at x<.5 by t=.08; exclude that region.
        if(x>.5)out.l1p+=fabs(pressure-(x<.7?1:11))/(.5*n);
        if(x>.74&&x<.84){out.post_p+=pressure;count+=1;}}
    out.post_p/=count;
    for(double x:s.floors.download())out.floor+=x;
    for(auto x:s.thermal.download())out.aux+=x.internal_energy;
    require(out.floor==0,"Parallel shock needed floors");
    std::cout<<std::setprecision(12)<<"parallel_shock n="<<n<<" Bn="<<bn<<" sync="<<(policy?"modern":"legacy")
        <<" steps="<<out.steps<<" L1p="<<out.l1p<<" post_p="<<out.post_p<<" exact_post_p=11 aux_dE="<<out.aux<<'\n';
    return out;
}
int main(int argc,char**argv) {
    try{World w(argc,argv);
        for(double bn:{0.,20.}) {
            auto coarse=run(w,96,bn,1),fine=run(w,192,bn,1);
            require(fine.l1p<coarse.l1p,"Within-scope shock does not converge with resolution");
            require(fabs(fine.post_p/11-1)<.04,"Within-scope shock lost downstream heat");
            auto conservative=run(w,192,bn,0);
            require(fabs(fine.post_p-conservative.post_p)<.02,"Modern changed reliable shock heating");
        }
        auto limited=run(w,96,200,1);auto reference=run(w,96,200,0);
        std::cout<<"LIMITATION extreme normal-field shock: modern_post_error="<<fabs(limited.post_p/11-1)
            <<" legacy_post_error="<<fabs(reference.post_p/11-1)<<"; positive pressure alone is not physical validation\n";
        std::cout<<"PASS within-scope analytic parallel shock heating and convergence\n";return 0;
    }catch(const std::exception&e){std::cerr<<"FAIL parallel shock: "<<e.what()<<'\n';return 1;}
}
