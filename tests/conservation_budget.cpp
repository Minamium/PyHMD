// End-to-end mass, angular momentum and E+rho*Phi budgets use the actual
// boundary numerical fluxes, cooling and accepted repair/projection terms.
// Include a density-floor case to expose its otherwise missing potential term.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>

static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
struct Integral {long double mass=0,energy=0,angular=0;};
static Integral integral(const Solver&s){Integral total;auto u=s.u.download();auto m=s.m;auto p=s.p;
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i){auto q=u[m.cell(i,j)];double vol=m.volume(i,j,p);total.mass+=q[RHO]*vol;total.energy+=(q[ENE]+q[RHO]*potential(m.xc[i],p))*vol;total.angular+=q[M3]*vol;}return total;}
static void run(const World&w,bool spherical,bool floor_case) {
    auto p=setup_case_defaults(UNIFORM);p.nx=24;p.ny=16;p.spherical=spherical;p.log_r=spherical;
    p.xlo=spherical?1:0;p.xhi=spherical?3:1;p.ylo=spherical?.3:0;p.yhi=spherical?1.2:1;
    p.lower1=p.upper1=p.lower2=p.upper2=OUTFLOW;p.gm=spherical?1:0;
    p.eta=.01;p.eta_model=0;p.cooling=1;p.ct_energy_correction=1;p.entropy_switch=0;
    if(floor_case)p.rho_floor=.01;
    Mesh mesh(p,w);Solver s(p,mesh.view,w);auto m=s.m;
    std::vector<double>a(s.a.n);std::vector<State>u(s.u.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=spherical?.02*sq(m.xf[i])*sq(sin(m.yf[j])):.02*m.yf[j];
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){Primitive v;v[RHO]=floor_case?.005:1+.1*sin(m.xc[i]);v[ENE]=1;
        v[M1]=floor_case?0:.15;v[M2]=floor_case?0:.03;v[M3]=floor_case?0:.2;
        v[B1]=.5*(field1(a.data(),i,j,m,p)+field1(a.data(),i+1,j,m,p));v[B2]=.5*(field2(a.data(),i,j,m,p)+field2(a.data(),i,j+1,m,p));v[B3]=.01*sin(m.xc[i])*sin(m.yc[j]);
        u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);}
    s.a.upload(a);s.u.upload(u);s.exchange(s.u.p,s.a.p);Integral before=integral(s);
    for(int n=0;n<5;++n)s.step(.04*s.timestep());Integral after=integral(s);
    long double massout=0,energyout=0,angularout=0,floorpot=0,floormass=0,floorenergy=0,boundary=0,ct=0,cooling=0,aux=0;
    for(auto b:s.flux_accounting.download()){massout+=b.mass_out;energyout+=b.energy_out;angularout+=b.angular_momentum_out;floorpot+=b.floor_potential;}
    for(double x:s.floor_mass.download())floormass+=x;
    for(double x:s.floor_energy.download())floorenergy+=x;
    for(double x:s.boundary_energy.download())boundary+=x;
    for(auto t:s.thermal.download()){ct+=t.ct_energy;cooling+=t.compatible_cooling;aux+=t.internal_energy+t.entropy_energy;}
    long double mr=after.mass-before.mass+massout-floormass;
    long double er=after.energy-before.energy+energyout+cooling-floorenergy-floorpot-boundary-ct-aux;
    long double lr=after.angular-before.angular+angularout;
    require(fabsl(mr)<2e-12L*std::max(1.L,fabsl(before.mass)),"Numerical mass budget does not close");
    require(fabsl(er)<2e-12L*std::max(1.L,fabsl(before.energy)),"Numerical total-energy budget does not close");
    require(fabsl(lr)<2e-12L*std::max(1.L,fabsl(before.angular)),"Numerical angular-momentum budget does not close");
    require(s.rejected_steps==0,"Budget fixture unexpectedly required retry");
    if(floor_case)require(floormass>0&&floorpot<0,"Density floor did not exercise negative gravitational-potential injection");
    std::cout<<std::setprecision(17)<<"PASS conservation spherical="<<spherical<<" floor="<<floor_case
             <<" residual_mass="<<static_cast<double>(mr)<<" energy="<<static_cast<double>(er)<<" angular="<<static_cast<double>(lr)
             <<" floor_potential="<<static_cast<double>(floorpot)<<'\n';
}
int main(int argc,char**argv){try{World w(argc,argv);require(w.size==1,"Conservation budget test requires one process");run(w,false,false);run(w,true,false);run(w,true,true);return 0;}catch(const std::exception&e){std::cerr<<"FAIL conservation budget: "<<e.what()<<'\n';return 1;}}
