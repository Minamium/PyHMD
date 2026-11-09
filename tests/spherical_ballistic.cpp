// Exact cold, rotating, magnetized spherical-coordinate trajectory.
// Cartesian horizontal ballistic expansion: D=1+Omega0^2 t^2,
// rho=rho0/D, p=p0/D^gamma, vR=Omega0^2 t R/D, vphi=Omega0 R/D,
// vz=constant, Bz=Bz0/D. Uniform pressure/B have zero force and curl B=0.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
static double run(const World&w,int n,int ct){auto p=setup_case_defaults(UNIFORM);p.nx=p.ny=n;p.spherical=p.log_r=1;p.xlo=1;p.xhi=5;p.ylo=.2;p.yhi=1.4;
    p.gm=p.eta=0;p.entropy_switch=3;p.internal_energy_sync=1;p.ct_emf_scheme=ct;p.ct_energy_correction=1;
    p.lower1=p.upper1=p.lower2=p.upper2=OUTFLOW;Mesh mesh(p,w);auto m=mesh.view;Solver s(p,m,w);
    constexpr double omega=.4,vz=.2,bz=.1,p0=1e-8,end=.2;
    std::vector<double>a(s.a.n);std::vector<State>u(s.u.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=.5*bz*sq(m.xf[i]*std::sin(m.yf[j]));
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){Primitive v;v[RHO]=1;v[ENE]=p0;v[M1]=vz*std::cos(m.yc[j]);v[M2]=-vz*std::sin(m.yc[j]);v[M3]=omega*m.radius(i,j,p);
        v[B1]=.5*(field1(a.data(),i,j,m,p)+field1(a.data(),i+1,j,m,p));v[B2]=.5*(field2(a.data(),i,j,m,p)+field2(a.data(),i,j+1,m,p));u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);}
    s.u.upload(u);s.a.upload(a);double t=0;int steps=0;while(t<end){t+=s.step(std::min(s.timestep(),end-t));if(++steps>100000)throw std::runtime_error("Ballistic CFL collapsed");}
    double d=1+omega*omega*t*t,exrho=1/d,expres=p0/std::pow(d,p.gamma),exbz=bz/d*std::cos((p.yhi-p.ylo)/(2*n));
    u=s.u.download();double rhoerr=0,perr=0,verr=0,berr=0;int count=0;
    for(int j=NG;j<NG+n;++j)for(int i=NG;i<NG+n;++i){double r=m.xc[i],th=m.yc[j];if(r<2||r>4||th<.6||th>1.)continue;
        auto q=u[m.cell(i,j)];auto v=primitive(q,m.radius(i,j,p),p);double R=r*std::sin(th),vR=omega*omega*t*R/d;
        rhoerr=std::max(rhoerr,fabs(v[RHO]/exrho-1));perr=std::max(perr,fabs(raw_pressure(q,R,p)/expres-1));
        verr=std::max(verr,std::max({fabs(v[M1]-vR*std::sin(th)-vz*std::cos(th)),fabs(v[M2]-vR*std::cos(th)+vz*std::sin(th)),fabs(v[M3]-omega*R/d)}));
        berr=std::max(berr,std::hypot(v[B1]-exbz*std::cos(th),v[B2]+exbz*std::sin(th))/bz);++count;}
    for(double f:s.floors.download())if(f)throw std::runtime_error("Ballistic trajectory needed a floor");
    if(!count||rhoerr>.01||perr>.02||verr>.02||berr>.02)throw std::runtime_error("Spherical exact trajectory error too large");
    double error=rhoerr+perr+verr+berr;
    std::cout<<"spherical_ballistic n="<<n<<" ct="<<ct<<" steps="<<steps<<" rho_rel="<<rhoerr<<" p_rel="<<perr<<" velocity_abs="<<verr<<" B_rel="<<berr<<'\n';return error;
}
int main(int argc,char**argv){try{World w(argc,argv);for(int ct:{0,1}){double c=run(w,24,ct),f=run(w,48,ct);if(f>.65*c)throw std::runtime_error("Cold rotating spherical trajectory does not converge");}
    std::cout<<"PASS analytic cold rotating magnetized spherical trajectory and refinement\n";return 0;}catch(const std::exception&e){std::cerr<<"FAIL spherical ballistic: "<<e.what()<<'\n';return 1;}}
