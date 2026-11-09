// Oblique circular Alfven wave and weak magnetic-loop advection.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
static void require(bool x,const char*s){if(!x)throw std::runtime_error(s);}
static Parameters square(int n){auto p=setup_case_defaults(UNIFORM);p.nx=p.ny=n;p.spherical=p.log_r=0;p.xlo=p.ylo=0;p.xhi=p.yhi=1;
    p.gm=p.eta=0;p.entropy_switch=3;p.internal_energy_sync=1;p.ct_energy_correction=1;p.ct_emf_scheme=1;
    p.lower1=p.upper1=p.lower2=p.upper2=PERIODIC;return p;}
static double alfven(const World&w,int n){auto p=square(n);Mesh mesh(p,w);auto m=mesh.view;Solver s(p,m,w);
    constexpr double base=.7071067811865475244,amp=.1,pressure=1e-6,end=.2;
    double delta=PI/n,ae=amp*std::sin(delta)/delta*std::cos(delta);
    std::vector<double>a(s.a.n);std::vector<State>u(s.u.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=base*(m.yf[j]-m.xf[i])-amp/(2*PI*sqrt(2.))*std::sin(2*PI*(m.xf[i]+m.yf[j]));
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){Primitive v;v[RHO]=1;v[ENE]=pressure;
        v[B1]=.5*(field1(a.data(),i,j,m,p)+field1(a.data(),i+1,j,m,p));v[B2]=.5*(field2(a.data(),i,j,m,p)+field2(a.data(),i,j+1,m,p));
        v[B3]=ae*std::sin(2*PI*(m.xc[i]+m.yc[j]));v[M1]=base-v[B1];v[M2]=base-v[B2];v[M3]=-v[B3];u[m.cell(i,j)]=conserved(v,1,p);}
    s.u.upload(u);s.a.upload(a);double t=0;int steps=0;while(t<end){t+=s.step(std::min(s.timestep(),end-t));require(++steps<10000,"Wave CFL collapsed");}
    u=s.u.download();double error=0,rho=0,minp=INFINITY,maxp=0;
    for(int j=NG;j<NG+n;++j)for(int i=NG;i<NG+n;++i){auto q=u[m.cell(i,j)];double ph=2*PI*(m.xc[i]+m.yc[j]-sqrt(2.)*t);
        error+=sq(q[B1]-base+base*ae*std::cos(ph))+sq(q[B2]-base-base*ae*std::cos(ph))+sq(q[B3]-ae*std::sin(ph));
        rho=std::max(rho,fabs(q[RHO]-1));double ps=raw_pressure(q,1,p);minp=std::min(minp,ps);maxp=std::max(maxp,ps);}
    error=sqrt(error/(3*n*n))/ae;for(double f:s.floors.download())require(f==0,"Alfven wave needed floors");
    require(rho<.02&&minp>0&&maxp<1.03*pressure,"Oblique low-beta wave thermodynamics drifted");
    std::cout<<"oblique_alfven n="<<n<<" normalized_RMS="<<error<<" rho_error="<<rho<<" min_p="<<minp<<" max_p="<<maxp<<'\n';return error;
}
static double loop(const World&w,int n){auto p=square(n);Mesh mesh(p,w);auto m=mesh.view;Solver s(p,m,w);constexpr double amp=1e-3,rad=.2,vx=.3,vy=.2,end=.5;
    std::vector<double>a(s.a.n);std::vector<State>u(s.u.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i){double x=m.xf[i]-.5,y=m.yf[j]-.5;x-=std::round(x);y-=std::round(y);a[m.vertex(i,j)]=amp*std::max(0.,rad-std::hypot(x,y));}
    double initial_energy=0;
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){Primitive v;v[RHO]=v[ENE]=1;v[M1]=vx;v[M2]=vy;
        v[B1]=.5*(field1(a.data(),i,j,m,p)+field1(a.data(),i+1,j,m,p));v[B2]=.5*(field2(a.data(),i,j,m,p)+field2(a.data(),i,j+1,m,p));u[m.cell(i,j)]=conserved(v,1,p);
        if(i>=NG&&i<NG+n&&j>=NG&&j<NG+n)initial_energy+=(sq(v[B1])+sq(v[B2]))/2;}
    s.u.upload(u);s.a.upload(a);double t=0;while(t<end)t+=s.step(std::min(s.timestep(),end-t));u=s.u.download();a=s.a.download();double magnetic=0,error=0;
    for(int j=NG;j<NG+n;++j)for(int i=NG;i<NG+n;++i){auto q=u[m.cell(i,j)];magnetic+=(sq(q[B1])+sq(q[B2]))/2;double x=m.xf[i]-.5-vx*t,y=m.yf[j]-.5-vy*t;x-=std::round(x);y-=std::round(y);
        error+=fabs(a[m.vertex(i,j)]-amp*std::max(0.,rad-std::hypot(x,y)));require(q[RHO]>0&&raw_pressure(q,1,p)>0,"Loop has invalid gas state");}
    error/=n*n*amp;require(magnetic<1.05*initial_energy,"Advected magnetic loop amplified its energy");
    for(double f:s.floors.download())require(f==0,"Loop needed floors");
    std::cout<<"magnetic_loop n="<<n<<" normalized_L1="<<error<<" retained_magnetic_energy="<<magnetic/initial_energy<<'\n';return error;
}
int main(int argc,char**argv){try{World w(argc,argv);double a=alfven(w,24),b=alfven(w,48);require(b<.6*a,"Oblique wave failed grid convergence");a=loop(w,24);b=loop(w,48);require(b<.8*a,"Loop advection failed refinement");std::cout<<"PASS face UCT 2D oblique low-beta wave and loop advection\n";return 0;}catch(const std::exception&e){std::cerr<<"FAIL face UCT 2D: "<<e.what()<<'\n';return 1;}}
