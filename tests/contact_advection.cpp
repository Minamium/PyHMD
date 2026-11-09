// Pressure-equilibrated contact: advected entropy mixing must not add gas heat.
// Build: c++ -O3 -std=c++17 -Iinclude tests/contact_advection.cpp -o build/contact_advection
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>

static void require(bool condition,const char* message) {
    if(!condition)throw std::runtime_error(message);
}

constexpr double contact_pressure=.01,vx=.2,vy=.13,bz=1e-4;
static void initialize_contact(Solver& s,double toroidal_field=bz) {
    auto m=s.m;auto p=s.p;std::vector<State> u(s.u.n);std::vector<double> a(s.a.n,0);
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i) {
        Primitive v;v[RHO]=m.xc[i]<.5?1:.001;v[ENE]=contact_pressure;
        v[M1]=vx;v[M2]=vy;v[B3]=toroidal_field;u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
    }
    s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p);
}
struct ContactMetric {double pressure_error=0,velocity_error=0,mass=0,energy=0,max_pressure=0,entropy_events=0,floors=0;};
static ContactMetric measure(const Solver& s) {
    auto m=s.m;auto p=s.p;auto u=s.u.download();ContactMetric out;
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i) {
        auto q=u[m.cell(i,j)];double pressure=raw_pressure(q,m.radius(i,j,p),p);
        require(q[RHO]>0&&pressure>0&&finite_value(pressure),"Contact has invalid raw density/pressure");
        out.pressure_error=std::max(out.pressure_error,std::abs(pressure-contact_pressure)/contact_pressure);
        out.max_pressure=std::max(out.max_pressure,pressure);
        out.velocity_error=std::max(out.velocity_error,std::max(std::abs(q[M1]/q[RHO]-vx),std::abs(q[M2]/q[RHO]-vy)));
        double vol=m.volume(i,j,p);out.mass+=q[RHO]*vol;out.energy+=q[ENE]*vol;
    }
    for(auto t:s.thermal.download())out.entropy_events+=t.entropy_events;
    for(double x:s.floors.download())out.floors+=x;
    return out;
}

int main(int argc,char** argv) {
    try {
        World w(argc,argv);require(w.size==1,"contact_advection is single process");
        auto p=setup_case_defaults(UNIFORM);p.nx=128;p.ny=4;p.spherical=p.log_r=0;
        p.xlo=p.ylo=0;p.xhi=1;p.yhi=.125;p.gm=p.eta=0;p.eta_model=p.cooling=0;
        p.lower1=p.upper1=p.lower2=p.upper2=PERIODIC;
        p.entropy_switch=2;p.ct_energy_correction=1;
        Mesh mesh(p,w);Solver dual(p,mesh.view,w);initialize_contact(dual);auto before=measure(dual);
        double time=0;int steps=0;constexpr double end=.1;
        while(time<end){double dt=dual.step(std::min(dual.timestep(),end-time));time+=dt;++steps;require(steps<10000,"Contact timestep collapsed");}
        auto after=measure(dual);
        require(after.pressure_error<1e-10,"Dual energy spuriously heated/cooled a pressure-balanced contact");
        require(after.velocity_error<1e-11,"Contact changed the uniform velocity");
        require(after.entropy_events==0,"High-thermal-fraction contact incorrectly used entropy recovery");
        require(after.floors==0&&dual.rejected_steps==0,"Pressure-balanced contact required floors or failed-step retries");
        require(std::abs(after.mass-before.mass)<1e-12*before.mass,"Periodic contact mass is not conserved");
        require(std::abs(after.energy-before.energy)<1e-11*before.energy,"Periodic contact energy is not conserved");
        std::cout<<std::setprecision(17)<<"PASS dual-energy pressure-balanced contact steps="<<steps
                 <<" p_relative_error="<<after.pressure_error<<" velocity_error="<<after.velocity_error
                 <<" entropy_events="<<after.entropy_events<<" floors="<<after.floors
                 <<" relative_energy_change="<<(after.energy-before.energy)/before.energy<<'\n';

        // Diagnostic only: passive entropy is nonlinear in density at fixed p.
        // Smooth-cell-only entropy selection may therefore heat a numerically
        // mixed contact. One step exposes it without running a thermal runaway.
        p.entropy_switch=1;Solver selective(p,mesh.view,w);initialize_contact(selective);
        double dt=selective.step(selective.timestep());auto one=measure(selective);
        std::cout<<"SELECTIVE contact one-step diagnostic dt="<<dt
                 <<" max_pressure/initial="<<one.max_pressure/contact_pressure
                 <<" pressure_relative_error="<<one.pressure_error
                 <<" entropy_events="<<one.entropy_events<<'\n';
        return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL contact advection: "<<error.what()<<'\n';return 1;}
}
