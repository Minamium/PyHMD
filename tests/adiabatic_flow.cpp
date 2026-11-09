// Analytic low-beta expansion trajectory plus expansion/compression PdV RHS.
// Build: c++ -O3 -std=c++17 -Iinclude tests/adiabatic_flow.cpp -o build/adiabatic_flow
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>

static void require(bool condition,const char* message) {
    if(!condition)throw std::runtime_error(message);
}

static void adiabatic_check(const World& w,double initial_gradient) {
    auto p=setup_case_defaults(ALFVEN);p.nx=128;p.ny=4;
    p.xlo=-1;p.xhi=1;p.ylo=0;p.yhi=.125;p.gm=p.eta=0;p.cooling=0;
    p.lower1=p.upper1=OUTFLOW;p.lower2=p.upper2=PERIODIC;
    p.entropy_switch=3;p.ct_energy_correction=1;
    p.internal_energy_sync=1;
    p.ct_emf_scheme=1;
    Mesh mesh(p,w);Solver s(p,mesh.view,w);auto m=s.m;
    constexpr double initial_pressure=1e-8;
    const double end=initial_gradient>0?.05:.01;
    std::vector<State> state(s.u.n);std::vector<double> potential(s.a.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)potential[m.vertex(i,j)]=m.yf[j];
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i) {
        Primitive v;v[RHO]=1;v[ENE]=initial_pressure;v[M1]=initial_gradient*m.xc[i];v[B1]=1;
        state[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
    }
    s.u.upload(state);s.a.upload(potential);s.exchange(s.u.p,s.a.p);
    // Instantaneous central RHS checks both signs without boundary evolution.
    // u_t = -div(u v) - p div(v) = -gamma*u*a for uniform gas energy.
    s.evaluate(s.u.p,s.a.p);auto initial_rhs=s.rhs.download();double rhs_error=0;
    const double expected_rhs=-p.gamma*initial_pressure/(p.gamma-1)*initial_gradient;
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i)
        if(std::abs(m.xc[i])<.25)rhs_error=std::max(rhs_error,std::abs(initial_rhs[m.cell(i,j)][EINT]/expected_rhs-1));
    require(rhs_error<1e-8,"Auxiliary internal RHS violates -gamma*u*div(v)");
    if(initial_gradient<0) {
        // OUTFLOW ghosts forbid the inflow required by a complete homologous
        // compression trajectory. Only its exact initial interior RHS is a
        // valid check with these production boundary conditions.
        std::cout<<std::setprecision(17)<<"adiabatic compression instantaneous_rhs_relative_error="<<rhs_error<<'\n';
        return;
    }
    double time=0;int steps=0;
    while(time<end){time+=s.step(std::min(s.timestep(),end-time));require(++steps<10000,"Adiabatic timestep collapsed");}
    const double scale=1+initial_gradient*time,expected_rho=1/scale,expected_pressure=initial_pressure/std::pow(scale,p.gamma);
    auto result=s.u.download();double density_error=0,pressure_error=0,velocity_error=0,entropy_error=0;int count=0;
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i) {
        // Outflow ghosts do not prescribe the exact time-dependent linear flow.
        // Assess only a central region remote from both radial boundaries.
        if(std::abs(m.xc[i])>.25)continue;
        auto q=result[m.cell(i,j)];double pressure=raw_pressure(q,m.radius(i,j,p),p);
        require(q[RHO]>0&&pressure>0&&finite_value(pressure),"Adiabatic state has invalid raw gas pressure");
        density_error=std::max(density_error,std::abs(q[RHO]/expected_rho-1));
        pressure_error=std::max(pressure_error,std::abs(pressure/expected_pressure-1));
        entropy_error=std::max(entropy_error,std::abs(pressure/std::pow(q[RHO],p.gamma)/initial_pressure-1));
        velocity_error=std::max(velocity_error,std::abs(q[M1]/q[RHO]-initial_gradient*m.xc[i]/scale));++count;
    }
    double floors=0,internal_events=0;for(double x:s.floors.download())floors+=x;
    for(auto x:s.thermal.download())internal_events+=x.internal_events;
    std::cout<<std::setprecision(17)<<"adiabatic gradient="<<initial_gradient<<" steps="<<steps<<" assessed_cells="<<count
             <<" rho_relative_error="<<density_error<<" p_relative_error="<<pressure_error
             <<" entropy_relative_error="<<entropy_error<<" rhs_relative_error="<<rhs_error
             <<" velocity_error="<<velocity_error<<" floors="<<floors<<'\n';
    require(count>0,"Adiabatic assessment region is empty");
    require(density_error<2e-4,"Homologous density differs from rho0/(1+a0t)");
    require(pressure_error<2e-3,"Auxiliary pressure violates the adiabatic gamma law");
    require(entropy_error<2e-3,"Smooth adiabatic auxiliary evolution spuriously changes entropy");
    require(velocity_error<2e-4,"Homologous velocity differs from a0*x/(1+a0t)");
    require(floors==0&&s.rejected_steps==0,"Homologous flow required floors or retries");
    require(internal_events>0,"Homologous expansion did not exercise internal-energy recovery");
}

int main(int argc,char** argv) {
    try {
        World w(argc,argv);require(w.size==1,"adiabatic_flow is single process");
        adiabatic_check(w,.4);adiabatic_check(w,-.4);
        std::cout<<"PASS auxiliary internal energy: adiabatic expansion trajectory and expansion/compression instantaneous RHS\n";return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL adiabatic flow: "<<error.what()<<'\n';return 1;}
}
