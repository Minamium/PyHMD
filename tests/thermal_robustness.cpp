// Physical low-beta regressions, independent of the disk initializer.
// Build: c++ -O3 -std=c++17 -Iinclude tests/thermal_robustness.cpp -o build/thermal_robustness
// Run:   build/thermal_robustness
// These short CPU checks do not certify a long-time disk or a CUDA execution.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>

static void require(bool condition,const std::string& message) {
    if(!condition)throw std::runtime_error(message);
}

static Parameters periodic_square(int n) {
    Parameters p=setup_case_defaults(UNIFORM);
    p.nx=p.ny=n;p.spherical=p.log_r=0;p.gm=0;
    p.xlo=p.ylo=0;p.xhi=p.yhi=1;
    p.lower1=p.upper1=p.lower2=p.upper2=PERIODIC;
    p.eta_model=0;p.eta=0;p.cooling=0;
    p.rho_floor=1e-10;p.p_floor=1e-12;
    p.internal_energy_sync=1;
    p.ct_emf_scheme=1;
    return p;
}

static void initialize_force_free(Solver& s,double pressure,double vx,double vy) {
    // Psi=A(cos kx+cos ky), Bz=k Psi implies curl B=k B, J cross B=0.
    // A constant velocity translates it. With constant eta, its exact amplitude
    // decays as exp(-eta k^2 t). Full Ohmic cooling leaves gas entropy unchanged.
    auto m=s.m;auto p=s.p;constexpr double k=2*PI,amplitude=1./(2*PI);
    std::vector<double> a(s.a.n);std::vector<State> u(s.u.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)
        a[m.vertex(i,j)]=amplitude*(std::cos(k*m.xf[i])+std::cos(k*m.yf[j]));
    const double sx=std::sin(k*(p.xhi-p.xlo)/(2*p.nx))/(k*(p.xhi-p.xlo)/(2*p.nx));
    const double sy=std::sin(k*(p.yhi-p.ylo)/(2*p.ny))/(k*(p.yhi-p.ylo)/(2*p.ny));
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i) {
        Primitive v;v[RHO]=1;v[ENE]=pressure;v[M1]=vx;v[M2]=vy;
        v[B1]=.5*(field1(a.data(),i,j,m,p)+field1(a.data(),i+1,j,m,p));
        v[B2]=.5*(field2(a.data(),i,j,m,p)+field2(a.data(),i,j+1,m,p));
        v[B3]=k*amplitude*(sx*std::cos(k*m.xc[i])+sy*std::cos(k*m.yc[j]));
        u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
    }
    s.a.upload(a);s.u.upload(u);s.exchange(s.u.p,s.a.p);
}

struct PhysicalMetrics {
    double min_entropy=1e300,max_entropy=0,min_p=1e300,max_p=0;
    double density_error=0,mode=0,energy=0,relative_divb=0;
    double floors=0,floor_energy=0,entropy_energy=0,internal_energy=0,internal_events=0,ct_energy=0,boundary_energy=0;
    double compatible_cooling=0,physical_ohmic=0,negative_exchange=0;
};

static PhysicalMetrics metrics(const Solver& s,double time,double vx,double vy) {
    auto m=s.m;auto p=s.p;auto u=s.u.download();auto a=s.a.download();
    PhysicalMetrics out;double mode_numerator=0,mode_denominator=0;
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i) {
        auto q=u[m.cell(i,j)];double pressure=raw_pressure(q,m.radius(i,j,p),p);
        require(q[RHO]>0&&pressure>0&&finite_value(pressure),"Invalid raw force-free state");
        const double entropy=pressure/std::pow(q[RHO],p.gamma),vol=m.volume(i,j,p);
        out.min_entropy=std::min(out.min_entropy,entropy);out.max_entropy=std::max(out.max_entropy,entropy);
        out.min_p=std::min(out.min_p,pressure);out.max_p=std::max(out.max_p,pressure);
        out.density_error=std::max(out.density_error,std::abs(q[RHO]-1));
        double mode=std::cos(2*PI*(m.xc[i]-vx*time))+std::cos(2*PI*(m.yc[j]-vy*time));
        mode_numerator+=vol*q[B3]*mode;mode_denominator+=vol*mode*mode;out.energy+=vol*q[ENE];
        double fl=m.area1(i,j,p)*field1(a.data(),i,j,m,p),fr=m.area1(i+1,j,p)*field1(a.data(),i+1,j,m,p);
        double fd=m.area2(i,j,p)*field2(a.data(),i,j,m,p),fu=m.area2(i,j+1,p)*field2(a.data(),i,j+1,m,p);
        double scale=std::abs(fl)+std::abs(fr)+std::abs(fd)+std::abs(fu);
        if(scale>0)out.relative_divb=std::max(out.relative_divb,std::abs(fr-fl+fu-fd)/scale);
    }
    out.mode=mode_numerator/mode_denominator;
    for(double x:s.floors.download())out.floors+=x;
    for(double x:s.floor_energy.download())out.floor_energy+=x;
    for(double x:s.boundary_energy.download())out.boundary_energy+=x;
    for(auto x:s.thermal.download()) {
        out.entropy_energy+=x.entropy_energy;out.internal_energy+=x.internal_energy;out.internal_events+=x.internal_events;out.ct_energy+=x.ct_energy;
        out.compatible_cooling+=x.compatible_cooling;out.physical_ohmic+=x.physical_ohmic;
        out.negative_exchange+=x.negative_exchange;
    }
    return out;
}

static double force_free_check(const World& w,bool resistive,int nx) {
    auto p=periodic_square(nx);constexpr double pressure=1e-8,vx=.2,vy=.13,end=.02;
    p.entropy_switch=3;p.ct_energy_correction=1;
    if(resistive){p.eta=.01;p.cooling=1;}
    Mesh mesh(p,w);Solver s(p,mesh.view,w);initialize_force_free(s,pressure,vx,vy);
    auto before=metrics(s,0,vx,vy);double time=0;int steps=0;
    while(time<end){double dt=std::min(s.timestep(),end-time);dt=s.step(dt);time+=dt;++steps;}
    auto after=metrics(s,time,vx,vy);
    double expected=std::exp(-p.eta*4*PI*PI*time);
    double relative_mode_error=std::abs(after.mode/before.mode-expected)/expected;
    double correction=after.floor_energy+after.entropy_energy+after.internal_energy+after.ct_energy+after.boundary_energy;
    double budget_defect=after.energy-before.energy+after.compatible_cooling-correction;
    std::cout<<(resistive?"resistive":"ideal")<<"_force_free nx="<<nx<<" steps="<<steps
             <<" min_entropy_ratio="<<after.min_entropy/pressure
             <<" max_entropy_ratio="<<after.max_entropy/pressure
             <<" rho_error="<<after.density_error<<" mode_error="<<relative_mode_error
             <<" floors="<<after.floors<<" floor_energy="<<after.floor_energy
             <<" energy_change="<<after.energy-before.energy
             <<" entropy_energy="<<after.entropy_energy<<" internal_energy="<<after.internal_energy<<" ct_energy="<<after.ct_energy
             <<" cooling="<<after.compatible_cooling<<" physical_ohmic="<<after.physical_ohmic
             <<" negative_exchange="<<after.negative_exchange<<" budget_defect="<<budget_defect
             <<" relative_divB="<<after.relative_divb<<'\n';
    require(after.relative_divb<2e-12,"Force-free CT divergence exceeds roundoff");
    require(after.floors==0,"Low-beta force-free wave required a pressure/density floor");
    require(after.internal_events>0,"Low-beta force-free case did not exercise internal-energy recovery");
    require(after.min_entropy/pressure>.99,"Force-free gas entropy spuriously decreased by more than 1 percent");
    require(after.max_entropy/pressure<1.01,"Force-free gas entropy spuriously increased by more than 1 percent");
    require(after.density_error<.02,"Force-free wave density differs by more than 2 percent");
    require(relative_mode_error<.003,"Force-free magnetic translation/decay mode error exceeds 0.3 percent");
    require(std::abs(budget_defect)<5e-11*std::max(1.,std::abs(before.energy)),"Periodic energy change is not explained by cooling and recorded corrections");
    return std::abs(correction);
}

static double isolated_ohmic_thermal_residual(const World& w,int nx) {
    auto p=periodic_square(nx);p.ny=4;constexpr double eta=.02;
    Mesh mesh(p,w);Solver s(p,mesh.view,w);auto m=s.m;
    std::vector<double> a(s.a.n,0);std::vector<State> u(s.u.n);
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i) {
        Primitive v;v[RHO]=1;v[ENE]=1;v[B3]=std::sin(2*PI*m.xc[i]);
        u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
    }
    s.a.upload(a);s.u.upload(u);s.evaluate(s.u.p,s.a.p);auto ideal=s.rhs.download();
    s.p.eta=eta;s.p.cooling=1;s.evaluate(s.u.p,s.a.p);auto cooled=s.rhs.download();
    double l1=0,volume=0,minimum=1e300,maximum=-1e300;
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i) {
        int id=m.cell(i,j);double b=u[id][B3];
        // At v=0, isolating eta-on minus eta-off removes ideal work. If all
        // Ohmic heat is radiated, the continuum thermal derivative is zero.
        // This diagnostic intentionally measures its finite-grid residual;
        // entropy recovery must not be mistaken for an exact energy identity.
        double residual=(cooled[id][ENE]-ideal[id][ENE])-b*(cooled[id][B3]-ideal[id][B3]);
        double vol=m.volume(i,j,p);l1+=vol*std::abs(residual);volume+=vol;
        minimum=std::min(minimum,residual);maximum=std::max(maximum,residual);
    }
    l1/=volume;
    std::cout<<"ohmic_residual nx="<<nx<<" L1="<<l1<<" min="<<minimum<<" max="<<maximum<<'\n';
    require(finite_value(l1),"Isolated Ohmic thermal residual is nonfinite");
    return l1;
}

static void check_ohmic_refinement(const World& w) {
    double coarse=isolated_ohmic_thermal_residual(w,32),fine=isolated_ohmic_thermal_residual(w,64);
    require(fine<1e-10||fine<.4*coarse,"Ohmic thermal truncation residual does not decrease under refinement");
    std::cout<<"PASS Ohmic residual refinement (diagnostic; nonzero values are not exact heat cancellation) ratio="<<fine/coarse<<'\n';
}

int main(int argc,char** argv) {
    try {
        World world(argc,argv);require(world.size==1,"thermal_robustness is single process");
        std::cout<<std::setprecision(17);
        check_ohmic_refinement(world);
        for(bool resistive:{false,true}) {
            double coarse=force_free_check(world,resistive,32),fine=force_free_check(world,resistive,64);
            require(fine<.6*coarse||fine<1e-12,"Thermal correction budget fails to decrease under refinement");
            std::cout<<"PASS force-free correction refinement ratio="<<fine/coarse<<'\n';
        }
        std::cout<<"PASS low-beta 2D ideal and radiatively cooled resistive force-free checks\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL thermal robustness: "<<error.what()<<'\n';return 1;}
}
