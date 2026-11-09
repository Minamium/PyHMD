// Independent analytic/geometric checks for the staggered poloidal CT update.
// Build: c++ -O2 -std=c++17 -Iinclude tests/ct_invariants.cpp -o build/ct_invariants
// Run:   build/ct_invariants
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>

struct CTDiagnostic {double absolute=0,relative=0;};

static CTDiagnostic ct_divergence(const Solver& solver) {
    const auto potential=solver.a.download();const auto m=solver.m;const auto p=solver.p;
    CTDiagnostic result;
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i) {
        const double left=m.area1(i,j,p)*field1(potential.data(),i,j,m,p);
        const double right=m.area1(i+1,j,p)*field1(potential.data(),i+1,j,m,p);
        const double down=m.area2(i,j,p)*field2(potential.data(),i,j,m,p);
        const double up=m.area2(i,j+1,p)*field2(potential.data(),i,j+1,m,p);
        const double residual=right-left+up-down;
        const double flux_scale=std::abs(right)+std::abs(left)+std::abs(up)+std::abs(down);
        result.absolute=std::max(result.absolute,std::abs(residual)/m.volume(i,j,p));
        if(flux_scale>0)result.relative=std::max(result.relative,std::abs(residual)/flux_scale);
    }
    return result;
}

static void require(bool condition,const std::string& message) {
    if(!condition)throw std::runtime_error(message);
}

static void initialize_flux(Solver& solver,int radial_power,double amplitude,double inflow=0) {
    const auto m=solver.m;const auto p=solver.p;
    std::vector<double> potential(solver.a.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)
        potential[m.vertex(i,j)]=amplitude*std::pow(m.xf[i],radial_power)*sq(std::sin(m.yf[j]));
    std::vector<State> state(solver.u.n);
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i) {
        Primitive v;v[RHO]=1;v[ENE]=1;v[M1]=inflow;
        v[B1]=.5*(field1(potential.data(),i,j,m,p)+field1(potential.data(),i+1,j,m,p));
        v[B2]=.5*(field2(potential.data(),i,j,m,p)+field2(potential.data(),i,j+1,m,p));
        state[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
    }
    solver.a.upload(potential);solver.u.upload(state);
    solver.exchange(solver.u.p,solver.a.p);
}

static void no_floor_events(const Solver& solver) {
    for(double count:solver.floors.download())require(count==0,"CT check triggered a positivity floor");
}

static void check_vertical_field(const World& world) {
    Parameters p=setup_case_defaults(UNIFORM);
    p.nx=64;p.ny=64;p.log_r=1;p.xlo=1;p.xhi=4;p.ylo=0;p.yhi=PI/2;
    p.lower2=AXIS;p.upper2=EQUATOR;p.eta=0;
    Mesh mesh(p,world);Solver solver(p,mesh.view,world);
    constexpr double bz=.02;
    // Psi = Bz R^2/2 is the exact uniform vertical field flux function.
    initialize_flux(solver,2,bz/2);
    const auto m=solver.m;const auto a=solver.a.download();const auto state=solver.u.download();
    const int radial=NG+m.nx/2,midplane=NG+m.ny;
    double parity_error=0;
    for(int g=1;g<=NG;++g) {
        parity_error=std::max(parity_error,std::abs(a[m.vertex(radial,NG-g)]-a[m.vertex(radial,NG+g)]));
        parity_error=std::max(parity_error,std::abs(a[m.vertex(radial,midplane+g)]-a[m.vertex(radial,midplane-g)]));
    }
    require(parity_error<1e-13,"Flux-function parity failed at the axis/equator");
    // Arithmetic averaging of face fields gives the known cos(dtheta/2)
    // factor, rather than an exactly point-sampled Cartesian Bz.
    const double expected_vertical=bz*std::cos((p.yhi-p.ylo)/(2*p.ny));
    double field_error=0;
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i) {
        const auto v=primitive(state[m.cell(i,j)],m.radius(i,j,p),p);
        const double radial_cartesian=v[B1]*std::sin(m.yc[j])+v[B2]*std::cos(m.yc[j]);
        const double vertical=v[B1]*std::cos(m.yc[j])-v[B2]*std::sin(m.yc[j]);
        field_error=std::max(field_error,std::max(std::abs(radial_cartesian),std::abs(vertical-expected_vertical)));
    }
    require(field_error<1e-12,"Uniform vertical flux did not reconstruct the analytically averaged field");
    const auto before=ct_divergence(solver);
    require(before.relative<5e-13,"Initial CT divergence exceeds roundoff scale");
    solver.step(std::min(1e-4,.1*solver.timestep()));solver.exchange(solver.u.p,solver.a.p);
    const auto after=ct_divergence(solver);
    require(after.relative<5e-13,"Evolved CT divergence exceeds roundoff scale");
    no_floor_events(solver);
    std::cout<<"PASS vertical-field parity="<<parity_error<<" field_error="<<field_error
             <<" initial_relative_divergence="<<before.relative<<" final_relative_divergence="<<after.relative<<"\n";
}

static void check_midplane_diffusion(const World& world) {
    Parameters p=setup_case_defaults(UNIFORM);
    p.nx=64;p.ny=64;p.log_r=1;p.xlo=1;p.xhi=4;p.ylo=0;p.yhi=PI/2;
    p.lower2=AXIS;p.upper2=EQUATOR;p.eta=.02;
    Mesh mesh(p,world);Solver solver(p,mesh.view,world);
    constexpr double amplitude=.01;
    // Psi = C r sin^2(theta): Jr=Jtheta=0, Jphi=2 C sin(theta)/r^2.
    // At the equator, Psi_t = -r eta Jphi = -2 C eta/r != 0.
    // This detects an incorrectly frozen equatorial flux function.
    initialize_flux(solver,1,amplitude);
    const auto m=solver.m;const auto before=solver.a.download();
    const double dt=std::min(1e-4,.1*solver.timestep());
    solver.step(dt);solver.exchange(solver.u.p,solver.a.p);
    const auto after=solver.a.download();const int j=NG+m.ny;
    double maximum_relative_rate_error=0,minimum_change=std::numeric_limits<double>::max();
    for(int i=NG+m.nx/4;i<NG+3*m.nx/4;++i) {
        const double expected=-2*amplitude*p.eta/m.xf[i];
        const double change=after[m.vertex(i,j)]-before[m.vertex(i,j)];
        const double measured=change/dt;
        require(change<0,"Resistive equatorial flux did not evolve in the expected direction");
        maximum_relative_rate_error=std::max(maximum_relative_rate_error,std::abs(measured-expected)/std::abs(expected));
        minimum_change=std::min(minimum_change,std::abs(change));
    }
    require(maximum_relative_rate_error<.05,"Resistive midplane rate differs from analytic -2 C eta/r by more than 5 percent");
    const auto divergence=ct_divergence(solver);
    require(divergence.relative<5e-13,"Resistive CT update violates divergence invariant");
    no_floor_events(solver);
    std::cout<<"PASS midplane-diffusion relative_rate_error="<<maximum_relative_rate_error
             <<" minimum_flux_change="<<minimum_change<<" relative_divergence="<<divergence.relative<<"\n";
}

int main(int argc,char** argv) {
    try {
        World world(argc,argv);
        require(world.size==1,"ct_invariants is a CPU single-process analytic test");
        std::cout<<std::setprecision(17);
        check_vertical_field(world);check_midplane_diffusion(world);
        return 0;
    }catch(const std::exception& error) {
        std::cerr<<"FAIL CT invariant: "<<error.what()<<"\n";
        return 1;
    }
}
