// Manufactured spherical resistive operators, including anisotropic eta(R).
// This evaluates the actual flux/CT/current operators without evolving gas.
// Analytic halos avoid conflating a prescribed boundary transient with an
// interior consistency error. Run on the single-process CPU backend.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <array>
#include <iomanip>

namespace {
constexpr double amplitude=.03,eta0=.07,anisotropy=.4;
enum Manufactured {POLOIDAL,TOROIDAL,POLOIDAL_VARIABLE,TOROIDAL_VARIABLE};
struct Norm {
    double error=0,reference=0;
    void add(double value,double exact,double scale,double volume) {
        error+=sq(value-exact)*volume;reference+=sq(scale)*volume;
    }
    double relative()const{return std::sqrt(error/reference);}
};
struct Result {double induction=0,heat=0,physical_heat=0,energy=0,cooled_energy=0;};
void require(bool condition,const std::string& what){if(!condition)throw std::runtime_error(what);}

Result evaluate(const World& world,int resolution,Manufactured test,int ct_correction) {
    auto p=setup_case_defaults(UNIFORM);
    p.nx=p.ny=resolution;p.xlo=1;p.xhi=4;p.ylo=0;p.yhi=PI/2;p.log_r=1;
    p.spherical=1;p.gm=0;p.eta=eta0;p.chi=anisotropy;
    p.lower2=AXIS;p.upper2=EQUATOR;p.entropy_switch=0;p.cooling=0;
    p.ct_energy_correction=ct_correction;
    const bool poloidal=test==POLOIDAL||test==POLOIDAL_VARIABLE;
    const bool variable=test==POLOIDAL_VARIABLE||test==TOROIDAL_VARIABLE;
    p.eta_model=variable?1:0;p.alpha_m=1;p.h=.4;
    Mesh mesh(p,world);Solver solver(p,mesh.view,world);const auto m=solver.m;
    std::vector<double> potential(solver.a.n,0);
    if(poloidal)for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)
        potential[m.vertex(i,j)]=amplitude*m.xf[i]*sq(std::sin(m.yf[j]));
    std::vector<State> state(solver.u.n);
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i) {
        Primitive v;v[RHO]=1;v[ENE]=1;
        v[B1]=.5*(field1(potential.data(),i,j,m,p)+field1(potential.data(),i+1,j,m,p));
        v[B2]=.5*(field2(potential.data(),i,j,m,p)+field2(potential.data(),i,j+1,m,p));
        if(!poloidal)v[B3]=amplitude*m.xc[i]*std::sin(m.yc[j]);
        state[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
    }
    solver.u.upload(state);solver.a.upload(potential);
    if(variable) {
        // On z/H<1, alpha=1, H=2hR, VA_mid=eta0/(2h) imply eta_P=eta0*R.
        // This intentionally prescribes a coefficient independent of the
        // manufactured field, rather than invoking a physical disk closure.
        std::vector<Midplane> mid(p.nx);
        for(auto& q:mid){q.va=eta0/(2*p.h);q.cs=.1;q.mu=.01;}
        solver.mid_global.upload(mid);solver.p.midplane=solver.mid_global.p;
    }
    p=solver.p;
    std::vector<State> full1(solver.flux1.n),full2(solver.flux2.n),res1(full1.size()),res2(full2.size());
    std::vector<double> full_a(solver.a.n),res_a(solver.a.n);
    for(int k=0;k<(m.nx+1)*m.ny;++k)face_flux(k,0,solver.u.p,solver.a.p,full1.data(),m,p,res1.data());
    for(int k=0;k<m.nx*(m.ny+1);++k)face_flux(k,1,solver.u.p,solver.a.p,full2.data(),m,p,res2.data());
    for(int k=0;k<(m.nx+1)*(m.ny+1);++k)corner_emf(k,solver.u.p,solver.a.p,full_a.data(),m,p,res_a.data());
    std::vector<State> rhs(solver.u.n),cooled(solver.u.n);
    std::vector<double> heat(m.nx*m.ny),physical(heat.size());
    Parameters cool=p;cool.cooling=1;
    for(int k=0;k<m.nx*m.ny;++k) {
        cell_rhs(k,solver.u.p,solver.a.p,res1.data(),res2.data(),rhs.data(),m,p,
                 res1.data(),res2.data(),res_a.data(),heat.data(),physical.data());
        cell_rhs(k,solver.u.p,solver.a.p,res1.data(),res2.data(),cooled.data(),m,cool,
                 res1.data(),res2.data(),res_a.data());
    }
    Norm ni,nh,np,ne,nc;
    for(int j=NG+3;j<NG+m.ny-3;++j)for(int i=NG+3;i<NG+m.nx-3;++i) {
        const double r=m.xc[i],theta=m.yc[j],st=std::sin(theta),ct=std::cos(theta);
        // A fixed physical band stays inside F=1 even for all neighboring
        // stencils at the coarsest grid. It excludes boundary-sensitive cells.
        if(r<1.4||r>3.0||theta<1.05||theta>1.35)continue;
        const int k=(j-NG)*m.nx+i-NG,id=m.cell(i,j);const double vol=m.volume(i,j,p);
        double induction,exact_induction,induction_scale,q,edot,mdot;
        if(poloidal) {
            // Psi=C r sin^2(theta), Jphi=2 C sin(theta)/r^2.
            // Its CT derivative is evaluated at the vertex, not at cell center.
            const double rv=m.xf[i],tv=m.yf[j];
            induction=res_a[m.vertex(i,j)];
            exact_induction=-2*amplitude*eta0*sq(std::sin(tv))/rv;
            induction_scale=exact_induction;
            q=4*sq(amplitude)*eta0*sq(st)/std::pow(r,4);
            mdot=sq(amplitude)*eta0/std::pow(r,4)*(-8*sq(ct)+2*sq(st));
            edot=mdot+q;
            if(variable) {
                // eta_P=eta0*R changes the angular shape of Psi_dot:
                // Psi_dot=-2 C eta0 sin^3(theta), Btheta_dot=0.
                exact_induction=-2*amplitude*eta0*std::pow(std::sin(tv),3);
                induction_scale=exact_induction;
                q=4*sq(amplitude)*eta0*std::pow(st,3)/std::pow(r,3);
                mdot=-12*sq(amplitude)*eta0*st*sq(ct)/std::pow(r,3);
                edot=mdot+q;
            }
        }else {
            // Bphi=kR has uniform axial current Jz=2k.
            const double eta=test==TOROIDAL_VARIABLE?eta0*r*st:eta0;
            const double tensor_factor=test==TOROIDAL_VARIABLE?anisotropy:1.;
            induction=rhs[id][B3];
            exact_induction=test==TOROIDAL_VARIABLE?2*amplitude*eta0*tensor_factor:0;
            induction_scale=test==TOROIDAL_VARIABLE?exact_induction:amplitude*eta0/r;
            q=4*sq(amplitude)*eta*tensor_factor;
            mdot=test==TOROIDAL_VARIABLE?2*sq(amplitude)*eta*tensor_factor:0;
            edot=mdot+q;
        }
        ni.add(induction,exact_induction,induction_scale,vol);
        nh.add(heat[k],q,q,vol);np.add(physical[k],q,q,vol);
        ne.add(rhs[id][ENE],edot,q,vol);nc.add(cooled[id][ENE],mdot,q,vol);
        require(std::abs(rhs[id][EINT]-heat[k])<1e-12*std::max(q,1e-30),"uncooled auxiliary heating differs from compatible exchange");
        require(std::abs(cooled[id][EINT])<1e-25,"instant Ohmic cooling leaves an auxiliary heat source");
    }
    require(ni.reference>0,"empty manufactured test mask");
    return {ni.relative(),nh.relative(),np.relative(),ne.relative(),nc.relative()};
}

void check_sequence(const World& world,Manufactured test,int correction) {
    std::array<Result,3> results;
    const char* name=test==POLOIDAL?"poloidal":test==TOROIDAL?"toroidal_harmonic":test==POLOIDAL_VARIABLE?"poloidal_eta_R":"toroidal_eta_R";
    for(int level=0;level<3;++level) {
        const int n=32<<level;results[level]=evaluate(world,n,test,correction);const auto r=results[level];
        std::cout<<name<<" ct_correction="<<correction<<" n="<<n
                 <<" induction="<<r.induction<<" heat="<<r.heat<<" physical_heat="<<r.physical_heat
                 <<" energy="<<r.energy<<" cooled_energy="<<r.cooled_energy<<'\n';
        for(double error:{r.induction,r.heat,r.physical_heat,r.energy,r.cooled_energy})
            require(std::isfinite(error)&&error<.15,std::string(name)+" operator has a large analytic error");
    }
    for(int component=0;component<5;++component) {
        auto value=[&](const Result& r){return component==0?r.induction:component==1?r.heat:component==2?r.physical_heat:component==3?r.energy:r.cooled_energy;};
        const double coarse=value(results[0]),medium=value(results[1]),fine=value(results[2]);
        require(fine<.005,std::string(name)+" finest operator error exceeds 0.5 percent");
        // Exact cancellations can hit roundoff; otherwise require measured
        // convergence toward the analytic resistive equation.
        if(coarse>1e-10&&medium>1e-10)
            require(medium<.45*coarse&&fine<.45*medium,std::string(name)+" does not converge near second order");
    }
}
} // namespace

int main(int argc,char** argv) {
    try {
        World world(argc,argv);require(world.size==1,"resistive_geometry requires one CPU process");
        std::cout<<std::setprecision(12);
        for(int correction:{0,1})for(auto test:{POLOIDAL,TOROIDAL,POLOIDAL_VARIABLE,TOROIDAL_VARIABLE})check_sequence(world,test,correction);
        std::cout<<"PASS spherical resistive geometry and heat convergence\n";return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL resistive geometry: "<<e.what()<<'\n';return 1;}
}
