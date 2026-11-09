// A cold rarefaction at the normal CFL step fails a high-order Euler stage,
// while donor-cell HLL accepts the same step. Spatial fallback must therefore
// be tried before reducing dt; the legacy halving policy remains selectable.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>

static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
static void initialize(Solver&s){auto m=s.m;auto p=s.p;std::vector<State>u(s.u.n);std::vector<double>a(s.a.n,0);
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){Primitive v;v[RHO]=1;v[ENE]=1e-10;v[M1]=.1*(m.xc[i]-.5);u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);}
    s.u.upload(u);s.a.upload(a);s.exchange(s.u.p,s.a.p);
}
static double compare(const Solver&a,const Solver&b){auto x=a.u.download(),y=b.u.download();auto m=a.m;double maximum=0;
    for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i)for(int c=0;c<NV;++c)maximum=std::max(maximum,std::abs(x[m.cell(i,j)][c]-y[m.cell(i,j)][c]));
    require(a.a.download()==b.a.download(),"Retry changed the flux function relative to the clean flat step");
    auto ab=a.flux_accounting.download(),bb=b.flux_accounting.download();
    for(size_t k=0;k<ab.size();++k)require(ab[k].mass_out==bb[k].mass_out&&ab[k].energy_out==bb[k].energy_out&&ab[k].angular_momentum_out==bb[k].angular_momentum_out&&ab[k].floor_potential==bb[k].floor_potential,"Rejected spatial trial polluted the accepted flux budget");
    return maximum;
}
int main(int argc,char**argv){try{
    World w(argc,argv);require(w.size==1,"Same-dt retry test requires one process");
    auto p=setup_case_defaults(UNIFORM);p.nx=32;p.ny=4;p.spherical=p.log_r=0;p.xlo=p.ylo=0;p.xhi=1;p.yhi=.125;
    p.lower1=p.upper1=OUTFLOW;p.lower2=p.upper2=PERIODIC;p.gm=p.eta=0;p.eta_model=p.cooling=p.entropy_switch=p.ct_energy_correction=0;
    p.max_step_retries=0;Mesh mesh(p,w);Solver ppm(p,mesh.view,w);initialize(ppm);double requested=ppm.timestep();
    bool rejected=false;try{ppm.step(requested);}catch(const std::exception&){rejected=true;}
    require(rejected&&ppm.first_rejected_stage.stage==1&&ppm.first_rejected_stage.reason==THERMAL_NONPOSITIVE_ENERGY_PRESSURE,"Fixture did not exercise a PPM pressure failure");
    auto flat_p=p;flat_p.flat_reconstruction=1;Solver flat(flat_p,mesh.view,w);initialize(flat);
    require(flat.step(requested)==requested&&flat.rejected_steps==0,"Clean flat solver cannot accept the same CFL step");

    auto retry_p=p;retry_p.max_step_retries=8;require(retry_p.retry_flat_same_dt==1,"Same-dt flat fallback is not the default");
    Solver retry(retry_p,mesh.view,w);initialize(retry);double accepted=retry.step(requested);
    require(accepted==requested&&retry.rejected_steps==1&&retry.flat_steps==1&&retry.last_step_retries==1,"PPM failure unnecessarily reduced dt before its successful flat fallback");
    require(retry.p.flat_reconstruction==0,"Flat fallback leaked into the configured high-order scheme");
    double error=compare(retry,flat);require(error==0,"Same-dt fallback differs from a clean flat solve");

    retry_p.retry_flat_same_dt=0;Solver legacy(retry_p,mesh.view,w);initialize(legacy);double old_accepted=legacy.step(requested);
    require(old_accepted==requested*.5&&legacy.rejected_steps==1,"Legacy halving-before-flat policy is not preserved");
    Solver flat_half(flat_p,mesh.view,w);initialize(flat_half);require(flat_half.step(old_accepted)==old_accepted,"Clean half-dt flat solve failed");
    require(compare(legacy,flat_half)==0,"Legacy retry differs from the clean half-dt flat reference");
    std::cout<<std::setprecision(17)<<"PASS same-dt flat retry requested="<<requested<<" accepted="<<accepted<<" legacy_accepted="<<old_accepted<<" clean_flat_state_error="<<error<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<"FAIL same-dt flat retry: "<<e.what()<<'\n';return 1;}}
