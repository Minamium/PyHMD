// Low-beta pressure-balanced contact with contact-protected dual recovery.
// Build: c++ -O3 -std=c++17 -Iinclude tests/strong_contact.cpp -o build/strong_contact
#define main weak_contact_acceptance_main
#include "contact_advection.cpp"
#undef main

struct ContactRun {ContactMetric end;double entropy_energy=0,internal_energy=0,internal_events=0,ct_energy=0,budget_defect=0;int steps=0;};
static ContactRun advance_contact(const World& w,double field,int mode) {
    auto p=setup_case_defaults(UNIFORM);p.nx=128;p.ny=4;p.spherical=p.log_r=0;
    p.xlo=p.ylo=0;p.xhi=1;p.yhi=.125;p.gm=p.eta=0;p.eta_model=p.cooling=0;
    p.lower1=p.upper1=p.lower2=p.upper2=PERIODIC;
    p.entropy_switch=mode;p.ct_energy_correction=1;
    Mesh mesh(p,w);Solver s(p,mesh.view,w);initialize_contact(s,field);auto before=measure(s);
    constexpr double end=.001;double time=0;ContactRun result;
    while(time<end) {
        time+=s.step(std::min(s.timestep(),end-time));
        require(++result.steps<10000,"Strong contact timestep collapsed");
    }
    result.end=measure(s);
    for(auto t:s.thermal.download()){result.entropy_energy+=t.entropy_energy;result.internal_energy+=t.internal_energy;result.internal_events+=t.internal_events;result.ct_energy+=t.ct_energy;}
    double boundary=0,floor=0;for(double x:s.boundary_energy.download())boundary+=x;for(double x:s.floor_energy.download())floor+=x;
    result.budget_defect=result.end.energy-before.energy-result.entropy_energy-result.internal_energy-result.ct_energy-boundary-floor;
    // Double-precision subtraction at beta=2e-6 amplifies an energy ulp in p.
    // The fixed pressure tolerance is compared with a total-energy reference,
    // rather than accepting arbitrary heating because the magnetic E is large.
    require(result.end.pressure_error<5e-8,"Strong contact spuriously heated/cooled beyond the pressure tolerance");
    require(result.end.velocity_error<1e-9,"Strong contact changed the uniform velocity");
    require(result.end.floors==0&&s.rejected_steps==0,"Strong contact required a floor or failed-step retry");
    require(std::abs(result.end.mass-before.mass)<1e-11*before.mass,"Strong periodic contact mass is not conserved");
    require(std::abs(result.budget_defect)<1e-12*std::abs(before.energy),"Strong contact energy change is not explained by recorded corrections and rounding tolerance");
    const double initial_gas_energy=contact_pressure/(p.gamma-1)*(p.xhi-p.xlo)*(p.yhi-p.ylo);
    require(std::abs(result.entropy_energy)<1e-8*initial_gas_energy,"Strong contact entropy correction is large relative to gas internal energy");
    require(std::abs(result.internal_energy)<1e-8*initial_gas_energy,"Strong contact auxiliary correction is large relative to gas internal energy");
    require(result.ct_energy==0,"Uniform toroidal field unexpectedly caused a CT energy correction");
    require(mode!=3||result.internal_events>0,"Strong contact did not exercise internal-energy recovery");
    return result;
}

int main(int argc,char** argv) {
    try {
        World w(argc,argv);require(w.size==1,"strong_contact is single process");std::cout<<std::setprecision(17);
        for(double field:{10.,100.}) {
            auto reference=advance_contact(w,field,0);
            for(int mode:{2,3}) {
            auto dual=advance_contact(w,field,mode);
            require(dual.end.pressure_error<10*reference.end.pressure_error+1e-10,"Contact-protected dual pressure error greatly exceeds total-energy reference");
            std::cout<<"PASS strong contact Bphi="<<field<<" beta="<<2*contact_pressure/(field*field)
                     <<" mode="<<mode<<" dual_steps="<<dual.steps<<" p_error="<<dual.end.pressure_error
                     <<" totalE_reference_p_error="<<reference.end.pressure_error
                     <<" velocity_error="<<dual.end.velocity_error<<" entropy_energy="<<dual.entropy_energy
                     <<" internal_energy="<<dual.internal_energy
                     <<" budget_defect="<<dual.budget_defect<<'\n';
            }
        }
        return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL strong contact: "<<error.what()<<'\n';return 1;}
}
