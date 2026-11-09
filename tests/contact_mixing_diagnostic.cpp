// Diagnostic, not an acceptance test: expose entropy mixing at low-beta contacts.
// Build: c++ -O3 -std=c++17 -Iinclude tests/contact_mixing_diagnostic.cpp -o build/contact_mixing_diagnostic
#define main weak_contact_acceptance_main
#include "contact_advection.cpp"
#undef main

static void diagnostic(const World& w,double field,int mode,double end,bool guard) {
    auto p=setup_case_defaults(UNIFORM);p.nx=128;p.ny=4;p.spherical=p.log_r=0;
    p.xlo=p.ylo=0;p.xhi=1;p.yhi=.125;p.gm=p.eta=0;p.eta_model=p.cooling=0;
    p.lower1=p.upper1=p.lower2=p.upper2=PERIODIC;
    p.entropy_switch=mode;p.ct_energy_correction=1;
    if(!guard)p.entropy_contact_threshold=0;
    Mesh mesh(p,w);Solver s(p,mesh.view,w);initialize_contact(s,field);auto before=measure(s);
    double time=0;
    for(int step=1;end>0?time<end:step<=20;++step) {
        require(step<100000,"Contact diagnostic timestep collapsed");
        double dt=s.timestep();if(end>0)dt=std::min(dt,end-time);
        dt=s.step(dt);time+=dt;
        if(step!=1&&step!=20&&!(end>0&&time>=end))continue;
        auto after=measure(s);double entropy_energy=0,internal_energy=0,ct_energy=0;
        for(auto t:s.thermal.download()){entropy_energy+=t.entropy_energy;internal_energy+=t.internal_energy;ct_energy+=t.ct_energy;}
        std::cout<<"CONTACT_DIAGNOSTIC Bphi="<<field<<" beta="<<2*contact_pressure/(field*field)
                 <<" entropy_mode="<<mode<<" step="<<step<<" time="<<time
                 <<" contact_guard="<<(guard?1:0)
                 <<" max_pressure/initial="<<after.max_pressure/contact_pressure
                 <<" p_relative_error="<<after.pressure_error<<" velocity_error="<<after.velocity_error
                 <<" entropy_events="<<after.entropy_events<<" floors="<<after.floors
                 <<" rejected="<<s.rejected_steps<<" entropy_energy="<<entropy_energy
                 <<" internal_energy="<<internal_energy<<" ct_energy="<<ct_energy<<" energy_change="<<after.energy-before.energy<<'\n';
    }
}

int main(int argc,char** argv) {
    try {
        World w(argc,argv);require(w.size==1,"contact_mixing_diagnostic is single process");
        double end=argc>1?std::stod(argv[1]):0;
        bool guard=!(argc>2&&std::string(argv[2])=="--unguarded");
        std::cout<<std::setprecision(17)<<"DIAGNOSTIC ONLY: pressure errors below must be assessed; this executable does not assert accuracy.\n";
        for(double field:{10.,100.})for(int mode:{0,2,3})diagnostic(w,field,mode,end,guard);
        return 0;
    }catch(const std::exception& error){std::cerr<<"CONTACT_DIAGNOSTIC STOPPED: "<<error.what()<<'\n';return 1;}
}
