// A literature-comparison pressure-floor policy is explicit, nonconservative,
// and budgeted. Strict retry remains default; invalid density/nonfinite data
// and invalid actively selected auxiliary heat must never be hidden.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>

static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
int main(int argc,char**argv){try{
    World w(argc,argv);require(w.size==1,"Pressure-policy test requires one process");
    auto p=setup_case_defaults(UNIFORM);p.spherical=p.log_r=0;p.nx=p.ny=8;p.xlo=p.ylo=0;p.xhi=p.yhi=1;p.entropy_switch=0;
    require(p.negative_pressure_floor==0,"Negative pressure flooring is not default-strict");
    Mesh mesh(p,w);auto m=mesh.view;int i=NG+4,j=NG+4,id=m.cell(i,j),k=4*m.nx+4;constexpr double weight=.25;
    Primitive v;v[RHO]=2;v[ENE]=1;v[M1]=.3;v[M2]=-.2;v[M3]=.1;v[B1]=3;v[B2]=-2;v[B3]=1;
    std::vector<State>source(m.tx*m.ty,conserved(v,1,p)),dest=source;
    std::vector<double>floors(m.nx*m.ny),energy(floors.size());std::vector<ThermalAccounting>accounting(floors.size());
    auto reset=[&](){dest=source;std::fill(floors.begin(),floors.end(),0);std::fill(energy.begin(),energy.end(),0);std::fill(accounting.begin(),accounting.end(),ThermalAccounting{});};
    auto recover=[&](int&reason){return recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight,&reason);};
    reset();dest[id][ENE]=1;int reason=THERMAL_OK;State original=dest[id];
    require(raw_pressure(original,1,p)<0,"Fixture does not have finite negative gas pressure");
    require(!recover(reason)&&reason==THERMAL_NONPOSITIVE_ENERGY_PRESSURE,"Default policy silently floored a negative pressure");
    require(floors[k]==0&&energy[k]==0,"Rejected strict pressure changed the floor ledger");
    for(int c=0;c<NV;++c)require(dest[id][c]==original[c],"Rejected strict pressure mutated the candidate");

    p.negative_pressure_floor=1;
    for(double candidate_energy:{1.,7.14}) {
        reset();dest[id][ENE]=candidate_energy;double old=dest[id][ENE];
        require(raw_pressure(dest[id],1,p)<=0,"Floor fixture is not nonpositive");
        require(recover(reason)&&reason==THERMAL_OK,"Explicit floor policy rejected a finite negative/zero pressure");
        double actual=raw_pressure(dest[id],1,p),expected=weight*(dest[id][ENE]-old)*m.volume(i,j,p);
        require(actual>=p.p_floor*(1-1e-12),"Explicit floor failed the configured pressure postcondition");
        require(floors[k]==1&&energy[k]==expected,"Explicit floor energy ledger is not the actual weighted E change");
        require(accounting[k].entropy_events==0&&accounting[k].internal_events==0,"Energy floor was mislabeled as auxiliary recovery");
        require(dest[id][EINT]==actual/(p.gamma-1),"Energy-authoritative floor did not synchronize auxiliary heat");
        require(dest[id][ENTR]==actual*pow(dest[id][RHO],1-p.gamma),"Energy-authoritative floor did not synchronize diagnostic entropy");
    }
    for(int component:{ENE,EINT,RHO}) {
        reset();dest[id][component]=NAN;require(!recover(reason)&&reason==THERMAL_NONFINITE_STATE,"Floor policy hid a nonfinite state component");
        require(floors[k]==0&&energy[k]==0,"Nonfinite rejection polluted the floor ledger");
    }
    reset();dest[id][RHO]=0;require(!recover(reason)&&reason==THERMAL_NONPOSITIVE_DENSITY,"Floor policy hid a nonpositive density");
    reset();dest[id][RHO]=1e-300;dest[id][M1]=1e300;require(!recover(reason)&&reason==THERMAL_NONPOSITIVE_ENERGY_PRESSURE,"Floor policy hid nonfinite derived pressure");

    // An explicitly selected auxiliary channel is still required to be
    // positive. An energy floor is not a substitute for a valid smooth heat.
    p.entropy_switch=3;reset();dest[id][ENE]=-1;dest[id][EINT]=-1;
    require(!recover(reason)&&reason==THERMAL_NONPOSITIVE_AUXILIARY,"Floor policy hid an invalid selected internal-energy channel");

    // Conversely, a shock-authoritative cell may explicitly use the energy
    // floor and synchronize its untrusted (but finite) auxiliary variable.
    Primitive l=v,r=v;l[M1]=1;r[M1]=-1;r[ENE]=10;
    source[m.cell(i-1,j)]=conserved(l,1,p);source[m.cell(i+1,j)]=conserved(r,1,p);reset();dest[id][ENE]=-1;dest[id][EINT]=-1;
    require(excluded_thermal_cell(source.data(),i,j,m,p)!=0,"Shock fixture did not select total energy");
    require(recover(reason)&&floors[k]==1&&dest[id][EINT]>0,"Explicit shock-authoritative pressure floor was not applied/synchronized");
    std::cout<<std::setprecision(17)<<"PASS strict/default and explicit negative-pressure floor, exact E accounting, shock authority, invalid rho/nonfinite/active-aux rejection\n";return 0;
}catch(const std::exception&e){std::cerr<<"FAIL negative-pressure policy: "<<e.what()<<'\n';return 1;}}
