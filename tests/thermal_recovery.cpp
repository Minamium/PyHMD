// Entropy/total-energy selection and nonoverlapping correction accounting.
// Build: c++ -O2 -std=c++17 -Iinclude tests/thermal_recovery.cpp -o build/thermal_recovery
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>

static void require(bool condition,const char* message) {
    if(!condition)throw std::runtime_error(message);
}
static bool close(double x,double y) {return std::abs(x-y)<1e-11*std::max(1.,std::max(std::abs(x),std::abs(y)));}
static bool budget_close(double x,double y) {return std::abs(x-y)<1e-10*std::max(std::abs(x),std::abs(y))+1e-30;}

int main(int argc,char** argv) {
    try {
        World w(argc,argv);require(w.size==1,"thermal_recovery is single process");
        auto p=setup_case_defaults(UNIFORM);p.nx=p.ny=8;p.spherical=p.log_r=0;
        p.xlo=p.ylo=0;p.xhi=p.yhi=1;p.entropy_switch=1;p.gm=0;
        Mesh mesh(p,w);auto m=mesh.view;int i=NG+4,j=NG+4,id=m.cell(i,j),k=4*m.nx+4;
        Primitive v;v[RHO]=1;v[ENE]=1;
        std::vector<State> source(m.tx*m.ty,conserved(v,1,p)),dest=source;
        std::vector<double> floors(m.nx*m.ny,0),energy(floors.size(),0);
        std::vector<ThermalAccounting> accounting(floors.size());
        auto reset=[&](){dest=source;std::fill(floors.begin(),floors.end(),0);std::fill(energy.begin(),energy.end(),0);std::fill(accounting.begin(),accounting.end(),ThermalAccounting{});};
        constexpr double weight=.25;
        dest[id][ENE]=-100;dest[id][ENTR]=2;double old=dest[id][ENE];
        require(recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Smooth entropy recovery rejected valid scalar");
        require(close(raw_pressure(dest[id],1,p),2),"Smooth recovery did not use advected entropy");
        require(accounting[k].entropy_events==1&&floors[k]==0,"Smooth recovery was misclassified as a floor");
        double change=weight*(dest[id][ENE]-old)*m.volume(i,j,p);
        require(budget_close(change,accounting[k].entropy_energy+energy[k]),"Smooth entropy energy budget does not close");

        reset();dest[id][ENE]=-p.p_floor/(p.gamma-1);dest[id][ENTR]=1e-16;old=dest[id][ENE];
        require(recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Entropy-floor recovery rejected valid state");
        require(raw_pressure(dest[id],1,p)>=p.p_floor*(1-1e-12),"Entropy floor postcondition failed");
        require(floors[k]==1&&accounting[k].entropy_events==1,"Entropy floor event counters are wrong");
        change=weight*(dest[id][ENE]-old)*m.volume(i,j,p);
        require(budget_close(change,accounting[k].entropy_energy+energy[k]),"Entropy/floor energy budgets overlap or do not close");

        // A compressive pressure jump is a shock detector input. Its gas heat
        // must come from total energy; the low advected entropy is reseeded.
        Primitive left=v,right=v;left[M1]=1;right[M1]=-1;right[ENE]=10;
        source[m.cell(i-1,j)]=conserved(left,1,p);source[m.cell(i+1,j)]=conserved(right,1,p);
        reset();v[RHO]=2;v[ENE]=3;dest[id]=conserved(v,1,p);dest[id][ENTR]=.01;
        require(recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Shock energy recovery rejected valid state");
        require(close(raw_pressure(dest[id],1,p),3),"Entropy selection erased shock gas heat");
        require(close(dest[id][ENTR],3*std::pow(2.,1-p.gamma)),"Shock entropy was not synchronized from gas energy");
        require(accounting[k].entropy_events==0,"Shock was incorrectly counted as entropy recovery");

        // Pressure-balanced density contacts must preserve total-energy gas
        // heat even at low beta. An invalid candidate must still be rejected.
        Primitive contact;contact[RHO]=1;contact[ENE]=1;contact[B3]=100;
        source.assign(source.size(),conserved(contact,1,p));contact[RHO]=.001;
        source[m.cell(i+1,j)]=conserved(contact,1,p);p.entropy_switch=2;
        reset();contact[RHO]=1;contact[ENE]=3;dest[id]=conserved(contact,1,p);dest[id][ENTR]=123;
        require(recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Contact total-energy recovery rejected a valid state");
        require(close(raw_pressure(dest[id],1,p),3),"Contact safeguard replaced gas heat by mixed entropy");
        require(close(dest[id][ENTR],3),"Contact entropy was not synchronized from total energy");
        require(accounting[k].entropy_events==0,"Contact incorrectly used entropy recovery");
        reset();dest[id][ENE]=-1;
        require(!recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Invalid contact energy was hidden by entropy recovery");

        // Internal-energy dual mode uses uAux as its low-thermal channel.
        // Finite negative values in an unused diagnostic channel are reseeded;
        // nonfinite or negative values in the selected channel are rejected.
        Primitive smooth;smooth[RHO]=1;smooth[ENE]=1;
        p.entropy_switch=3;source.assign(source.size(),conserved(smooth,1,p));reset();
        dest[id][ENE]=-100;dest[id][EINT]=2/(p.gamma-1);dest[id][ENTR]=-1;old=dest[id][ENE];
        require(recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Internal-energy mode failed to reseed unused negative entropy");
        require(close(raw_pressure(dest[id],1,p),2)&&close(dest[id][ENTR],2),"Internal-energy mode did not recover pressure and diagnostic entropy");
        require(accounting[k].internal_events==1&&accounting[k].entropy_events==0,"Internal-energy recovery was misclassified as entropy recovery");
        change=weight*(dest[id][ENE]-old)*m.volume(i,j,p);
        require(budget_close(change,accounting[k].internal_energy+energy[k]),"Internal-energy correction budget does not close");
        reset();dest[id][ENE]=-p.p_floor/(p.gamma-1);dest[id][EINT]=1e-16/(p.gamma-1);old=dest[id][ENE];
        require(recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Internal-energy floor recovery rejected valid channel");
        change=weight*(dest[id][ENE]-old)*m.volume(i,j,p);
        require(budget_close(change,accounting[k].internal_energy+energy[k]),"Internal-energy/floor correction budgets overlap");
        reset();dest[id][EINT]=-1;dest[id][ENTR]=-1;
        require(recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Valid total energy failed to reseed unused finite negative channels");
        require(dest[id][ENTR]>0&&dest[id][EINT]>0,"Unused channels were not synchronized from total energy");
        reset();dest[id][ENE]=-1;dest[id][EINT]=-1;
        require(!recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Selected negative internal-energy channel was accepted");
        reset();dest[id][EINT]=NAN;
        require(!recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Nonfinite auxiliary channel was hidden");
        p.entropy_switch=2;reset();dest[id][ENE]=-100;dest[id][ENTR]=2;dest[id][EINT]=-1;
        require(recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Entropy mode failed to reseed unused negative internal energy");
        require(close(raw_pressure(dest[id],1,p),2)&&close(dest[id][EINT],2/(p.gamma-1)),"Entropy recovery did not synchronize unused internal energy");

        p.entropy_switch=1;source.assign(source.size(),conserved(left,1,p));reset();dest[id][ENTR]=-1;
        require(!recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Negative entropy scalar was accepted in smooth selective mode");
        reset();dest[id][ENE]=INFINITY;
        require(!recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Nonfinite candidate was hidden by entropy recovery");
        reset();dest[id][RHO]=-1;
        require(!recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Negative density was hidden by entropy recovery");
        reset();p.entropy_switch=0;dest[id][ENE]=-1;
        require(!recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),accounting.data(),weight),"Negative total-energy pressure was silently floored instead of rejected");
        std::cout<<"PASS smooth/shock/contact thermal selection, correction accounting, and invalid-stage rejection\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL thermal recovery: "<<error.what()<<'\n';return 1;}
}
