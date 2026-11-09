// Published Enzo-E modern synchronization: all three gates, not a shock-veto
// override. Independently check energy accounting and invalid active channels.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include "pymhd/checkpoint.hpp"
static void require(bool x,const char* message){if(!x)throw std::runtime_error(message);}
static bool close(double a,double b){return std::abs(a-b)<1e-11*std::max({1.,std::abs(a),std::abs(b)});}
int main(int argc,char** argv) {
    try {
        World w(argc,argv);auto p=setup_case_defaults(UNIFORM);p.nx=p.ny=8;p.spherical=p.log_r=0;
        p.xlo=p.ylo=0;p.xhi=p.yhi=1;p.entropy_switch=3;p.internal_energy_sync=1;p.gm=p.eta=0;
        Mesh mesh(p,w);auto m=mesh.view;int i=NG+4,j=NG+4,id=m.cell(i,j),k=4*m.nx+4;
        Primitive v;v[RHO]=1;v[ENE]=1;v[B1]=100;
        auto q=conserved(v,1,p);
        require(!modern_energy_reliable(q,1,p),"Magnetic reliability gate missing");
        v[B1]=0;v[M1]=100;q=conserved(v,1,p);
        require(!modern_energy_reliable(q,1,p),"Kinetic reliability gate missing");
        v[M1]=0;q=conserved(v,1,p);q[EINT]=4*q[EINT];
        require(!modern_energy_reliable(q,1,p),"Half-aux reliability gate missing");
        q=conserved(v,1,p);require(modern_energy_reliable(q,1,p),"Reliable thermal state rejected");
        std::vector<State> source(m.tx*m.ty,q),dest=source;
        std::vector<double> floors(m.nx*m.ny),energy(floors.size());
        std::vector<ThermalAccounting> account(floors.size());
        std::vector<int> flags(m.tx*m.ty,1);p.shock_flags=flags.data();
        // The original neighbour-veto failure: positive transported heat, tiny
        // negative E-K-M. Modern must select aux regardless of stale flags.
        Primitive low;low[RHO]=1;low[ENE]=1e-8;low[B1]=10;low[M1]=1;
        q=conserved(low,1,p);source.assign(source.size(),q);dest=source;
        dest[id][ENE]-=2e-8/(p.gamma-1);double old=dest[id][ENE];int reason=0;
        require(recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),account.data(),.25,&reason),"Modern failed the neighbour-veto regression");
        require(std::abs(raw_pressure(dest[id],1,p)-1e-8)<2e-14,"Modern pressure differs from transported heat");
        require(account[k].internal_events==1&&floors[k]==0,"Modern recovery was hidden as a floor");
        require(close(account[k].internal_energy,.25*(dest[id][ENE]-old)*m.volume(i,j,p)),"Modern energy ledger is inconsistent");
        p.internal_energy_sync=0;dest=source;dest[id][ENE]-=2e-8/(p.gamma-1);
        require(!recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),account.data(),.25,&reason),"Legacy control did not reject the captured failure path");
        p.internal_energy_sync=1;dest=source;dest[id][EINT]=-1;
        require(!recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),account.data(),.25,&reason),"Selected invalid auxiliary energy was accepted");
        dest=source;dest[id][EINT]=NAN;
        require(!recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),account.data(),.25,&reason),"Nonfinite state was concealed");
        // A resolved thermal-energy-dominated shock must reseed aux from E.
        v[RHO]=3;v[ENE]=11;v[B1]=20;v[M1]=sqrt(15.)/3;dest[id]=conserved(v,1,p);dest[id][EINT]=1;
        require(modern_energy_reliable(dest[id],1,p),"Within-scope Mach3 postshock failed reliability gates");
        require(recover_thermal(k,dest.data(),source.data(),m,p,floors.data(),energy.data(),account.data(),.25,&reason),"Resolved shock heat was rejected");
        require(close((p.gamma-1)*dest[id][EINT],11),"Resolved shock heating was erased");
        v[B1]=200;q=conserved(v,1,p);
        require(!modern_energy_reliable(q,1,p),"Extreme low-beta limitation changed without validation");
        auto sig=pymhd_checkpoint::signature(p);p.internal_energy_sync=0;
        require(sig!=pymhd_checkpoint::signature(p),"Restart did not record synchronization policy");
        p.internal_energy_sync=1;p.dual_energy_eta*=2;
        require(sig!=pymhd_checkpoint::signature(p),"Restart did not record modern eta");
        std::cout<<"PASS modern dual-energy three gates, neighbour-veto regression, invalid channels, shock heat and exact correction accounting\n"
                 <<"LIMITATION: modern eta=.001 does not guarantee correct shock heating when downstream beta < .0012 (gamma=5/3)\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL modern dual energy: "<<e.what()<<'\n';return 1;}
}
