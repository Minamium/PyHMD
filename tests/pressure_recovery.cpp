// Postcondition tests for representable gas pressure at disparate energy scales.
// Build: c++ -O2 -std=c++17 -Iinclude tests/pressure_recovery.cpp -o build/pressure_recovery
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>
#include <random>

static void require(bool condition,const char* message) {
    if(!condition)throw std::runtime_error(message);
}

static void check_case(double rho,double kinetic,double magnetic,double radius,double initial_energy,double pressure) {
    Parameters p;
    State q;q[RHO]=rho;q[M1]=std::sqrt(2*rho*kinetic/3);
    q[M2]=q[M1];q[M3]=radius*q[M1];q[B1]=std::sqrt(2*magnetic);
    q[ENE]=initial_energy;q[ENTR]=pressure*std::pow(rho,1-p.gamma);
    require(set_pressure(q,radius,pressure,p),"Finite pressure reconstruction failed");
    const double recovered=raw_pressure(q,radius,p);
    require(finite_value(q[ENE])&&finite_value(recovered)&&recovered>0,"Recovery reported success with invalid raw pressure");
    require(recovered>=pressure*(1-1e-12),"Recovered pressure lies below the requested representable lower bound");
    // Recovery must not alter mass, momentum or magnetic fields.
    require(q[RHO]==rho&&q[M1]==q[M2]&&q[M3]==radius*q[M1],"Pressure recovery changed mass or momentum");
    require(q[B1]==std::sqrt(2*magnetic),"Pressure recovery changed the magnetic field");
}

int main() {
    try {
        // The old E += deficit, then single-nextafter implementation fails this
        // finite counterexample because subtracting -1e20 loses the O(1) part.
        check_case(1,1,0,1,-1e20,1e-12);
        check_case(1,1e20,1e20,.3,-1e20,1e-12);
        check_case(1e-8,1e-8,1e-12,1e-3,-1e-7,1e-12);
        check_case(1,.01,.1,.4,.1100000000001,1e-12);
        std::mt19937_64 generator(20261105);
        std::uniform_real_distribution<double> exponent(-15,15),unit(.1,1);
        for(int n=0;n<10000;++n) {
            double rho=std::pow(10.,exponent(generator));
            double kinetic=std::pow(10.,exponent(generator)),magnetic=std::pow(10.,exponent(generator));
            double radius=unit(generator),pressure=std::pow(10.,exponent(generator));
            check_case(rho,kinetic,magnetic,radius,-std::pow(10.,exponent(generator)),pressure);
        }
        Parameters p;State invalid;invalid[RHO]=1;invalid[M1]=INFINITY;
        require(!set_pressure(invalid,1,p.p_floor,p),"Nonfinite kinetic state was accepted");
        invalid=State{};invalid[RHO]=-1;
        require(!set_pressure(invalid,1,p.p_floor,p),"Negative density was accepted");
        invalid=State{};invalid[RHO]=1;
        require(!set_pressure(invalid,1,NAN,p),"Nonfinite requested pressure was accepted");
        std::cout<<"PASS pressure reconstruction: 10004 finite cases and 3 invalid-input checks\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL pressure reconstruction: "<<error.what()<<'\n';return 1;}
}
