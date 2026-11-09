// Independent checks of the paper initializer. No time evolution is exercised.
// Build: c++ -std=c++17 -O2 -Iinclude tests/paper_initial.cpp -o build/test_paper_initial
#include "pymhd/paper_initial.hpp"
#include <iomanip>
#include <iostream>
#include <stdexcept>

static void require(bool ok,const char* message) {
    if(!ok)throw std::runtime_error(message);
}

int main() {
    try {
        Parameters p;PaperInitial disk(p);
        auto at=[&](double s,double z){return disk.at(std::hypot(s,z),std::atan2(s,z));};
        // Differentiate the scalar flux directly, independently of the
        // initializer's closed-form f', f'' and current formulas.
        auto flux=[&](double s,double z) {
            return (4./3)*p.b0*std::pow(s,.75)*std::pow(p.bend,1.25)
                  /std::pow(p.bend*p.bend+(z/s)*(z/s),.625);
        };
        auto bs=[&](double s,double z) {
            double h=1e-5*s;
            return -(flux(s,z+h)-flux(s,z-h))/(2*h*s);
        };
        auto bz=[&](double s,double z) {
            double h=1e-5*s;
            return (flux(s+h,z)-flux(s-h,z))/(2*h*s);
        };
        double max_residual=0;
        for(double s:{1.,3.,20.})for(double x:{0.,.03,.07,.10,.14,.16,.18}) {
            double z=x*s,h=1e-4*s,r=std::hypot(s,z);auto v=at(s,z);
            require(v[M3]>0,"force-balance sample lies outside rotating disk");
            double j=(bs(s,z+h)-bs(s,z-h))/(2*h)
                    -(bz(s+h,z)-bz(s-h,z))/(2*h);
            double dpz=(at(s,z+h)[ENE]-at(s,z-h)[ENE])/(2*h);
            double dps=(at(s+h,z)[ENE]-at(s-h,z)[ENE])/(2*h);
            double rz=dpz+v[RHO]*p.gm*z/std::pow(r,3)+j*bs(s,z);
            double rs=dps+v[RHO]*p.gm*s/std::pow(r,3)-j*bz(s,z)
                      -v[RHO]*v[M3]*v[M3]/s;
            double force_unit=p.gm*std::pow(s,-3.5);
            max_residual=std::max(max_residual,std::max(std::abs(rs),std::abs(rz))/force_unit);
        }
        require(max_residual<2e-6,"magnetohydrostatic force residual too large");
        auto mid=at(1,0);
        double expected_v2=p.gm*(1-2.5*p.h*p.h)
            -p.b0*p.b0*(1.25+5./(3*p.bend*p.bend));
        require(std::abs(mid[RHO]-1)<1e-12,"incorrect midplane density");
        require(std::abs(mid[ENE]-p.h*p.h*p.gm)<1e-12,"incorrect midplane pressure");
        require(std::abs(mid[M3]*mid[M3]-expected_v2)<1e-12,"incorrect midplane rotation");
        auto axis=disk.at(2,0);
        require(std::isfinite(axis[RHO]) && std::isfinite(axis[ENE]) && axis[M3]==0,
                "invalid axis corona");
        require(std::abs(axis[RHO]-p.corona*std::pow(2.,-1.5))<1e-14,
                "incorrect coronal density");

        // The zero-field limit has a closed-form enthalpy profile. Also vary
        // GM to check the dimensional normalization independently of the
        // reference setting GM=1.
        Parameters hp=p;hp.b0=0;hp.gm=4;PaperInitial hydro(hp);
        for(double x:{0.,.03,.1,.15}) {
            auto v=hydro.at(std::hypot(1.,x),std::atan2(1.,x));
            double w=1+(.4/(hp.h*hp.h))*(1/std::sqrt(1+x*x)-1);
            require(std::abs(v[RHO]-std::pow(w,1.5))<1e-7,"incorrect hydrostatic limit");
            require(std::abs(v[ENE]-hp.h*hp.h*hp.gm*std::pow(v[RHO],hp.gamma))<1e-12,
                    "incorrect polytropic relation");
            require(std::abs(v[M3]*v[M3]-hp.gm*(1-2.5*hp.h*hp.h))<2e-8,
                    "incorrect zero-field rotation");
        }
        // Explicit hydrostatic choice must retain this analytic gas profile
        // even when the separately initialized CT magnetic field is nonzero.
        hp.initial_equilibrium=0;hp.b0=.03;PaperInitial hydro_with_field(hp);
        for(double x:{0.,.03,.1,.15}) {
            const double r=std::hypot(1.,x),theta=std::atan2(1.,x);
            auto a=hydro.at(r,theta),b=hydro_with_field.at(r,theta);
            for(int k:{RHO,ENE,M3})require(a[k]==b[k],"hydrostatic gas profile depends on separately imposed B0");
        }
        bool rejected=false;try {Parameters bad=p;bad.gamma=1.4;PaperInitial invalid(bad);}
        catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"unsupported gamma was accepted");
        std::cout<<std::setprecision(10)<<"PASS paper_initial: max normalized force residual="
                 <<max_residual<<"; midplane, corona, hydro limit, GM scaling, gamma validation\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL paper_initial: "<<e.what()<<'\n';return 1;}
}
