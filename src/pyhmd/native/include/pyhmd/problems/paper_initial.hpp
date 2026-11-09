#pragma once
#include "pyhmd/problems/builtin.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// Host-only, self-similar magnetohydrostatic initial disk.  See
// docs/numerical_model.md: SF2014's prose does not uniquely specify whether its
// original initializer included the Lorentz force.  This implements the force
// equilibrium stated in its Eq. (7), using the Zanni et al. (2007) construction.
// CT, independently, initializes B from the common analytic flux function.
class PaperInitial {
    Parameters p_;
    double dx_, p0_, k_;
    std::vector<double> density_;

    struct FieldShape { double f, fp, j; };
    FieldShape shape(double x) const {
        const double a=p_.bend*p_.bend+x*x;
        const double f=std::pow(p_.bend,1.25)*std::pow(a,-.625);
        const double fp=-1.25*x*f/a;
        const double fpp=f*(-1.25/a+4.0625*x*x/(a*a));
        return {f,fp,.9375*f-1.5*x*fp-(1+x*x)*fpp};
    }
    double derivative(double x,double d) const {
        if(!(d>0) || !std::isfinite(d))return NAN;
        const auto b=shape(x);
        const double force=-p_.gm*d*x/std::pow(1+x*x,1.5)
                           +k_*k_*b.j*b.fp;
        return force/(p0_*p_.gamma*std::pow(d,p_.gamma-1));
    }
    double rotation_squared(double x,double d) const {
        const auto b=shape(x);
        return p_.gm/std::sqrt(1+x*x)
             -2.5*p0_*std::pow(d,p_.gamma-1)
             -.75*k_*k_*b.j*b.f/d;
    }
    double corona_pressure_at_unit_cylindrical_radius(double x) const {
        return (p_.gamma-1)/p_.gamma*p_.gm*p_.corona
             *std::pow(1+x*x,-1.25);
    }
public:
    explicit PaperInitial(Parameters p):p_(p) {
        if(!(p_.h>0 && p_.corona>0 && p_.bend>0 && p_.gm>0)
           || !std::isfinite(p_.h) || !std::isfinite(p_.corona)
           || !std::isfinite(p_.bend) || !std::isfinite(p_.gm)
           || !std::isfinite(p_.b0))
            throw std::invalid_argument("PaperInitial requires positive h, corona, bend and GM, and finite B0");
        if(std::abs(p_.gamma-5./3)>1e-12)
            throw std::invalid_argument("PaperInitial self-similar polytrope requires gamma=5/3");
        p0_=p_.h*p_.h*p_.gm;
        // The prose and Eq.(7) of SF2014 do not uniquely identify the original
        // initializer. Keep both choices explicit, with the same CT field.
        k_=p_.initial_equilibrium==0?0.:4./3*p_.b0;
        if(rotation_squared(0,1)<=0)
            throw std::invalid_argument("PaperInitial has no real midplane rotation for these parameters");
        if(p0_<=corona_pressure_at_unit_cylindrical_radius(0))
            throw std::invalid_argument("PaperInitial corona pressure exceeds disk midplane pressure");

        // Fine uniform x=z/s table; all four RK stages evaluate the same
        // analytic Lorentz force. End at the disk/corona pressure crossing.
        const double xmax=std::max(1.,8*p_.h);
        constexpr int n=65536;
        dx_=xmax/n;
        density_.reserve(n+1);
        density_.push_back(1.);
        for(int i=0;i<n;++i) {
            const double x=i*dx_,d=density_.back();
            const double k1=derivative(x,d);
            const double k2=derivative(x+dx_/2,d+dx_*k1/2);
            const double k3=derivative(x+dx_/2,d+dx_*k2/2);
            const double k4=derivative(x+dx_,d+dx_*k3);
            const double next=d+dx_*(k1+2*k2+2*k3+k4)/6;
            if(!(next>0) || !std::isfinite(next))
                throw std::runtime_error("PaperInitial ODE lost positivity before corona match");
            density_.push_back(next);
            const double xn=x+dx_;
            if(p0_*std::pow(next,p_.gamma)
                    <=corona_pressure_at_unit_cylindrical_radius(xn)
                || rotation_squared(xn,next)<=0)return;
        }
        throw std::runtime_error("PaperInitial disk/corona crossing outside ODE table");
    }

    Primitive at(double spherical_r,double theta) const {
        if(!(spherical_r>0) || !std::isfinite(spherical_r) || !std::isfinite(theta))
            throw std::invalid_argument("PaperInitial requires positive radius and finite theta");
        Primitive v;
        v[RHO]=p_.corona*std::pow(spherical_r,-1.5);
        v[ENE]=(p_.gamma-1)/p_.gamma*p_.gm*p_.corona
               *std::pow(spherical_r,-2.5);
        const double st=std::abs(std::sin(theta));
        const double s=spherical_r*st;
        if(st>1e-14) {
            const double x=std::abs(std::cos(theta))/st;
            const double fi=x/dx_;
            if(fi<static_cast<double>(density_.size()-1)) {
                const auto i=static_cast<std::size_t>(fi);
                const double w=fi-i;
                const double d=(1-w)*density_[i]+w*density_[i+1];
                const double pd=p0_*std::pow(s,-2.5)*std::pow(d,p_.gamma);
                const double v2=rotation_squared(x,d)/s;
                if(pd>=v[ENE] && v2>0) {
                    v[RHO]=std::pow(s,-1.5)*d;
                    v[ENE]=pd;
                    v[M3]=std::sqrt(v2);
                }
            }
        }
        // B is deliberately zero here; the host must insert the authoritative
        // CT face-derived cell fields before converting to conserved energy.
        v[RHO]=std::max(v[RHO],p_.rho_floor);
        v[ENE]=std::max(v[ENE],p_.p_floor);
        return v;
    }
};
