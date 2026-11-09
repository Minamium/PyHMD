#pragma once
#include "pyhmd/geometry/fields.hpp"
#include "pyhmd/problems/benchmarks.hpp"

HD inline Primitive initial(double r,double th,Parameters p) {
    if(benchmark_problem(p.problem))return benchmark_primitive(r,th,0,p);
    Primitive v; v[RHO]=1;v[ENE]=1;
    if(p.problem==ALFVEN) {v[B1]=1;v[B2]=p.wave_amplitude*sin(2*PI*r);v[B3]=p.wave_amplitude*cos(2*PI*r);v[M2]=-v[B2];v[M3]=-v[B3];}
    if(p.problem==DIFFUSION)v[B3]=p.wave_amplitude*sin(2*PI*r);
    if(p.problem==BRIO_WU){bool l=r<.5;v[RHO]=l?1:.125;v[ENE]=l?1:.1;v[B1]=.75;v[B2]=l?1:-1;}
    if(p.problem==DISK) {
        // The paper-specific polytropic equilibrium is filled by the host initializer.
        double rr=fmax(r*sin(th),1e-12),z=r*cos(th);
        double f=fmax(0.,1.-(p.gamma-1)/(p.gamma*p.h*p.h)*(1.-sin(th)));
        double rho_d=pow(rr,-1.5)*pow(f,1/(p.gamma-1));
        double rho_c=p.corona*pow(r,-1.5),pc=(p.gamma-1)/p.gamma*p.corona*pow(r,-2.5);
        double pd=p.h*p.h*pow(rr,-2.5)*pow(f,p.gamma/(p.gamma-1));
        if(pd>=pc && rho_d>rho_c) {v[RHO]=rho_d;v[ENE]=pd;v[M3]=sqrt(fmax(0.,1.-2.5*p.h*p.h))/sqrt(rr);}
        else {v[RHO]=rho_c;v[ENE]=pc;v[M3]=0;}
        (void)z;
        v[RHO]=fmax(v[RHO],p.rho_floor);v[ENE]=fmax(v[ENE],p.p_floor);
    }
    return v;
}
HD inline double initial_potential(double x,double y,Parameters p) {
    if(!p.spherical) {
        if(benchmark_problem(p.problem))return benchmark_potential(x,y,0,p);
        if(p.problem==ALFVEN)return y+p.wave_amplitude/(2*PI)*cos(2*PI*x);
        if(p.problem==BRIO_WU)return .75*y-(x<=.5?x:1.-x);
        return 0.;
    }
    if(p.problem!=DISK)return 0.;
    double rc=fmax(fabs(x*sin(y)),1e-30), ct=cos(y)/fmax(fabs(sin(y)),1e-30);
    double ap=(4./3)*p.b0*pow(rc,-.25)*pow(p.bend,1.25)/pow(p.bend*p.bend+ct*ct,.625);
    return rc*ap;
}
