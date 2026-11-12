#pragma once
#include "pymhd/core/types.hpp"

HD inline double sq(double a){return a*a;}
HD inline double clamp(double a,double lo,double hi){return fmax(lo,fmin(hi,a));}
HD inline double potential(double r,Parameters p){return p.spherical?-p.gm/r:0.;}
HD inline Primitive primitive(State u,double radius,Parameters p) {
    Primitive v; v[RHO]=fmax(u[RHO],p.rho_floor);
    double safe_radius=fabs(radius)>1e-30?radius:(radius<0?-1e-30:1e-30);
    v[M1]=u[M1]/v[RHO]; v[M2]=u[M2]/v[RHO]; v[M3]=u[M3]/(v[RHO]*safe_radius);
    for(int k=B1;k<=B3;++k)v[k]=u[k];
    v[ENE]=fmax((p.gamma-1)*(u[ENE]-.5*v[RHO]*(sq(v[M1])+sq(v[M2])+sq(v[M3]))-.5*(sq(v[B1])+sq(v[B2])+sq(v[B3]))),p.p_floor);
    v[ENTR]=u[ENTR]/v[RHO];
    v[EINT]=u[EINT];
    if(p.entropy_switch==3&&u[EINT]>0&&u[EINT]<INFINITY)v[ENE]=fmax((p.gamma-1)*u[EINT],p.p_floor);
    return v;
}
HD inline State conserved(Primitive v,double radius,Parameters p) {
    State u;u[RHO]=v[RHO];u[M1]=v[RHO]*v[M1];u[M2]=v[RHO]*v[M2];u[M3]=radius*v[RHO]*v[M3];
    for(int k=B1;k<=B3;++k)u[k]=v[k];
    u[ENE]=v[ENE]/(p.gamma-1)+.5*v[RHO]*(sq(v[M1])+sq(v[M2])+sq(v[M3]))+.5*(sq(v[B1])+sq(v[B2])+sq(v[B3]));
    u[ENTR]=v[ENTR]>0?v[RHO]*v[ENTR]:v[ENE]*pow(v[RHO],1-p.gamma);
    u[EINT]=v[EINT]>0?v[EINT]:v[ENE]/(p.gamma-1);
    return u;
}
HD inline double fast_speed(Primitive v,int dir,Parameters p) {
    double a2=p.gamma*v[ENE]/v[RHO], b2=(sq(v[B1])+sq(v[B2])+sq(v[B3]))/v[RHO],bn2=sq(v[B1+dir])/v[RHO];
    return sqrt(.5*(a2+b2+sqrt(fmax(0.,sq(a2+b2)-4*a2*bn2))));
}
HD inline State ideal_flux(Primitive v,int dir,Parameters p) {
    State u=conserved(v,1,p),f;double vn=v[M1+dir],bn=v[B1+dir],pt=v[ENE]+.5*(sq(v[B1])+sq(v[B2])+sq(v[B3]));
    f[RHO]=v[RHO]*vn;
    for(int a=0;a<3;++a) {f[M1+a]=v[RHO]*vn*v[M1+a]-bn*v[B1+a]+(a==dir?pt:0.);f[B1+a]=vn*v[B1+a]-bn*v[M1+a];}
    f[ENE]=(u[ENE]+pt)*vn-bn*(v[M1]*v[B1]+v[M2]*v[B2]+v[M3]*v[B3]);
    f[ENTR]=u[ENTR]*vn;f[EINT]=u[EINT]*vn;return f;
}
HD inline State hll(Primitive l,Primitive r,int dir,double bn,Parameters p,double* interface_velocity=nullptr,FaceSignal* signal=nullptr,HllTransport* transport=nullptr) {
    // CT owns the unique normal face field; rebuild energies with unchanged gas pressure.
    l[B1+dir]=r[B1+dir]=bn;
    l[EINT]=l[ENE]/(p.gamma-1);r[EINT]=r[ENE]/(p.gamma-1);
    State fl=ideal_flux(l,dir,p),fr=ideal_flux(r,dir,p),ul=conserved(l,1,p),ur=conserved(r,1,p),f;
    double sl=fmin(l[M1+dir]-fast_speed(l,dir,p),r[M1+dir]-fast_speed(r,dir,p));
    double sr=fmax(l[M1+dir]+fast_speed(l,dir,p),r[M1+dir]+fast_speed(r,dir,p));
    if(signal) {
        double ap=fmax(0.,sr),am=fmax(0.,-sl),den=fmax(ap+am,1e-30);
        signal->left=ap/den;signal->right=am/den;signal->diffusion=ap*am/den;
        signal->transverse_velocity=signal->left*l[M1+1-dir]+signal->right*r[M1+1-dir];
    }
    for(int k=0;k<NV;++k)f[k]=sl>=0?fl[k]:(sr<=0?fr[k]:(sr*fl[k]-sl*fr[k]+sl*sr*(ur[k]-ul[k]))/(sr-sl));
    if(transport) {
        double wl=sl>=0?1.:(sr<=0?0.:sr/(sr-sl));
        double wr=sl>=0?0.:(sr<=0?1.:-sl/(sr-sl));
        double jump=sl<0&&sr>0?sl*sr/(sr-sl):0.;
        transport->mass_adv=wl*fl[RHO]+wr*fr[RHO];
        transport->mass_jump=jump*(ur[RHO]-ul[RHO]);
        transport->angular_adv=wl*l[RHO]*l[M1+dir]*l[M3]+wr*r[RHO]*r[M1+dir]*r[M3];
        transport->angular_mag=-bn*(wl*l[B3]+wr*r[B3]);
        transport->angular_jump=jump*(ur[M3]-ul[M3]);
    }
    if(interface_velocity) {
        double density=sr*ur[RHO]-sl*ul[RHO]-(fr[RHO]-fl[RHO]);
        double momentum=sr*ur[M1+dir]-sl*ul[M1+dir]-(fr[M1+dir]-fl[M1+dir]);
        *interface_velocity=sl>=0?l[M1+dir]:(sr<=0?r[M1+dir]:(density>0?momentum/density:.5*(l[M1+dir]+r[M1+dir])));
    }
    f[B1+dir]=0;return f;
}
