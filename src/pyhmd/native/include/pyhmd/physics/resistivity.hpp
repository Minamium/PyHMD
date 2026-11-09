#pragma once
#include "pyhmd/numerics/reconstruction.hpp"

HD inline double eta_at(const State* u,int i,double radius,double theta,MeshView m,Parameters p) {
    if(p.eta_model==0)return p.eta;
    double h=2*p.h*fmax(radius*fabs(sin(theta)),1e-30),z=fabs(radius*cos(theta));
    double profile=z<=h?1.:exp(-2*sq((z-h)/h));
    if(profile<1e-30)return 0.;
    int jj=NG+m.ny-1;
    Primitive v=get_primitive(u,i,jj,m,p);
    double va=sqrt((sq(v[B1])+sq(v[B2]))/v[RHO]),cs=sqrt(v[ENE]/v[RHO]),mu=(sq(v[B1])+sq(v[B2]))/(2*v[ENE]);
    if(p.midplane) {
        // Lookup at cylindrical s=r sin(theta), including cross-rank radii.
        double s=radius*fabs(sin(theta));
        double q=p.log_r?log(s/p.xlo)/log(p.xhi/p.xlo)*p.nx:(s-p.xlo)/(p.xhi-p.xlo)*p.nx;
        int index=static_cast<int>(floor(q-.5));index=index<0?0:(index>p.nx-2?p.nx-2:index);
        double faces[3];for(int d=0;d<3;++d)faces[d]=p.log_r?p.xlo*exp(log(p.xhi/p.xlo)*(index+d)/p.nx):p.xlo+(p.xhi-p.xlo)*(index+d)/p.nx;
        double x0=p.spherical?cbrt(.5*(pow(faces[0],3)+pow(faces[1],3))):.5*(faces[0]+faces[1]);
        double x1=p.spherical?cbrt(.5*(pow(faces[1],3)+pow(faces[2],3))):.5*(faces[1]+faces[2]);
        double f=(s-x0)/(x1-x0);Midplane l=p.midplane[index],r=p.midplane[index+1];
        va=fmax(0.,l.va+f*(r.va-l.va));cs=fmax(0.,l.cs+f*(r.cs-l.cs));mu=fmax(0.,l.mu+f*(r.mu-l.mu));
    }
    double alpha=p.alpha_m;
    if(p.eta_model==2) {
        // Sec.5 strong diffusion: alpha_ss = alpha_m sqrt(2 mu_ref) (mu/mu_ref)^2.
        return alpha*sqrt(2*p.mu_ref)*sq(mu/p.mu_ref)*cs*h*profile;
    }
    return alpha*va*h*profile;
}
HD inline EtaGeometry make_eta_geometry(double radius,double theta,Parameters p) {
    EtaGeometry g;g.height=2*p.h*fmax(radius*fabs(sin(theta)),1e-30);
    double z=fabs(radius*cos(theta));g.profile=z<=g.height?1.:exp(-2*sq((z-g.height)/g.height));
    if(g.profile<1e-30){g.profile=0;return g;}
    double s=radius*fabs(sin(theta));
    double q=p.log_r?log(s/p.xlo)/log(p.xhi/p.xlo)*p.nx:(s-p.xlo)/(p.xhi-p.xlo)*p.nx;
    int index=static_cast<int>(floor(q-.5));index=index<0?0:(index>p.nx-2?p.nx-2:index);
    double faces[3];for(int d=0;d<3;++d)faces[d]=p.log_r?p.xlo*exp(log(p.xhi/p.xlo)*(index+d)/p.nx):p.xlo+(p.xhi-p.xlo)*(index+d)/p.nx;
    double x0=p.spherical?cbrt(.5*(pow(faces[0],3)+pow(faces[1],3))):.5*(faces[0]+faces[1]);
    double x1=p.spherical?cbrt(.5*(pow(faces[1],3)+pow(faces[2],3))):.5*(faces[1]+faces[2]);
    g.index=index;g.fraction=(s-x0)/(x1-x0);return g;
}
HD inline double eta_from_geometry(EtaGeometry g,Parameters p) {
    // Cached geometry can outlive a runtime operator override used by
    // diagnostics. A constant/ideal model must not inherit SF2014 eta.
    if(p.eta_model==0)return p.eta;
    if(g.profile==0)return 0;
    Midplane l=p.midplane[g.index],r=p.midplane[g.index+1];
    if(p.eta_model==2) {
        double cs=fmax(0.,l.cs+g.fraction*(r.cs-l.cs)),mu=fmax(0.,l.mu+g.fraction*(r.mu-l.mu));
        return p.alpha_m*sqrt(2*p.mu_ref)*sq(mu/p.mu_ref)*cs*g.height*g.profile;
    }
    double va=fmax(0.,l.va+g.fraction*(r.va-l.va));
    return p.alpha_m*va*g.height*g.profile;
}
HD inline double eta_cell_value(const State* u,int i,int j,MeshView m,Parameters p) {
    if(p.reconstruction_cache&&p.midplane&&p.eta_cell_geometry)return eta_from_geometry(p.eta_cell_geometry[m.cell(i,j)],p);
    return eta_at(u,i,m.xc[i],m.yc[j],m,p);
}
HD inline double current_phi(const double* a,int i,int j,MeshView m,Parameters p) {
    double dr=m.xc[i]-m.xc[i-1],dy=m.yc[j]-m.yc[j-1];
    if(p.spherical)return ((m.xc[i]*field2(a,i,j,m,p)-m.xc[i-1]*field2(a,i-1,j,m,p))/dr
        -(field1(a,i,j,m,p)-field1(a,i,j-1,m,p))/dy)/m.xf[i];
    return (field2(a,i,j,m,p)-field2(a,i-1,j,m,p))/dr-(field1(a,i,j,m,p)-field1(a,i,j-1,m,p))/dy;
}
HD inline Current cell_current(const State* u,const double* a,int i,int j,MeshView m,Parameters p) {
    if(p.current_cache)return p.current_cache[m.cell(i,j)];
    Current c;double r=m.xc[i];
    double bleft=u[m.cell(i-1,j)][B3],bright=u[m.cell(i+1,j)][B3];
    double bdown=u[m.cell(i,j-1)][B3],bup=u[m.cell(i,j+1)][B3];
    if(p.spherical) {
        c[0]=(m.sine_center(j+1)*bup-m.sine_center(j-1)*bdown)/(r*m.sine_center(j)*(m.yc[j+1]-m.yc[j-1]));
        c[1]=-(m.xc[i+1]*bright-m.xc[i-1]*bleft)/(r*(m.xc[i+1]-m.xc[i-1]));
    } else {c[0]=(bup-bdown)/(m.yc[j+1]-m.yc[j-1]);c[1]=-(bright-bleft)/(m.xc[i+1]-m.xc[i-1]);}
    c[2]=.25*(current_phi(a,i,j,m,p)+current_phi(a,i+1,j,m,p)+current_phi(a,i,j+1,m,p)+current_phi(a,i+1,j+1,m,p));
    return c;
}
