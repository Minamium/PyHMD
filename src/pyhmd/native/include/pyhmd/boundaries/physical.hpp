#pragma once
#include "pyhmd/physics/resistivity.hpp"

HD inline void theta_boundary(int k,State* u,double* a,MeshView m,Parameters p) {
    int i=k%(m.tx+1),g=k/(m.tx+1);if(g>=NG)return;
    if(i<m.tx)for(int upper=0;upper<2;++upper) {
        int bc=upper?p.upper2:p.lower2,j=upper?NG+m.ny+g:NG-1-g;
        int src=bc==PERIODIC?(upper?NG+g:NG+m.ny-1-g):(bc==OUTFLOW?(upper?NG+m.ny-1:NG):(upper?NG+m.ny-1-g:NG+g));
        Primitive v=get_primitive(u,i,src,m,p);
        if(bc==REFLECT||bc==AXIS){v[M2]*=-1;v[B2]*=-1;if(bc==AXIS){v[M3]*=-1;v[B3]*=-1;}}
        if(bc==EQUATOR){v[M2]*=-1;v[B1]*=-1;v[B3]*=-1;}
        if(bc==OUTFLOW)v[M2]=upper?fmax(v[M2],0.):fmin(v[M2],0.);
        u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
    }
    // Potential vertices: axis/equator reflection, or linear continuation/outflow.
    for(int upper=0;upper<2;++upper) {
        int bc=upper?p.upper2:p.lower2,edge=upper?NG+m.ny:NG,vertex=upper?edge+g+1:edge-g-1;
        int src=bc==PERIODIC?(upper?NG+g+1:NG+m.ny-g-1):(upper?edge-g-1:edge+g+1);
        double val=a[m.vertex(i,src)];
        if(bc==OUTFLOW)val=a[m.vertex(i,edge)]+(a[m.vertex(i,edge)]-a[m.vertex(i,upper?edge-1:edge+1)])*(g+1);
        if(bc==PERIODIC) {
            // A flux function is quasi-periodic when the domain carries net
            // normal magnetic flux. Infer its jump from the physical vertices;
            // hard-coded test-case fields corrupt arbitrary periodic B fields.
            double jump=a[m.vertex(i,NG+m.ny)]-a[m.vertex(i,NG)];
            val+=(upper?1.:-1.)*jump;
        }
        // The axisymmetric flux function is EVEN at the axis (Br even, Btheta odd).
        a[m.vertex(i,vertex)]=val;
        if(bc==AXIS && g==0)a[m.vertex(i,edge)]=0;
    }
}
HD inline void radial_boundary(int k,State* u,MeshView m,Parameters p,int lower,int upper) {
    int j=k/m.tx,i=k%m.tx;if(j>=m.ty)return;
    if((i>=NG&&i<NG+m.nx)||(!lower&&i<NG)||(!upper&&i>=NG+m.nx))return;
    bool is_upper=i>=NG+m.nx;int bc=is_upper?p.upper1:p.lower1;
    int src=bc==PERIODIC?(is_upper?i-m.nx:i+m.nx):(is_upper?NG+m.nx-1:NG);
    Primitive v=get_primitive(u,src,j,m,p);
    if(bc==PAPER_INNER||bc==PAPER_OUTER) {
        double fac=m.xc[i]/m.xc[src];v[RHO]*=pow(fac,-1.5);v[ENE]*=pow(fac,-2.5);v[B3]/=fac;
        bool disk=m.yc[j]>PI/2-2*p.h;
        if(bc==PAPER_INNER) {
            if(disk){v[M1]=fmin(v[M1]*pow(fac,-.5),0.);v[M2]=0;v[M3]*=pow(fac,-.5);}
            else {
                Primitive mid=get_primitive(u,NG,NG+m.ny-1,m,p);
                double angle=(70*PI/180)/(1+exp(-(m.yc[j]-PI/4)/(PI/12)));
                v[RHO]=fmax(p.rho_floor,p.corona*mid[RHO]*pow(m.xc[i]/m.xc[NG],-1.5));
                v[ENE]=(p.gamma-1)/p.gamma*p.gm*v[RHO]/m.xc[i];
                v[M1]=.2*cos(angle);v[M2]=-.2*sin(angle);v[M3]=0;v[B3]=0;
            }
        } else v[M1]=disk?fmin(v[M1],0.):fmax(v[M1],0.);
        // Density/pressure are prescribed anew at the radial boundary.
        v[ENTR]=v[ENE]*pow(v[RHO],-p.gamma);
        v[EINT]=v[ENE]/(p.gamma-1);
    }
    if(bc==OUTFLOW)v[M1]=is_upper?fmax(v[M1],0.):fmin(v[M1],0.);
    u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
}
HD inline void potential_radial_boundary(int j,double* a,MeshView m,Parameters p,int lower,int upper) {
    if(j>m.ty)return;
    if(lower&&p.lower1!=PAPER_INNER)for(int g=1;g<=NG;++g) {
        int i=NG-g;
        if(p.lower1==PERIODIC)a[m.vertex(i,j)]=a[m.vertex(i+m.nx,j)]-(a[m.vertex(NG+m.nx,j)]-a[m.vertex(NG,j)]);
        else a[m.vertex(i,j)]=a[m.vertex(NG,j)]+(a[m.vertex(NG+1,j)]-a[m.vertex(NG,j)])*(m.xf[i]-m.xf[NG])/(m.xf[NG+1]-m.xf[NG]);
    }
    if(upper)for(int g=1;g<=NG;++g) {
        int edge=NG+m.nx,i=edge+g;
        if(p.upper1==PERIODIC)a[m.vertex(i,j)]=a[m.vertex(i-m.nx,j)]+(a[m.vertex(edge,j)]-a[m.vertex(NG,j)]);
        else {
            double xi=p.spherical?sq(m.xf[i]):m.xf[i],xe=p.spherical?sq(m.xf[edge]):m.xf[edge],xl=p.spherical?sq(m.xf[edge-1]):m.xf[edge-1];
            a[m.vertex(i,j)]=a[m.vertex(edge,j)]+(a[m.vertex(edge,j)]-a[m.vertex(edge-1,j)])*(xi-xe)/(xe-xl);
        }
    }
}
HD inline void paper_inner_potential(int unused,double* a,MeshView m,Parameters p) {
    (void)unused;
    // y'=k(theta)(source-y). Freeze k at each angular interval midpoint and
    // integrate a linearly varying source exactly. The nonnegative exponential
    // weights avoid the stiff, sign-alternating trapezoidal boundary mode.
    // The physical normal field at i=NG is a solution variable evolved by
    // Faraday's law. A ghost boundary must not reproject that active face.
    for(int i=p.inner_flux_project?NG:NG-1;i>=0;--i) {
        a[m.vertex(i,NG)]=0;
        double prev=0,source_prev=a[m.vertex(i+1,NG)];
        for(int j=NG+1;j<=NG+m.ny;++j) {
            double t=.5*(m.yf[j]+m.yf[j-1]),ang=(70*PI/180)/(1+exp(-(t-PI/4)/(PI/12)));
            double kval=m.xf[i]/((m.xf[i+1]-m.xf[i])*fmax(tan(ang),1e-10));
            double src=a[m.vertex(i+1,j)],dy=m.yf[j]-m.yf[j-1];
            double z=kval*dy,decay=exp(-z);
            double avg=z<1e-5?1-z/2+z*z/6-z*z*z/24:-expm1(-z)/z;
            double value=decay*prev+(avg-decay)*source_prev+(1-avg)*src;
            a[m.vertex(i,j)]=value;prev=value;source_prev=src;
        }
    }
}
