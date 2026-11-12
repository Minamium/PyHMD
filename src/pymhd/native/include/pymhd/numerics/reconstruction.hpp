#pragma once
#include "pymhd/geometry/fields.hpp"
#include "pymhd/numerics/metric_reconstruction.hpp"

HD inline bool finite_value(double x) {
#ifdef __CUDA_ARCH__
    return isfinite(x);
#else
    return std::isfinite(x);
#endif
}

// Each operation is shared by the C++ reference backend and the CUDA backend.
HD inline Primitive get_primitive(const State* u,int i,int j,MeshView m,Parameters p) {
    if(p.primitive_cache)return p.primitive_cache[m.cell(i,j)];
    return primitive(u[m.cell(i,j)],m.radius(i,j,p),p);
}
HD inline void limit_ppm_edges(double qc,double& left,double& right,int i,int j,int dir,MeshView m,Parameters p,bool staggered=false) {
    const PpmGeometry* geometry=dir==0?(staggered?m.gx_staggered:m.gx):m.gy;
    if(p.metric_ppm&&geometry){metric_ppm_limit(qc,left,right,geometry[dir==0?i:j]);return;}
    double dq=right-left,q6=6*(qc-.5*(left+right));
    if((right-qc)*(qc-left)<=0)left=right=qc;
    else if(dq*q6>dq*dq)left=3*qc-2*right;
    else if(dq*q6<-dq*dq)right=3*qc-2*left;
}
HD inline double ppm_edge(const State* u,int i,int j,int dir,int side,int variable,MeshView m,Parameters p) {
    int c=dir==0?i:j; const Weights* weights=dir==0?m.wx:m.wy;
    double edge[2];
    for(int s=0;s<2;++s) {
        int face=c+s;double q=0;
        for(int a=0;a<4;++a) {
            int offset=face-2+a;
            q+=weights[face].w[a]*get_primitive(u,dir==0?offset:i,dir==1?offset:j,m,p)[variable];
        }
        double q0=get_primitive(u,dir==0?face-1:i,dir==1?face-1:j,m,p)[variable];
        double q1=get_primitive(u,dir==0?face:i,dir==1?face:j,m,p)[variable];
        edge[s]=clamp(q,fmin(q0,q1),fmax(q0,q1));
    }
    double qc=get_primitive(u,i,j,m,p)[variable];limit_ppm_edges(qc,edge[0],edge[1],i,j,dir,m,p);
    return edge[side];
}
HD inline void reconstructed_pair(const State* u,int i,int j,int dir,MeshView m,Parameters p,Primitive& left,Primitive& right) {
    if(p.flat_reconstruction){left=right=get_primitive(u,i,j,m,p);return;}
    // Both edges use the same five-cell stencil. Convert each state once, not
    // separately for every primitive component and each repeated face sample.
    Primitive stencil[5];
    for(int a=0;a<5;++a)stencil[a]=get_primitive(u,i+(dir==0?a-2:0),j+(dir==1?a-2:0),m,p);
    int c=dir==0?i:j;const Weights* weights=dir==0?m.wx:m.wy;
    for(int k=0;k<NV;++k) {
        double edge[2];
        for(int s=0;s<2;++s) {
            double value=0;
            for(int a=0;a<4;++a)value+=weights[c+s].w[a]*stencil[a+s][k];
            edge[s]=clamp(value,fmin(stencil[1+s][k],stencil[2+s][k]),fmax(stencil[1+s][k],stencil[2+s][k]));
        }
        double qc=stencil[2][k];limit_ppm_edges(qc,edge[0],edge[1],i,j,dir,m,p);
        left[k]=edge[0];right[k]=edge[1];
    }
    left[RHO]=fmax(left[RHO],p.rho_floor);left[ENE]=fmax(left[ENE],p.p_floor);
    right[RHO]=fmax(right[RHO],p.rho_floor);right[ENE]=fmax(right[ENE],p.p_floor);
}
HD inline Primitive reconstructed(const State* u,int i,int j,int dir,int side,MeshView m,Parameters p) {
    const Primitive* cache=dir==0?(side==0?p.edge_xl:p.edge_xr):(side==0?p.edge_yl:p.edge_yr);
    if(cache)return cache[m.cell(i,j)];
    Primitive left,right;reconstructed_pair(u,i,j,dir,m,p,left,right);return side==0?left:right;
}
// GPU PPM preparation: one thread owns one primitive component, retaining only
// five scalar samples. Adjacent threads access adjacent AoS components.
HD inline void reconstructed_component_pair(const Primitive* cells,int i,int j,int dir,int variable,MeshView m,Parameters p,double& left,double& right) {
    if(p.flat_reconstruction){left=right=cells[m.cell(i,j)][variable];return;}
    double stencil[5];
    for(int a=0;a<5;++a)stencil[a]=cells[m.cell(i+(dir==0?a-2:0),j+(dir==1?a-2:0))][variable];
    int c=dir==0?i:j;const Weights* weights=dir==0?m.wx:m.wy;double edge[2];
    for(int s=0;s<2;++s) {
        double value=0;for(int a=0;a<4;++a)value+=weights[c+s].w[a]*stencil[a+s];
        edge[s]=clamp(value,fmin(stencil[1+s],stencil[2+s]),fmax(stencil[1+s],stencil[2+s]));
    }
    double qc=stencil[2];limit_ppm_edges(qc,edge[0],edge[1],i,j,dir,m,p);
    if(variable==RHO){edge[0]=fmax(edge[0],p.rho_floor);edge[1]=fmax(edge[1],p.rho_floor);}
    if(variable==ENE){edge[0]=fmax(edge[0],p.p_floor);edge[1]=fmax(edge[1],p.p_floor);}
    left=edge[0];right=edge[1];
}
HD inline Primitive corner_state(const State* u,int i,int j,int sx,int sy,MeshView m,Parameters p) {
    Primitive c=get_primitive(u,i,j,m,p),x=reconstructed(u,i,j,0,sx,m,p),y=reconstructed(u,i,j,1,sy,m,p),v;
    for(int k=0;k<NV;++k)v[k]=x[k]+y[k]-c[k];
    // A diagonal corner can overshoot even after separate 1D limiting.
    if(v[RHO]<=p.rho_floor || v[ENE]<=p.p_floor)v=c;
    if(v[ENTR]<=0)v[ENTR]=c[ENTR];
    if(v[EINT]<=0)v[EINT]=c[EINT];
    return v;
}
HD inline double staggered_edge(const double* a,int i,int j,int dir,int side,MeshView m,Parameters p) {
    if(p.flat_reconstruction)return dir==0?field2(a,i,j,m,p):field1(a,i,j,m,p);
    int c=dir==0?i:j;const Weights* weights=dir==0?(p.metric_ppm&&m.wx_staggered?m.wx_staggered:m.wx):m.wy;double edges[2];
    for(int s=0;s<2;++s){int face=c+s;double q=0;
        for(int v=0;v<4;++v){int loc=face-2+v;double b=dir==0?field2(a,loc,j,m,p):field1(a,i,loc,m,p);q+=weights[face].w[v]*b;}
        double ql=dir==0?field2(a,face-1,j,m,p):field1(a,i,face-1,m,p);
        double qr=dir==0?field2(a,face,j,m,p):field1(a,i,face,m,p);
        edges[s]=clamp(q,fmin(ql,qr),fmax(ql,qr));
    }
    double qc=dir==0?field2(a,i,j,m,p):field1(a,i,j,m,p);limit_ppm_edges(qc,edges[0],edges[1],i,j,dir,m,p,true);
    return edges[side];
}
