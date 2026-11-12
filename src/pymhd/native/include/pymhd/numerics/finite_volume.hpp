#pragma once
#include "pymhd/physics/resistivity.hpp"

HD inline void update_cell_fields(int k,State* u,const double* a,MeshView m,Parameters p,int correct_energy,double* boundary_energy=nullptr,double weight=0) {
    int i=k%m.tx,j=k/m.tx;
    State v=u[k];double b1=.5*(field1(a,i,j,m,p)+field1(a,i+1,j,m,p));
    double b2=.5*(field2(a,i,j,m,p)+field2(a,i,j+1,m,p));
    if(correct_energy) {
        double delta=.5*(sq(b1)+sq(b2)-sq(v[B1])-sq(v[B2]));v[ENE]+=delta;
        if(boundary_energy&&weight!=0&&i>=NG&&i<NG+m.nx&&j>=NG&&j<NG+m.ny)
            boundary_energy[(j-NG)*m.nx+i-NG]+=weight*delta*m.volume(i,j,p);
    }
    v[B1]=b1;v[B2]=b2;u[k]=v;
}
template<bool Cached> HD inline void face_flux_impl(int k,int dir,const State* u,const double* a,State* flux,MeshView m,Parameters p,State* resistive_flux=nullptr,double* interface_velocity=nullptr,FaceSignal* signals=nullptr,TransportFace* transport=nullptr,double transport_weight=0) {
    int ni=dir==0?m.nx+1:m.nx,nj=dir==0?m.ny:m.ny+1;
    if(k>=ni*nj)return;
    int i=k%ni+NG,j=k/ni+NG;
    int il=i-(dir==0),jl=j-(dir==1);
    Primitive l,r;
    if constexpr(Cached) {
        l=(dir==0?p.edge_xr:p.edge_yr)[m.cell(il,jl)];
        r=(dir==0?p.edge_xl:p.edge_yl)[m.cell(i,j)];
    } else {l=reconstructed(u,il,jl,dir,1,m,p);r=reconstructed(u,i,j,dir,0,m,p);}
    double bn=dir==0?field1(a,i,j,m,p):field2(a,i,j,m,p);
    HllTransport parts;
    double velocity=0;State f=hll(l,r,dir,bn,p,&velocity,signals?signals+m.cell(i,j):nullptr,transport?&parts:nullptr),rf;
    if(interface_velocity)interface_velocity[k]=velocity;
    double radius=dir==0?m.xf[i]:m.xc[i],theta=dir==1?m.yf[j]:m.yc[j];
    const EtaGeometry* geometry=dir==0?p.eta_radial_geometry:p.eta_theta_geometry;
    double ep=p.reconstruction_cache&&p.midplane&&geometry?eta_from_geometry(geometry[k],p):.5*(eta_at(u,il,radius,theta,m,p)+eta_at(u,i,radius,theta,m,p));
    if(ep>0) {
        Current cl,cr,cur;
        if constexpr(Cached) {cl=p.current_cache[m.cell(il,jl)];cr=p.current_cache[m.cell(i,j)];}
        else {cl=cell_current(u,a,il,jl,m,p);cr=cell_current(u,a,i,j,m,p);}
        for(int d=0;d<3;++d)cur[d]=.5*(cl[d]+cr[d]);
        if(dir==0) {
            double dr=m.xc[i]-m.xc[i-1];
            cur[1]=p.spherical?-(m.xc[i]*u[m.cell(i,j)][B3]-m.xc[i-1]*u[m.cell(i-1,j)][B3])/(radius*dr):-(u[m.cell(i,j)][B3]-u[m.cell(i-1,j)][B3])/dr;
            cur[2]=.5*(current_phi(a,i,j,m,p)+current_phi(a,i,j+1,m,p));
        } else {
            double dy=m.yc[j]-m.yc[j-1];
            // Exact axis limit is immaterial to the zero-area flux; avoid division by sin(0).
            cur[0]=p.spherical?(fabs(m.sine_face(j))<1e-20?0.:(m.sine_center(j)*u[m.cell(i,j)][B3]-m.sine_center(j-1)*u[m.cell(i,j-1)][B3])/(radius*m.sine_face(j)*dy)):(u[m.cell(i,j)][B3]-u[m.cell(i,j-1)][B3])/dy;
            cur[2]=.5*(current_phi(a,i,j,m,p)+current_phi(a,i+1,j,m,p));
        }
        double er=ep*p.chi*cur[0],et=ep*p.chi*cur[1],ez=ep*cur[2];
        if(p.eta_model==0){er=ep*cur[0];et=ep*cur[1];}
        // Diffusive Poynting flux pairs with the adjacent-cell arithmetic B,
        // rather than a separately limited hyperbolic reconstruction.
        State ul=u[m.cell(il,jl)],ur=u[m.cell(i,j)];
        double br=dir==0?bn:.5*(ul[B1]+ur[B1]),bt=dir==1?bn:.5*(ul[B2]+ur[B2]),bz=.5*(ul[B3]+ur[B3]);
        if(dir==0){rf[B2]-=ez;rf[B3]+=et;rf[ENE]+=et*bz-ez*bt;}
        else {rf[B1]+=ez;rf[B3]-=er;rf[ENE]+=ez*br-er*bz;}
    }
    for(int c=0;c<NV;++c)f[c]+=rf[c];
    double face_radius=p.spherical?radius*(dir==0?m.sine_center(j):m.sine_face(j)):1.;f[M3]*=face_radius;
    flux[k]=f;
    if(resistive_flux)resistive_flux[k]=rf;
    if(transport&&transport_weight!=0) {
        double weight=transport_weight*(dir==0?m.area1(i,j,p):m.area2(i,j,p));
        auto& q=transport[k];
        q.q[TF_MASS]+=weight*f[RHO];q.q[TF_MASS_ADV]+=weight*parts.mass_adv;q.q[TF_MASS_JUMP]+=weight*parts.mass_jump;
        q.q[TF_ENERGY]+=weight*(f[ENE]+potential(radius,p)*f[RHO]);
        q.q[TF_ANGULAR]+=weight*f[M3];q.q[TF_ANG_ADV]+=weight*face_radius*parts.angular_adv;
        q.q[TF_ANG_MAG]+=weight*face_radius*parts.angular_mag;q.q[TF_ANG_JUMP]+=weight*face_radius*parts.angular_jump;
    }
}
HD inline void face_flux(int k,int dir,const State* u,const double* a,State* flux,MeshView m,Parameters p,State* resistive_flux=nullptr,double* interface_velocity=nullptr) {
    if(p.primitive_cache&&p.edge_xl&&p.edge_xr&&p.edge_yl&&p.edge_yr&&p.current_cache)
        face_flux_impl<true>(k,dir,u,a,flux,m,p,resistive_flux,interface_velocity);
    else face_flux_impl<false>(k,dir,u,a,flux,m,p,resistive_flux,interface_velocity);
}
template<bool Cached> HD inline void corner_emf_impl(int k,const State* u,const double* a,double* arhs,MeshView m,Parameters p,double* resistive_arhs=nullptr) {
    int i=k%(m.nx+1)+NG,j=k/(m.nx+1)+NG;
    double bw=staggered_edge(a,i-1,j,0,1,m,p),be=staggered_edge(a,i,j,0,0,m,p);
    double bs=staggered_edge(a,i,j-1,1,1,m,p),bn=staggered_edge(a,i,j,1,0,m,p);
    double ap0=0,am0=0,ap1=0,am1=0,e[4];
    // Keep one corner state live at a time. The previous four full Primitive
    // states plus a dynamically indexed array caused register spills on SM90.
    // Cached instantiation cannot enter the expensive uncached PPM path.
    for(int q=0;q<4;++q) {
        int east=q&1,north=q>>1,ci=i-1+east,cj=j-1+north;Primitive v;
        if constexpr(Cached) {
            int id=m.cell(ci,cj);const Primitive& c=p.primitive_cache[id];
            const Primitive& x=(east?p.edge_xl:p.edge_xr)[id];
            const Primitive& y=(north?p.edge_yl:p.edge_yr)[id];
            for(int n=0;n<NV;++n)v[n]=x[n]+y[n]-c[n];
            if(v[RHO]<=p.rho_floor||v[ENE]<=p.p_floor)v=c;
        } else v=corner_state(u,ci,cj,1-east,1-north,m,p);
        v[B1]=north?bn:bs;v[B2]=east?be:bw;
        double c0=fast_speed(v,0,p),c1=fast_speed(v,1,p);
        ap0=fmax(ap0,v[M1]+c0);am0=fmax(am0,c0-v[M1]);
        ap1=fmax(ap1,v[M2]+c1);am1=fmax(am1,c1-v[M2]);
        e[q]=v[M2]*v[B1]-v[M1]*v[B2];
    }
    double s0=fmax(ap0+am0,1e-30),s1=fmax(ap1+am1,1e-30);
    double ef=(ap0*ap1*e[0]+am0*ap1*e[1]+ap0*am1*e[2]+am0*am1*e[3])/(s0*s1);
    ef+=ap0*am0/s0*(be-bw);
    ef-=ap1*am1/s1*(bn-bs);
    double ep=p.reconstruction_cache&&p.midplane&&p.eta_corner_geometry?eta_from_geometry(p.eta_corner_geometry[k],p):.5*(eta_at(u,i-1,m.xf[i],m.yf[j],m,p)+eta_at(u,i,m.xf[i],m.yf[j],m,p));
    double resistive=ep*current_phi(a,i,j,m,p);ef+=resistive;
    double value=-m.edge_length(i,j,p)*ef;
    if(p.lower2==AXIS&&j==NG)value=0;
    arhs[m.vertex(i,j)]=value;
    if(resistive_arhs)resistive_arhs[m.vertex(i,j)]=(p.lower2==AXIS&&j==NG)?0.:-m.edge_length(i,j,p)*resistive;
}
HD inline void corner_emf(int k,const State* u,const double* a,double* arhs,MeshView m,Parameters p,double* resistive_arhs=nullptr) {
    if(p.primitive_cache&&p.edge_xl&&p.edge_xr&&p.edge_yl&&p.edge_yr)
        corner_emf_impl<true>(k,u,a,arhs,m,p,resistive_arhs);
    else corner_emf_impl<false>(k,u,a,arhs,m,p,resistive_arhs);
}
// Face-Riemann UCT follows the HLL a/d decomposition of Mignone & Del Zanna
// (2021). Transverse velocities are wave-weighted primitive velocities, not
// the conservative intermediate momentum/density. No vendor code is copied.
HD inline void ghost_face_signal(int k,int dir,const State* u,const double* a,FaceSignal* out,MeshView m,Parameters p) {
    int ni=dir==0?m.nx+1:2*NG,nj=dir==0?2*NG:m.ny+1;
    if(k>=ni*nj)return;
    int g=dir==0?k/ni:k%ni;
    int i=dir==0?k%ni+NG:(g<NG?g:m.nx+g);
    int j=dir==1?k/ni+NG:(g<NG?g:m.ny+g);
    int il=i-(dir==0),jl=j-(dir==1);
    // Normal reconstruction needs no stencil in the transverse direction.
    // Its outer ghost rows were not prepared by the shared 2D edge cache.
    p.edge_xl=p.edge_xr=p.edge_yl=p.edge_yr=nullptr;
    Primitive l=reconstructed(u,il,jl,dir,1,m,p),r=reconstructed(u,i,j,dir,0,m,p);
    double bn=dir==0?field1(a,i,j,m,p):field2(a,i,j,m,p);
    hll(l,r,dir,bn,p,nullptr,out+m.cell(i,j));
}
HD inline double face_velocity_edge(const FaceSignal* signal,int i,int j,int dir,int side,MeshView m,Parameters p) {
    double qc=signal[m.cell(i,j)].transverse_velocity;if(p.flat_reconstruction)return qc;
    int c=dir==0?i:j;const Weights* weights=dir==0?(p.metric_ppm&&m.wx_staggered?m.wx_staggered:m.wx):m.wy;double edge[2];
    for(int s=0;s<2;++s){int face=c+s;double value=0;
        for(int n=0;n<4;++n){int q=face-2+n;value+=weights[face].w[n]*signal[m.cell(dir==0?q:i,dir==1?q:j)].transverse_velocity;}
        double l=signal[m.cell(dir==0?face-1:i,dir==1?face-1:j)].transverse_velocity;
        double r=signal[m.cell(dir==0?face:i,dir==1?face:j)].transverse_velocity;
        edge[s]=clamp(value,fmin(l,r),fmax(l,r));}
    limit_ppm_edges(qc,edge[0],edge[1],i,j,dir,m,p,true);return edge[side];
}
HD inline void face_based_emf(int k,const State* u,const double* a,const FaceSignal* sx,const FaceSignal* sy,double* arhs,MeshView m,Parameters p,double* resistive_arhs,TransportVertex* transport=nullptr,double transport_weight=0) {
    int i=k%(m.nx+1)+NG,j=k/(m.nx+1)+NG;
    FaceSignal xs=sx[m.cell(i,j-1)],xn=sx[m.cell(i,j)],yw=sy[m.cell(i-1,j)],ye=sy[m.cell(i,j)];
    double axl=.5*(xs.left+xn.left),axr=.5*(xs.right+xn.right),dx=.5*(xs.diffusion+xn.diffusion);
    double ayl=.5*(yw.left+ye.left),ayr=.5*(yw.right+ye.right),dy=.5*(yw.diffusion+ye.diffusion);
    double by_w=staggered_edge(a,i-1,j,0,1,m,p),by_e=staggered_edge(a,i,j,0,0,m,p);
    double bx_s=staggered_edge(a,i,j-1,1,1,m,p),bx_n=staggered_edge(a,i,j,1,0,m,p);
    double vr_w=face_velocity_edge(sy,i-1,j,0,1,m,p),vr_e=face_velocity_edge(sy,i,j,0,0,m,p);
    double vt_s=face_velocity_edge(sx,i,j-1,1,1,m,p),vt_n=face_velocity_edge(sx,i,j,1,0,m,p);
    double e=-axl*vr_w*by_w-axr*vr_e*by_e+ayl*vt_s*bx_s+ayr*vt_n*bx_n+dx*(by_e-by_w)-dy*(bx_n-bx_s);
    double ep=p.reconstruction_cache&&p.midplane&&p.eta_corner_geometry?eta_from_geometry(p.eta_corner_geometry[k],p):.5*(eta_at(u,i-1,m.xf[i],m.yf[j],m,p)+eta_at(u,i,m.xf[i],m.yf[j],m,p));
    double res=ep*current_phi(a,i,j,m,p);e+=res;bool axis=p.lower2==AXIS&&j==NG;
    arhs[m.vertex(i,j)]=axis?0.:-m.edge_length(i,j,p)*e;
    if(resistive_arhs)resistive_arhs[m.vertex(i,j)]=axis?0.:-m.edge_length(i,j,p)*res;
    if(transport&&transport_weight!=0&&!axis) {
        double scale=-m.edge_length(i,j,p);
        double adv=scale*(-axl*vr_w*by_w-axr*vr_e*by_e+ayl*vt_s*bx_s+ayr*vt_n*bx_n);
        double jump=scale*(dx*(by_e-by_w)-dy*(bx_n-bx_s));
        double resistive=scale*res;
        auto& q=transport[m.vertex(i,j)];
        q.q[TV_ADV]+=transport_weight*adv;q.q[TV_JUMP]+=transport_weight*jump;q.q[TV_RES]+=transport_weight*resistive;
        q.q[TV_ROUND]+=transport_weight*(arhs[m.vertex(i,j)]-(adv+jump+resistive));
    }
}
HD inline void cell_rhs(int k,const State* u,const double* a,const State* f1,const State* f2,State* rhs,MeshView m,Parameters p,
                        const State* res1=nullptr,const State* res2=nullptr,const double* res_a=nullptr,double* ohmic_heat=nullptr,double* ohmic_physical=nullptr,
                        const double* velocity1=nullptr,const double* velocity2=nullptr) {
    int i=k%m.nx+NG,j=k/m.nx+NG,ii=i-NG,jj=j-NG;
    State left=f1[jj*(m.nx+1)+ii],right=f1[jj*(m.nx+1)+ii+1];
    State down=f2[jj*m.nx+ii],up=f2[(jj+1)*m.nx+ii],out;
    double al=m.area1(i,j,p),ar=m.area1(i+1,j,p),ad=m.area2(i,j,p),au=m.area2(i,j+1,p),vol=m.volume(i,j,p);
    for(int c=0;c<NV;++c)out[c]=-(ar*right[c]-al*left[c]+au*up[c]-ad*down[c])/vol;
    Primitive v=get_primitive(u,i,j,m,p);
    if(p.spherical) {
        double qr=m.qr[i],qt=qr*m.qt[j],pt=v[ENE]+.5*(sq(v[B1])+sq(v[B2])+sq(v[B3]));
        double trt=v[RHO]*v[M1]*v[M2]-v[B1]*v[B2],ttt=v[RHO]*sq(v[M2])+pt-sq(v[B2]),tpp=v[RHO]*sq(v[M3])+pt-sq(v[B3]);
        out[M1]+=qr*(ttt+tpp)-v[RHO]*p.gm/sq(m.xc[i]);
        out[M2]+=qt*tpp-qr*trt;
        // Toroidal induction geometry includes both ideal and resistive EMF flux.
        out[B3]+=qr*.5*(left[B3]+right[B3])+qt*.5*(down[B3]+up[B3]);
        // This provisional FV Btheta is also needed for CT energy correction.
        out[B2]+=qr*.5*(left[B2]+right[B2]);
        // Conservative gravitational work: Eg + rho Phi(cell) has only face fluxes.
        double ph=potential(m.xc[i],p);
        out[ENE]-=(ar*right[RHO]*(potential(m.xf[i+1],p)-ph)-al*left[RHO]*(potential(m.xf[i],p)-ph))/vol;
    }
    double heating=0;
    if(res1&&res2&&res_a) {
        State rl=res1[jj*(m.nx+1)+ii],rr=res1[jj*(m.nx+1)+ii+1],rd=res2[jj*m.nx+ii],ru=res2[(jj+1)*m.nx+ii];
        double edot=-(ar*rr[ENE]-al*rl[ENE]+au*ru[ENE]-ad*rd[ENE])/vol;
        double b3dot=-(ar*rr[B3]-al*rl[B3]+au*ru[B3]-ad*rd[B3])/vol;
        if(p.spherical)b3dot+=m.qr[i]*.5*(rl[B3]+rr[B3])+m.qr[i]*m.qt[j]*.5*(rd[B3]+ru[B3]);
        double b1dot=.5*(field1(res_a,i,j,m,p)+field1(res_a,i+1,j,m,p));
        double b2dot=.5*(field2(res_a,i,j,m,p)+field2(res_a,i,j+1,m,p));
        if(p.ct_energy_correction) {
            b1dot=-(ar*rr[B1]-al*rl[B1]+au*ru[B1]-ad*rd[B1])/vol;
            b2dot=-(ar*rr[B2]-al*rl[B2]+au*ru[B2]-ad*rd[B2])/vol;
            if(p.spherical)b2dot+=m.qr[i]*.5*(rl[B2]+rr[B2]);
        }
        heating=edot-v[B1]*b1dot-v[B2]*b2dot-v[B3]*b3dot;
    } else {
        // Standalone operator probes can omit the staggered heat buffers.
        double ep=eta_cell_value(u,i,j,m,p),et=p.eta_model==0?ep:p.chi*ep;
        if(ep>0){Current c=cell_current(u,a,i,j,m,p);heating=et*(sq(c[0])+sq(c[1]))+ep*sq(c[2]);}
    }
    if(p.cooling)out[ENE]-=heating;
    // Immediate radiation of all Ohmic heat has exactly zero net entropy
    // source. Without cooling, use the actual discrete magnetic/energy RHS.
    if(p.entropy_switch&&!p.cooling)out[ENTR]+=(p.gamma-1)*pow(v[RHO],1-p.gamma)*heating;
    // Auxiliary gas energy is a volume density, not an entropy mass scalar.
    // Its conservative HLL transport preserves a constant-pressure density
    // contact; the pressure-work source supplies adiabatic compression.
    Primitive vl=get_primitive(u,i-1,j,m,p),vr=get_primitive(u,i+1,j,m,p),vd=get_primitive(u,i,j-1,m,p),vu=get_primitive(u,i,j+1,m,p);
    double divv=(vr[M1]-vl[M1])/(m.xc[i+1]-m.xc[i-1])
        +(vu[M2]-vd[M2])/((m.yc[j+1]-m.yc[j-1])*(p.spherical?m.xc[i]:1.));
    if(p.spherical)divv+=(2*v[M1]+m.cosine_center(j)/m.sine_center(j)*v[M2])/m.xc[i];
    if(velocity1&&velocity2)divv=(ar*velocity1[jj*(m.nx+1)+ii+1]-al*velocity1[jj*(m.nx+1)+ii]
        +au*velocity2[(jj+1)*m.nx+ii]-ad*velocity2[jj*m.nx+ii])/vol;
    out[EINT]-=v[ENE]*divv;
    if(!p.cooling)out[EINT]+=heating;
    if(ohmic_heat)ohmic_heat[k]=heating;
    if(ohmic_physical) {
        double ep=eta_cell_value(u,i,j,m,p),et=p.eta_model==0?ep:p.chi*ep;
        double value=0;if(ep>0){Current c=cell_current(u,a,i,j,m,p);value=et*(sq(c[0])+sq(c[1]))+ep*sq(c[2]);}
        ohmic_physical[k]=value;
    }
    rhs[m.cell(i,j)]=out;
}
