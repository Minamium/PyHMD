#pragma once
#include "pyhmd/geometry/fields.hpp"
#include <cstdint>

// Analytic test definitions derived from the equations, not another code's source.
// Felker & Stone 2018, sections 5.1--5.7; see docs/phase0.md for conventions.
HD inline bool benchmark_problem(int problem){return problem>=OBLIQUE_ALFVEN&&problem<=BLAST_LOW_BETA;}
struct WaveFrame {double nx,ny,k,speed;};
HD inline WaveFrame wave_frame(Parameters p) {
    double kx=2*PI/(p.xhi-p.xlo),ky=2*PI/(p.yhi-p.ylo),k=sqrt(kx*kx+ky*ky);
    double speed=p.problem==LINEAR_WAVE?(p.wave_mode==0?2.:(p.wave_mode==2?.5:1.)):1.;
    return {kx/k,ky/k,k,speed};
}
HD inline void linear_eigenvector(int mode,double& dr,double& vn,double& vt,double& vz,double& dp,double& bt,double& bz) {
    dr=vn=vt=vz=dp=bt=bz=0;
    if(mode==3){dr=1;return;}
    if(mode==1){double norm=1.5;vt=-.5/norm;vz=sqrt(2.)/norm;bt=-vt;bz=-vz;return;}
    double speed=mode==0?2.:.5,den=speed*speed-1;
    vn=1;dr=1/speed;dp=1/speed;
    bt=sqrt(2.)*speed/den;bz=.5*speed/den;vt=-sqrt(2.)/den;vz=-.5/den;
}
HD inline Primitive benchmark_primitive(double x,double y,double time,Parameters p) {
    Primitive v;v[RHO]=1;v[ENE]=1;
    if(p.problem==OBLIQUE_ALFVEN||p.problem==LINEAR_WAVE) {
        auto f=wave_frame(p);double phase=f.k*(f.nx*(x-p.xlo)+f.ny*(y-p.ylo)-f.speed*time);
        double s=p.wave_amplitude*sin(phase),c=p.wave_amplitude*cos(phase);
        double bn=1,bt=0,bz=0,vn=0,vt=0,vz=0;
        if(p.problem==OBLIQUE_ALFVEN){v[ENE]=.1;bt=s;bz=c;vt=-s;vz=-c;}
        else {
            double dr,dvn,dvt,dvz,dp,dbt,dbz;linear_eigenvector(p.wave_mode,dr,dvn,dvt,dvz,dp,dbt,dbz);
            v[RHO]+=s*dr;v[ENE]=.6+s*dp;vn=(p.wave_mode==3?1.:0.)+s*dvn;
            vt=s*dvt;vz=s*dvz;bt=sqrt(2.)+s*dbt;bz=.5+s*dbz;
        }
        v[M1]=f.nx*vn-f.ny*vt;v[M2]=f.ny*vn+f.nx*vt;v[M3]=vz;
        v[B1]=f.nx*bn-f.ny*bt;v[B2]=f.ny*bn+f.nx*bt;v[B3]=bz;
    } else if(p.problem==FIELD_LOOP) {v[M1]=2;v[M2]=1;}
    else if(p.problem==RJ2A) {
        bool left=x<.5*(p.xlo+p.xhi);double unit=1/sqrt(PI);
        v[RHO]=left?1.08:1;v[ENE]=left?.95:1;v[M1]=left?1.2:0;v[M2]=left?.01:0;v[M3]=left?.5:0;
        v[B1]=unit;v[B2]=(left?1.8:2)*unit;v[B3]=unit;
    } else if(p.problem==ORSZAG_TANG) {
        double xx=x-p.xlo,yy=y-p.ylo,unit=1/sqrt(4*PI);
        v[RHO]=25/(36*PI);v[ENE]=5/(12*PI);v[M1]=-sin(2*PI*yy);v[M2]=sin(2*PI*xx);
        v[B1]=-unit*sin(2*PI*yy);v[B2]=unit*sin(4*PI*xx);
    } else if(p.problem==BLAST_LOW_BETA) {
        v[ENE]=x*x+y*y<.01?10:.1;v[B1]=v[B2]=1/sqrt(2.);
    }
    return v;
}
HD inline double benchmark_potential(double x,double y,double time,Parameters p) {
    if(p.problem==OBLIQUE_ALFVEN||p.problem==LINEAR_WAVE) {
        auto f=wave_frame(p);double bt=0,coefficient=1;
        if(p.problem==LINEAR_WAVE){double a,b,c,d,e,g;linear_eigenvector(p.wave_mode,a,b,c,d,e,coefficient,g);bt=sqrt(2.);}
        double bx=f.nx-f.ny*bt,by=f.ny+f.nx*bt;
        double phase=f.k*(f.nx*(x-p.xlo)+f.ny*(y-p.ylo)-f.speed*time);
        return bx*y-by*x+p.wave_amplitude*coefficient/f.k*cos(phase);
    }
    if(p.problem==FIELD_LOOP)return p.wave_amplitude*fmax(.3-sqrt(x*x+y*y),0.);
    if(p.problem==RJ2A){double mid=.5*(p.xlo+p.xhi),unit=1/sqrt(PI);return unit*y-unit*(x<mid?1.8*x:1.8*mid+2*(x-mid));}
    if(p.problem==ORSZAG_TANG){double unit=1/sqrt(4*PI);return unit*(cos(2*PI*(y-p.ylo))/(2*PI)+cos(4*PI*(x-p.xlo))/(4*PI));}
    if(p.problem==BLAST_LOW_BETA)return (y-x)/sqrt(2.);
    return 0;
}
inline State benchmark_average(double xl,double xr,double yl,double yr,double time,Parameters p) {
    // Eight-point Gaussian volume integration; independent error reporting also
    // retains the CT-synchronized B representation and its t=0 discrepancy.
    constexpr double node[4]={.1834346424956498,.525532409916329,.7966664774136267,.9602898564975363};
    constexpr double weight[4]={.362683783378362,.3137066458778873,.2223810344533745,.1012285362903763};
    State result;
    for(int a=0;a<4;++a)for(int sa:{-1,1})for(int b=0;b<4;++b)for(int sb:{-1,1}) {
        auto v=benchmark_primitive(.5*(xl+xr)+.5*(xr-xl)*sa*node[a],.5*(yl+yr)+.5*(yr-yl)*sb*node[b],time,p);
        auto u=conserved(v,1,p);double w=.25*weight[a]*weight[b];
        for(int c=0;c<NV;++c)result[c]+=w*u[c];
    }
    return result;
}
inline double initial_noise(int i,int j,int seed) {
    // Fixed 64-bit mixing; no rank, thread, or process-local RNG state.
    std::uint64_t z=static_cast<std::uint64_t>(static_cast<std::int64_t>(i))*0x9e3779b97f4a7c15ULL;
    z^=static_cast<std::uint64_t>(static_cast<std::int64_t>(j))*0xbf58476d1ce4e5b9ULL;
    z^=static_cast<std::uint64_t>(seed)+0x94d049bb133111ebULL;
    z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;z^=z>>31;
    return 2*static_cast<double>(z>>11)*0x1.0p-53-1;
}
