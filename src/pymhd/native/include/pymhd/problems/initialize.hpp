#pragma once
#include "pymhd/runtime/solver.hpp"
#include "pymhd/geometry/mesh.hpp"
#include "pymhd/problems/paper_initial.hpp"
#include <memory>

struct HostMesh {
    std::vector<double> xf,yf,xc,yc,qr,qt;MeshView view;
    explicit HostMesh(const Mesh& mesh):xf(mesh.xf.download()),yf(mesh.yf.download()),xc(mesh.xc.download()),yc(mesh.yc.download()),qr(mesh.qr.download()),qt(mesh.qt.download()),view(mesh.view) {
        view.xf=xf.data();view.yf=yf.data();view.xc=xc.data();view.yc=yc.data();view.qr=qr.data();view.qt=qt.data();
        view.radial_geometry=nullptr;view.polar_geometry=nullptr;
    }
};
static void initialize(Solver& s,const HostMesh& h) {
    auto m=h.view;auto p=s.p;std::vector<State> states(s.u.n);std::vector<double> aa(s.a.n);
    std::unique_ptr<PaperInitial> disk;if(p.problem==DISK)disk=std::make_unique<PaperInitial>(p);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)aa[m.vertex(i,j)]=initial_potential(m.xf[i],m.yf[j],p);
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i) {
        Primitive v=disk?disk->at(m.xc[i],m.yc[j]):initial(m.xc[i],m.yc[j],p);
        if(p.spherical&&m.yc[j]<0){v[M2]*=-1;v[M3]*=-1;v[B2]*=-1;v[B3]*=-1;}
        v[B1]=.5*(field1(aa.data(),i,j,m,p)+field1(aa.data(),i+1,j,m,p));
        v[B2]=.5*(field2(aa.data(),i,j,m,p)+field2(aa.data(),i,j+1,m,p));
        if(p.noise_amplitude>0) {
            v[RHO]*=1+p.noise_amplitude*initial_noise(m.start+i-NG,j-NG,p.noise_seed);
            v[ENTR]=v[ENE]*pow(v[RHO],-p.gamma);v[EINT]=v[ENE]/(p.gamma-1);
        }
        states[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
        if(p.problem==OBLIQUE_ALFVEN||p.problem==LINEAR_WAVE) {
            states[m.cell(i,j)]=benchmark_average(m.xf[i],m.xf[i+1],m.yf[j],m.yf[j+1],0,p);
            states[m.cell(i,j)][B1]=v[B1];states[m.cell(i,j)][B2]=v[B2];
        }
    }
    s.u.upload(states);s.a.upload(aa);s.exchange(s.u.p,s.a.p);
}
