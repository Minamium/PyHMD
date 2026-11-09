#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
int main(int argc,char** argv) {
    World world(argc,argv);double maximum=0;long comparisons=0;
    auto compare=[&](double actual,double expected) {
        if(!std::isfinite(actual)||!std::isfinite(expected))throw std::runtime_error("nonfinite eta cache comparison");
        double error=std::abs(actual-expected)/std::max(std::abs(expected),1e-300);
        maximum=std::max(maximum,error);++comparisons;
        if(error>1e-13)throw std::runtime_error("fixed eta geometry changed diffusivity");
    };
    for(int model:{1,2})for(int log_r:{0,1}) {
        auto p=setup_case_defaults(DISK);p.nx=64;p.ny=32;p.log_r=log_r;p.eta_model=model;
        if(!log_r)p.xlo=100; // Keep the three radial halo cells strictly positive.
        Mesh mesh(p,world);auto m=mesh.view;Solver solver(p,m,world);
        std::vector<State> states(solver.u.n);std::vector<double> potential_a(solver.a.n);
        for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)potential_a[m.vertex(i,j)]=initial_potential(m.xf[i],m.yf[j],p);
        for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i) {
            auto v=initial(m.xc[i],m.yc[j],p);
            v[B1]=.5*(field1(potential_a.data(),i,j,m,p)+field1(potential_a.data(),i+1,j,m,p));
            v[B2]=.5*(field2(potential_a.data(),i,j,m,p)+field2(potential_a.data(),i,j+1,m,p));
            states[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);
        }
        solver.u.upload(states);solver.a.upload(potential_a);
        std::vector<Midplane> mid(p.nx);
        for(int i=0;i<p.nx;++i){mid[i].va=.03+.01*sin(.3*i);mid[i].cs=.1+.02*cos(.2*i);mid[i].mu=.01+.005*sin(.1*i);}
        solver.mid_global.upload(mid);solver.p.midplane=solver.mid_global.p;
        Parameters reference=solver.p;reference.reconstruction_cache=0;Parameters optimized=solver.p;
        for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i)
            compare(eta_cell_value(solver.u.p,i,j,m,optimized),eta_at(solver.u.p,i,m.xc[i],m.yc[j],m,reference));
        for(int dir=0;dir<2;++dir) {
            int ni=dir==0?m.nx+1:m.nx,nj=dir==0?m.ny:m.ny+1;
            auto geometry=dir==0?optimized.eta_radial_geometry:optimized.eta_theta_geometry;
            std::vector<State> raw(ni*nj),cached(ni*nj);
            for(int k=0;k<ni*nj;++k) {
                int i=k%ni+NG,j=k/ni+NG,il=i-(dir==0);
                double radius=dir==0?m.xf[i]:m.xc[i],theta=dir==1?m.yf[j]:m.yc[j];
                compare(eta_from_geometry(geometry[k],optimized),.5*(eta_at(solver.u.p,il,radius,theta,m,reference)+eta_at(solver.u.p,i,radius,theta,m,reference)));
                face_flux(k,dir,solver.u.p,solver.a.p,raw.data(),m,reference);face_flux(k,dir,solver.u.p,solver.a.p,cached.data(),m,optimized);
                for(int q=0;q<NV;++q)compare(cached[k][q],raw[k][q]);
            }
        }
        std::vector<double> raw(solver.a.n,0),cached(solver.a.n,0);
        for(int k=0;k<(m.nx+1)*(m.ny+1);++k) {
            int i=k%(m.nx+1)+NG,j=k/(m.nx+1)+NG;
            compare(eta_from_geometry(optimized.eta_corner_geometry[k],optimized),eta_at(solver.u.p,i,m.xf[i],m.yf[j],m,reference));
            corner_emf(k,solver.u.p,solver.a.p,raw.data(),m,reference);corner_emf(k,solver.u.p,solver.a.p,cached.data(),m,optimized);
        }
        for(size_t k=0;k<raw.size();++k)compare(cached[k],raw[k]);
    }
    std::cout<<"eta/flux/EMF comparisons="<<comparisons<<" maximum relative difference="<<maximum<<'\n';
}
