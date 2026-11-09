// Cached SF2014 geometry must not change the meaning of a constant/ideal
// runtime operator override. Compare eta, full face fluxes, and corner EMFs.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>

static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
int main(int argc,char** argv) {
    try {
        World w(argc,argv);require(w.size==1,"eta override test requires one process");
        auto p=setup_case_defaults(DISK);p.nx=32;p.ny=32;p.xhi=4;p.flat_reconstruction=1;
        Mesh mesh(p,w);Solver s(p,mesh.view,w);auto m=s.m;
        std::vector<double>a(s.a.n);std::vector<State>u(s.u.n);
        for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=.1*m.xf[i]*sq(sin(m.yf[j]));
        for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){Primitive v;v[RHO]=1;v[ENE]=1;v[B1]=.5*(field1(a.data(),i,j,m,p)+field1(a.data(),i+1,j,m,p));v[B2]=.5*(field2(a.data(),i,j,m,p)+field2(a.data(),i,j+1,m,p));v[B3]=.02*m.xc[i]*sin(m.yc[j]);u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);}
        s.a.upload(a);s.u.upload(u);s.exchange(s.u.p,s.a.p);s.refresh_midplane(s.u.p,s.a.p);
        double max_error=0;int checks=0;
        for(double eta:{0.,.037}) {
            Parameters cached=s.p;cached.eta_model=0;cached.eta=eta;
            Parameters direct=cached;direct.reconstruction_cache=0;
            for(int j=NG;j<NG+m.ny;++j)for(int i=NG;i<NG+m.nx;++i){
                double x=eta_cell_value(s.u.p,i,j,m,cached),y=eta_at(s.u.p,i,m.xc[i],m.yc[j],m,direct);
                require(x==eta&&y==eta,"Cached geometry ignored constant eta model");++checks;
            }
            std::vector<State>fc1(s.flux1.n),fd1(fc1.size()),fc2(s.flux2.n),fd2(fc2.size());
            for(int dir=0;dir<2;++dir){auto& fc=dir?fc2:fc1;auto& fd=dir?fd2:fd1;
                for(size_t k=0;k<fc.size();++k){face_flux(static_cast<int>(k),dir,s.u.p,s.a.p,fc.data(),m,cached);face_flux(static_cast<int>(k),dir,s.u.p,s.a.p,fd.data(),m,direct);for(int c=0;c<NV;++c){max_error=std::max(max_error,std::abs(fc[k][c]-fd[k][c]));require(fc[k][c]==fd[k][c],"Cached override face flux differs from direct constant model");++checks;}}
            }
            std::vector<double>ac(s.a.n),ad(s.a.n);
            for(int k=0;k<(m.nx+1)*(m.ny+1);++k){corner_emf(k,s.u.p,s.a.p,ac.data(),m,cached);corner_emf(k,s.u.p,s.a.p,ad.data(),m,direct);int i=k%(m.nx+1)+NG,j=k/(m.nx+1)+NG;max_error=std::max(max_error,std::abs(ac[m.vertex(i,j)]-ad[m.vertex(i,j)]));require(ac[m.vertex(i,j)]==ad[m.vertex(i,j)],"Cached override corner EMF differs from direct constant model");++checks;}
            // A zero vertical profile must still yield the requested constant.
            EtaGeometry vacuum;vacuum.profile=0;require(eta_from_geometry(vacuum,cached)==eta,"Constant model inherited a geometric cutoff");
        }
        std::cout<<std::setprecision(17)<<"PASS eta_model override comparisons="<<checks<<" maximum_difference="<<max_error<<'\n';return 0;
    }catch(const std::exception&error){std::cerr<<"FAIL eta model override: "<<error.what()<<'\n';return 1;}
}
