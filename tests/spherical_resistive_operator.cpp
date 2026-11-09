// Manufactured spherical curl/curl tests independent of disk evolution.
// eta_P is an exactly constant frozen coefficient, eta_T=chi*eta_P.
// Analytic solutions independently check metric terms,
// anisotropy orientation, and poloidal/toroidal resistive signs.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
#include <iomanip>

static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
struct Error {double poloidal=0,toroidal=0;};
static Error run(const World& w,int n) {
    auto p=setup_case_defaults(UNIFORM);p.nx=p.ny=n;p.xlo=1;p.xhi=4;p.ylo=.9;p.yhi=1.5;
    p.spherical=p.log_r=1;p.gm=0;p.lower1=p.upper1=p.lower2=p.upper2=OUTFLOW;
    p.eta_model=1;p.chi=.5;p.alpha_m=1.3;p.h=.5;p.reconstruction_cache=1;p.flat_reconstruction=1;
    Mesh mesh(p,w);auto m=mesh.view;constexpr double eta=.03,b=.01;
    std::vector<Midplane> mid(n);
    for(auto& point:mid)point.va=eta/p.alpha_m;
    p.midplane=mid.data();
    // Inject a known constant coefficient into the cached-geometry interface.
    // This deliberately tests the resistive operator independently of the
    // SF2014 closure and its radial interpolation. Unit height/profile and a
    // constant lookup ensure eta is exactly the same at every cell/face/edge.
    EtaGeometry uniform;uniform.index=0;uniform.height=uniform.profile=1;
    std::vector<EtaGeometry>gc(m.tx*m.ty,uniform),gr((m.nx+1)*m.ny,uniform),gt(m.nx*(m.ny+1),uniform),gv((m.nx+1)*(m.ny+1),uniform);
    p.eta_cell_geometry=gc.data();p.eta_radial_geometry=gr.data();p.eta_theta_geometry=gt.data();p.eta_corner_geometry=gv.data();
    std::vector<double>a((m.tx+1)*(m.ty+1));std::vector<State>u(m.tx*m.ty);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=b*pow(m.xf[i],3)*sq(sin(m.yf[j]));
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){Primitive v;v[RHO]=1;v[ENE]=1;
        v[B1]=.5*(field1(a.data(),i,j,m,p)+field1(a.data(),i+1,j,m,p));
        v[B2]=.5*(field2(a.data(),i,j,m,p)+field2(a.data(),i,j+1,m,p));
        v[B3]=b*sq(m.xc[i])*sin(m.yc[j]);u[m.cell(i,j)]=conserved(v,m.radius(i,j,p),p);}
    std::vector<State>rf1((m.nx+1)*m.ny),rf2(m.nx*(m.ny+1)),flux(std::max(rf1.size(),rf2.size())),rhs(u.size());
    std::vector<double>arhs(a.size()),res_a(a.size());
    for(size_t k=0;k<rf1.size();++k)face_flux(static_cast<int>(k),0,u.data(),a.data(),flux.data(),m,p,rf1.data());
    for(size_t k=0;k<rf2.size();++k)face_flux(static_cast<int>(k),1,u.data(),a.data(),flux.data(),m,p,rf2.data());
    for(int k=0;k<(m.nx+1)*(m.ny+1);++k)corner_emf(k,u.data(),a.data(),arhs.data(),m,p,res_a.data());
    for(int k=0;k<m.nx*m.ny;++k)cell_rhs(k,u.data(),a.data(),rf1.data(),rf2.data(),rhs.data(),m,p);
    double errp=0,normp=0,errt=0,normt=0;
    for(int j=NG+3;j<NG+m.ny-3;++j)for(int i=NG+3;i<NG+m.nx-3;++i) {
        if(m.xc[i]<2||m.xc[i]>3.4)continue;
        // Psi=b*r^3*sin²theta: Psi_t=4 eta b r sin²theta.
        double expected_p=4*eta*b*m.xf[i]*sq(sin(m.yf[j]));
        double measured_p=res_a[m.vertex(i,j)];errp+=sq(measured_p-expected_p);normp+=sq(expected_p);
        // Bphi=b*r²*sin(theta): Bphi_t=4 chi eta b sin(theta).
        double expected_t=4*p.chi*eta*b*sin(m.yc[j]);
        double measured_t=rhs[m.cell(i,j)][B3];errt+=sq(measured_t-expected_t);normt+=sq(expected_t);
    }
    return {sqrt(errp/normp),sqrt(errt/normt)};
}
int main(int argc,char** argv) {
    try {
        World w(argc,argv);require(w.size==1,"Manufactured spherical operator is a CPU single-process test");
        auto a=run(w,32),b=run(w,64),c=run(w,128);
        std::cout<<std::setprecision(17)<<"spherical operator poloidal="<<a.poloidal<<','<<b.poloidal<<','<<c.poloidal
                 <<" toroidal="<<a.toroidal<<','<<b.toroidal<<','<<c.toroidal<<'\n';
        require(a.poloidal/b.poloidal>3&&b.poloidal/c.poloidal>3,"Poloidal resistive curl failed second-order refinement");
        require(a.toroidal/b.toroidal>3&&b.toroidal/c.toroidal>3,"Toroidal anisotropic resistive geometry failed second-order refinement");
        require(c.poloidal<.003&&c.toroidal<.003,"Manufactured spherical operator error exceeds tolerance");
        std::cout<<std::setprecision(17)<<"PASS spherical resistive poloidal errors="<<a.poloidal<<','<<b.poloidal<<','<<c.poloidal
                 <<" toroidal errors="<<a.toroidal<<','<<b.toroidal<<','<<c.toroidal<<" etaT/etaP=0.5\n";return 0;
    }catch(const std::exception&error){std::cerr<<"FAIL spherical resistive operator: "<<error.what()<<'\n';return 1;}
}
