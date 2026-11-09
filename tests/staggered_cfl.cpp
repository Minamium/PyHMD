// Divergence-free checkerboard Psi: opposite strong face fields cancel in
// BOTH cell-centered components. CFL must still see the CT face wave speeds.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
int main(int argc,char**argv) {
    try{World w(argc,argv);auto p=setup_case_defaults(UNIFORM);p.nx=p.ny=16;p.spherical=p.log_r=0;
        p.xlo=p.ylo=0;p.xhi=p.yhi=1;p.gm=p.eta=0;p.entropy_switch=0;
        Mesh mesh(p,w);auto m=mesh.view;std::vector<State>u(m.tx*m.ty);std::vector<double>a((m.tx+1)*(m.ty+1));
        for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)a[m.vertex(i,j)]=((i+j)%2?-.5:.5);
        Primitive v;v[RHO]=1;v[ENE]=1e-6;std::fill(u.begin(),u.end(),conserved(v,1,p));
        int i=NG+7,j=NG+7,k=7*m.nx+7;
        double div=(m.area1(i+1,j,p)*field1(a.data(),i+1,j,m,p)-m.area1(i,j,p)*field1(a.data(),i,j,m,p)
            +m.area2(i,j+1,p)*field2(a.data(),i,j+1,m,p)-m.area2(i,j,p)*field2(a.data(),i,j,m,p))/m.volume(i,j,p);
        if(div!=0)throw std::runtime_error("Checkerboard fixture is not discretely divergence-free");
        p.cfl_staggered=0;auto old=rate_at(k,u.data(),a.data(),m,p);p.cfl_staggered=1;auto fixed=rate_at(k,u.data(),a.data(),m,p);
        if(!(fixed.adv>1000*old.adv&&fixed.bad==0))throw std::runtime_error("CFL missed cancelling strong face fields");
        a[m.vertex(i,j)]=NAN;auto bad=rate_at(k,u.data(),a.data(),m,p);
        if(bad.bad!=1)throw std::runtime_error("NaN staggered field hidden by max operation");
        std::cout<<"PASS divergence-free checkerboard CT CFL: centered_rate="<<old.adv<<" face_envelope_rate="<<fixed.adv<<" ratio="<<fixed.adv/old.adv<<"; nonfinite Psi rejected\n";return 0;
    }catch(const std::exception&e){std::cerr<<"FAIL staggered CFL: "<<e.what()<<'\n';return 1;}
}
