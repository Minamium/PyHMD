// Both Cartesian 1D limits of face-Riemann HLL UCT, including supersonic flow.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
static void run(const World&w,int dir,double velocity,int flat,int cache) {
    auto p=setup_case_defaults(UNIFORM);p.nx=p.ny=24;p.spherical=p.log_r=0;p.xlo=p.ylo=0;p.xhi=p.yhi=1;
    p.gm=p.eta=0;p.entropy_switch=0;p.ct_emf_scheme=1;p.flat_reconstruction=flat;p.reconstruction_cache=cache;
    p.lower1=p.upper1=p.lower2=p.upper2=OUTFLOW;
    Mesh mesh(p,w);auto m=mesh.view;Solver s(p,m,w);std::vector<State>u(s.u.n);std::vector<double>a(s.a.n);
    for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i){double x=m.xf[i],y=m.yf[j];
        a[m.vertex(i,j)]=dir==0?.7*y+.05*cos(2*PI*x)/(2*PI):-.7*x+.05*sin(2*PI*y)/(2*PI);}
    for(int j=0;j<m.ty;++j)for(int i=0;i<m.tx;++i){double x=dir==0?m.xc[i]:m.yc[j];Primitive v;
        v[RHO]=1+.1*sin(2*PI*x);v[ENE]=1+.03*cos(2*PI*x);v[M1+dir]=velocity+.1*sin(2*PI*x);v[M1+1-dir]=.2*cos(2*PI*x);v[B3]=.03*sin(2*PI*x);
        v[B1]=.5*(field1(a.data(),i,j,m,p)+field1(a.data(),i+1,j,m,p));v[B2]=.5*(field2(a.data(),i,j,m,p)+field2(a.data(),i,j+1,m,p));
        u[m.cell(i,j)]=conserved(v,1,p);}
    s.u.upload(u);s.a.upload(a);s.evaluate(s.u.p,s.a.p);auto emf=s.arhs.download();auto f=dir==0?s.flux1.download():s.flux2.download();double error=0,scale=1e-30;
    for(int j=NG+3;j<NG+p.ny-3;++j)for(int i=NG+3;i<NG+p.nx-3;++i){int index=dir==0?(j-NG)*(p.nx+1)+i-NG:(j-NG)*p.nx+i-NG;
        double expected=dir==0?-f[index][B2]:f[index][B1],actual=-emf[m.vertex(i,j)];
        error=std::max(error,fabs(expected-actual));scale=std::max(scale,fabs(expected));}
    if(error/scale>2e-12)throw std::runtime_error("Face UCT lost the 1D HLL flux limit");
    std::cout<<"PASS face UCT 1D dir="<<dir<<" v="<<velocity<<" flat="<<flat<<" cache="<<cache<<" relative_error="<<error/scale<<'\n';
}
int main(int argc,char**argv){try{World w(argc,argv);for(int d:{0,1})for(double v:{-8.,0.,8.})for(int f:{0,1})for(int c:{0,1})run(w,d,v,f,c);return 0;}catch(const std::exception&e){std::cerr<<"FAIL face UCT: "<<e.what()<<'\n';return 1;}}
