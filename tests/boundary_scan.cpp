// CPU simulation of the exact affine prefix schedule used by the CUDA kernel.
// GPU execution remains a required validation on a CUDA system.
#include "pymhd/setup.hpp"
#include "pymhd/solver.hpp"
int main(int argc,char** argv) {
    World world(argc,argv);double maximum=0;
    for(int ny:{16,128,257})for(int project:{0,1}) {
        auto p=setup_case_defaults(DISK);p.nx=64;p.ny=ny;p.inner_flux_project=project;Mesh mesh(p,world);auto m=mesh.view;
        std::vector<double> original((m.tx+1)*(m.ty+1));
        for(int j=0;j<=m.ty;++j)for(int i=0;i<=m.tx;++i)
            original[m.vertex(i,j)]=initial_potential(m.xf[i],m.yf[j],p)*(1+.1*sin(3*m.yf[j]));
        auto scan=original;paper_inner_potential(0,original.data(),m,p);
        for(int i=project?NG:NG-1;i>=0;--i) {
            scan[m.vertex(i,NG)]=0;double carry=0;
            for(int start=0;start<m.ny;start+=128) {
                PaperBoundaryAffine maps[128];
                for(int t=0;t<128;++t){int j=NG+start+t+1;maps[t]=j<=NG+m.ny?paper_boundary_affine(scan.data(),i,j,m):PaperBoundaryAffine{1,0};}
                for(int offset=1;offset<128;offset*=2) {
                    PaperBoundaryAffine next[128];
                    for(int t=0;t<128;++t)next[t]=t>=offset?compose_paper_boundary(maps[t],maps[t-offset]):maps[t];
                    std::copy(next,next+128,maps);
                }
                double next_carry=0;
                for(int t=0;t<128 && start+t<m.ny;++t){int j=NG+start+t+1;double value=maps[t].b+maps[t].a*carry;scan[m.vertex(i,j)]=value;next_carry=value;}
                carry=next_carry;
            }
        }
        for(int i=0;i<=NG;++i)for(int j=NG;j<=NG+m.ny;++j) {
            double expected=original[m.vertex(i,j)],difference=std::abs(expected-scan[m.vertex(i,j)]);
            maximum=std::max(maximum,difference);
            if(!std::isfinite(scan[m.vertex(i,j)]) || difference>2e-13*std::max(1.,std::abs(expected)))throw std::runtime_error("parallel prefix changed paper boundary");
        }
        if(!project)for(int j=NG;j<=NG+m.ny;++j) {
            double expected=initial_potential(m.xf[NG],m.yf[j],p)*(1+.1*sin(3*m.yf[j]));
            if(original[m.vertex(NG,j)]!=expected)throw std::runtime_error("Ghost fill overwrote physical normal flux");
        }
    }
    std::cout<<"paper boundary scan maximum absolute difference="<<maximum<<'\n';
}
